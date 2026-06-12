#!/usr/bin/env python3

"""
.. module:: Decomposer
   :synopsis: Decomposition of SLHA events and creation of TopologyLists.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>

"""

import time
from typing import Iterator
from smodels.decomposition.theorySMS import TheorySMS
from smodels.decomposition.topologyDict import TopologyDict
from smodels.base.particleNode import ParticleNode
from smodels.base.physicsUnits import fb, GeV
from smodels.decomposition.exceptions import SModelSDecompositionError as SModelSError
from smodels.base.smodelsLogging import logger
from itertools import product


def decompose(model, sigmacut=0 * fb, massCompress=True, invisibleCompress=True,
              minmassgap = 0*GeV, minmassgapISR = 0*GeV):
    """
    Perform decomposition using the information stored in model.

    :param sigmacut: minimum sigma*BR to be generated, by default sigcut = 0.1 fb
    :param massCompress: turn mass compression on/off
    :param invisibleCompress: turn invisible compression on/off
    :param minmassgap: maximum value (in GeV) for considering two BSM particles
                       degenerate (only revelant for massCompress=True )
    :param minmassgapISR: maximum value (in GeV) for mass compression leading to pure
                       ISR signature, i.e. PV > MET + MET + ... MET,
                       (only revelant for massCompress=True )
    :returns: list of topologies (TopologyList object)

    """
    t0 = time.time()
    t1= time.time()

    xSectionList = model.xsections
    if massCompress and minmassgap / GeV < 0.:
        logger.error("Asked for compression without specifying minmassgap. Please set minmassgap.")
        raise SModelSError()

    if isinstance(sigmacut, (float, int)):
        sigmacut = float(sigmacut) * fb
    sigmacutFB = sigmacut.asNumber(fb)  # sigmacut in fb (faster comparison)

    xSectionList.removeLowerOrder()
    # Order xsections by highest xsec value to improve performance
    xSectionList.sort()

    # Generate all primary nodes (e.g. PV > X+Y)
    # and assign the nodeWeight to the cross-section list
    productionSMS = []
    for pdgs in xSectionList.getPIDpairs():
        weight = xSectionList.getXsecsFor(pdgs)
        maxWeight = weight.getMaxXsec().asNumber(fb)
        if  maxWeight < sigmacutFB:
            continue
        pv = ParticleNode(model.getParticle(label='PV'))
        primaryMothers = [ParticleNode(model.getParticle(pdg=pdg)) for pdg in pdgs]
        newSMS = TheorySMS()
        newSMS.maxWeight = maxWeight
        newSMS.prodXSec = weight
        pvIndex = newSMS.add_node(pv)
        motherIndices = newSMS.add_nodes_from(primaryMothers)
        newSMS.add_edges_from(product([pvIndex],motherIndices))
        productionSMS.append(newSMS)

    # Sort production SMS by their maximum weights
    productionSMS = sorted(productionSMS,
                             key=lambda sms: sms.maxWeight,
                             reverse=True)
    logger.debug(f"{len(productionSMS)} primary production trees generated in {time.time() - t1:.2f} s.")
    t1 = time.time()

    # Create elements for each tree and combine equal elements
    smsTopDict = TopologyDict()
    nCascadeTrees = 0

    # For each production tree, produce all allowed cascade decays (above sigmacut):
    for sms in productionSMS:
        for decayedSMS in iterCascadeDecay(sms, sigmacutFB=sigmacutFB):
            decayedSMS.setGlobalProperties()  # Set global properties for each tree
            decayedSMS.ancestors = [decayedSMS]  # Set ancestors (before compression)
            smsTopDict.addSMS(decayedSMS)
            nCascadeTrees += 1
    logger.debug(f"{nCascadeTrees} cascade topologies trees generated and added to TopoDict in {time.time() - t1:.2f} s.")
    t1 = time.time()
    
    if massCompress or invisibleCompress:
        smsTopDict.compress(massCompress, invisibleCompress, 
                            minmassgap, minmassgapISR)
        
    logger.debug(f"Compression done in {time.time() - t1:.2f} s.")
    t1 = time.time()
    # Sort the topology dictionary according to the canonical names
    smsTopDict.sort()
    # Set the SMS IDs
    smsTopDict.setSMSIds()


    logger.debug(f"decomposer done in {time.time() - t0:.2f} s.")

    return smsTopDict


def getDecayNodes(mother):
    """
    Generates a simple list of trees with all the decay channels
    for the mother. In each tree the mother appears as the root
    and each of its decays as daughters.
    (The node numbering for the root/mother node is kept equal,
    while the numbering of the daughters is automatically assigned to
    avoid overlap with any previously created nodes, so the
    decay tree can be directly merged to any other tree.)


    :param mother: Mother for which the decay trees will be generated (ParticleNode object)

    :return: A list with simple tuples ((mom,daughters,BRs)) where
             the first entry is the new mother ParticleNode,
             the second is a list of daughter ParticleNode objects and
             the third the corresponding BR.
    """


    # If the trees were already computed, store them in the particle object
    if hasattr(mother.particle,'_decayTrees'):
        return mother.particle._decayTrees

    # Otherwise, compute them
    decayTrees = []

    # Sort decays:
    decays = []
    for decay in mother.decays:
        if decay is not None:
            decays.append(decay)
        else:
            # Include possibility of mother appearing as a final state
            mom = mother.copy()
            mom.isFinalState = True  # Forbids further node decays
            decayTrees.append((mom, [], 1.0))

    decays = sorted(decays, key=lambda dec: dec.br, reverse=True)
    # Loop over decays of the daughter
    for decay in decays:
        if not decay.br:
            continue  # Skip decays with zero BRs
        daughters = []
        mom = mother.copy()
        for ptc in decay.daughters:
            ptcNode = ParticleNode(particle=ptc)
            daughters.append(ptcNode)

        decayTrees.append((mom, daughters, decay.br))

    # Store the decays in the particle object
    mother.particle._decayTrees = decayTrees

    return decayTrees


def getMaxDecayBR(mother):
    """
    Return the maximum BR available for decaying ``mother`` once.

    For stable/final-state-like nodes (or nodes with no available decays), this
    returns 1.0 because no additional BR suppression is expected.
    """

    if mother.isFinalState:
        return 1.0
    if mother.particle.isStable():
        return 1.0
    if not hasattr(mother, 'decays'):
        return 1.0
    if not mother.decays:
        return 1.0

    # Cache per particle, since all nodes of the same species share decay data.
    if hasattr(mother.particle, '_maxDecayBR'):
        return mother.particle._maxDecayBR

    decayNodesList = getDecayNodes(mother)
    if not decayNodesList:
        maxBR = 1.0
    else:
        maxBR = max(decayNodes[2] for decayNodes in decayNodesList)

    mother.particle._maxDecayBR = maxBR
    return maxBR


def getMaxRemainingBRFactor(sms):
    """
    Compute an upper bound for additional BR suppression from all undecayed leaves.

    The returned factor is the product of the best one-step BR choices for each
    unresolved leaf and therefore provides a conservative upper bound for any
    fully decayed continuation of ``sms``.
    """

    factor = 1.0
    for nodeIndex in sms.nodeIndices:
        if sms.out_degree(nodeIndex) != 0:
            continue
        mom = sms.indexToNode(nodeIndex)
        factor *= getMaxDecayBR(mom)
        if factor == 0.0:
            return 0.0

    return factor


def iterCascadeDecay(tree, sigmacutFB=0.0) -> Iterator[TheorySMS]:
    """
    Yield final trees obtained by fully cascading all allowed decays.

    Uses a DFS stack that expands one unstable leaf at a time.  Peak memory
    scales as O(depth × max_branching) rather than O(max_branching^num_leaves),
    which avoids the combinatorial blow-up of the previous BFS approach when
    many leaves have multiple decay channels.

    :param tree: Tree (TheorySMS object) for which to add the decays
    :param sigmacutFB: Cut on the tree weight (xsec*BR) in fb. Any tree with
                     weights smaller than sigmacutFB will be ignored.

    :return: Iterator of fully-decayed TheorySMS trees.
    """

    stack = [tree]
    while stack:
        T = stack.pop()

        if T.maxWeight < sigmacutFB:
            continue

        # Find the first leaf that still needs to be decayed
        unstableIndex = None
        for nodeIndex in T.nodeIndices:
            if T.out_degree(nodeIndex) != 0:
                continue  # Internal node — skip
            mom = T.indexToNode(nodeIndex)
            if mom.isFinalState:
                continue
            if mom.particle.isStable():
                mom.isFinalState = True
                continue
            if not hasattr(mom, 'decays') or not mom.decays:
                mom.isFinalState = True
                continue
            unstableIndex = nodeIndex
            break

        if unstableIndex is None:
            # Every leaf is a final state — tree is fully decayed
            yield T
            continue

        # Branch-and-bound: even taking the best BR at every remaining leaf,
        # this subtree cannot reach sigmacutFB — prune it.
        if T.maxWeight * getMaxRemainingBRFactor(T) < sigmacutFB:
            continue

        mom = T.indexToNode(unstableIndex)
        decayNodesList = getDecayNodes(mom)
        if not decayNodesList:
            mom.isFinalState = True
            stack.append(T)
            continue

        for decayNodes in decayNodesList:
            br = decayNodes[2]
            if T.maxWeight * br < sigmacutFB:
                break  # decayNodesList is sorted by BR desc; remainder also fails
            newT = T.attachDecay(unstableIndex, decayNodes, br=br, copy=True)
            stack.append(newT)
