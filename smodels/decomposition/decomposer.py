#!/usr/bin/env python3

"""
.. module:: Decomposer
   :synopsis: Decomposition of SLHA events and creation of TopologyLists.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>

"""

import time
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
    t1 = time.time()

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

    # For each production tree, produce all allowed cascade decays (above sigmacut):
    allSMS = []
    for sms in productionSMS:
        allSMS += cascadeDecay(sms, sigmacutFB=sigmacutFB)

    # Create elements for each tree and combine equal elements
    smsTopDict = TopologyDict()

    for sms in allSMS:
        sms.ancestors = [sms]  # Set ancestors (before compression)
        # Sort SMS, compute canonical name and its total weight
        sms.setGlobalProperties()
        smsTopDict.addSMS(sms)

    if massCompress or invisibleCompress:
        smsTopDict.compress(massCompress, invisibleCompress, 
                            minmassgap, minmassgapISR)
    # Sort the topology dictionary according to the canonical names
    smsTopDict.sort()
    # Set the SMS IDs
    smsTopDict.setSMSIds()


    logger.debug(f"decomposer done in {time.time() - t1:.2f} s.")

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


def addOneStepDecays(sms, sigmacutFB=0.0):
    """
    Given a tree, generates a list of new trees (Tree objects),
    where all the (unstable) nodes appearing at the end of the original tree
    have been decayed. Each entry in the list corresponds to a different combination
    of decays. If no decays were possible, return an empty list.

    :param tree: Tree (Tree object) for which to add the decays
    :param sigmacutFB: Cut on the tree weight (xsec*BR) in fb. Any tree with weights
                     smaller than sigmacutFB will be ignored.


    :return: List of trees with all possible 1-step decays added.
    """

    smsList = [sms]
    # Get all (current) final states which are the mothers
    # of the decays to be added:
    motherIndices = [n for n in sms.nodeIndices if sms.out_degree(n) == 0]
    mothers = sms.indexToNode(motherIndices)
    for motherIndex,mom in zip(motherIndices,mothers):
        # Check if mom should decay:
        if mom.isFinalState:
            continue
        if mom.particle.isStable():
            mom.isFinalState = True
            continue  # Skip if particle is stable
        # Skip if particle has no decays
        if not hasattr(mom, 'decays'):
            mom.isFinalState = True
            continue
        if not mom.decays:
            mom.isFinalState = True
            continue

        # Get a list of decay nodes for final state (sorted by highest BR):
        decayNodesList = getDecayNodes(mom)
        if not decayNodesList:
            mom.isFinalState = True
            continue

        # Add all decay channels to all the trees
        newSMSList = []
        for T in smsList:
            tweight = T.maxWeight
            if tweight < sigmacutFB:
                continue
            for decayNodes in decayNodesList:
                br = decayNodes[2]
                if tweight*br < sigmacutFB:
                    break  # Since the decays are sorted, the next ones will also fall below sigmacut

                # Attach decay to original tree
                # (the mother node gets replaced by node from the decay dict)
                newSMS = T.attachDecay(motherIndex, decayNodes, br=br, copy=True)
                newSMSList.append(newSMS)

        if not newSMSList:
            continue
        smsList = sorted(newSMSList, key=lambda t: t.maxWeight,
                          reverse=True)

    if len(smsList) == 1 and smsList[0] is sms:
        return []
    else:
        return smsList


def cascadeDecay(tree, sigmacutFB=0.0):
    """
    Given a tree, generates a list of new trees (Tree objects),
    where all the particles have cascade decayed to stable final states.

    :param tree: Tree (Tree object) for which to add the decays
    :param sigmacutFB: Cut on the tree weight (xsec*BR) inf b. Any tree with weights
                     smaller than sigmacutFB will be ignored.

    :return: List of trees with all possible decays added.
    """

    treeList = [tree]
    finalTrees = []
    while treeList:
        newTrees = []
        for T in treeList:
            newT = addOneStepDecays(T, sigmacutFB)
            if not newT:
                finalNodes = [n for n in T.nodeIndices if T.out_degree(n) == 0]
                # Make sure all the final states have decayed
                # (newT can be empty if there is no allowed decay above sigmacutFB)
                if any(T.indexToNode(fn).isFinalState is False for fn in finalNodes):
                    continue
                finalTrees.append(T)  # It was not possible to add any new decay to the tree
            else:
                newTrees += newT  # Add decayed trees to the next iteration

        treeList = newTrees[:]

    return finalTrees
