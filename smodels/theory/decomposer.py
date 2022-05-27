#!/usr/bin/env python3

"""
.. module:: Decomposer
   :synopsis: Decomposition of SLHA events and creation of TopologyLists.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>

"""

import time
from smodels.theory.element import Element
from smodels.theory.topology import TopologyDict
from smodels.theory.tree import Tree, ParticleNode
from smodels.tools.physicsUnits import fb, GeV
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger


def decompose(model, sigmacut=0 * fb, massCompress=True, invisibleCompress=True,
              minmassgap=0 * GeV):
    """
    Perform decomposition using the information stored in model.

    :param sigmacut: minimum sigma*BR to be generated, by default sigcut = 0.1 fb
    :param massCompress: turn mass compression on/off
    :param invisibleCompress: turn invisible compression on/off
    :param minmassgap: maximum value (in GeV) for considering two R-odd particles
                       degenerate (only revelant for massCompress=True )
    :returns: list of topologies (TopologyList object)

    """
    t1 = time.time()

    xSectionList = model.xsections
    if massCompress and minmassgap / GeV < 0.:
        logger.error("Asked for compression without specifying minmassgap. Please set minmassgap.")
        raise SModelSError()

    if isinstance(sigmacut, (float, int)):
        sigmacut = float(sigmacut) * fb

    xSectionList.removeLowerOrder()
    # Order xsections by highest xsec value to improve performance
    xSectionList.sort()

    # Generate all primary nodes (e.g. PV > X+Y)
    # and assign the nodeWeight as the maximum cross-section
    productionTrees = []
    for pid in xSectionList.getPIDpairs():
        weight = xSectionList.getXsecsFor(pid)
        if weight < sigmacut:
            continue
        pv = ParticleNode(model.getParticlesWith(label='PV')[0], 0, nodeWeight=weight)
        pv.xsection = xSectionList.getXsecsFor(pid)
        primaryMothers = [ParticleNode(model.getParticlesWith(pdg=pdg)[0], i + 1)
                          for i, pdg in enumerate(pid)]
        productionTrees.append(Tree({pv: primaryMothers}))

    # Sort production trees
    productionTrees = sorted(productionTrees,
                             key=lambda t: t.getTreeWeight().getMaxXsec(),
                             reverse=True)

    # For each production tree, produce all allowed cascade decays (above sigmacut):
    allTrees = []
    for tree in productionTrees:
        allTrees += cascadeDecay(tree, sigmacut=sigmacut)

    # Create elements for each tree and combine equal elements
    smsTopDict = TopologyDict()

    for tree in allTrees:
        newElement = Element(tree)
        newElement.weight = tree.getTreeWeight()
        smsTopDict.addElement(newElement)

    if massCompress or invisibleCompress:
        smsTopDict.compressElements(massCompress, invisibleCompress, minmassgap)
    smsTopDict._setElementIds()

    logger.debug("decomposer done in %.2f s." % (time.time() - t1))

    return smsTopDict


def getDecayTrees(mother):
    """
    Generates a simple list of trees with all the decay channels
    for the mother. In each tree the mother appears as the root
    and each of its decays as daughters.
    The  mother node weight is set to the respective decay branching ratio.
    (The node numbering for the root/mother node is kept equal,
    while the numbering of the daughters is automatically assigned to
    avoid overlap with any previously created nodes, so the
    decay tree can be directly merged to any other tree.)


    :param mother: Mother for which the decay trees will be generated (ParticleNode object)

    :return: List of simple trees representing the decay channels of daughter sorted in reverse
             order according to the BR (largest BR first)
    """

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
            decayTrees.append(Tree({mom: []}))

    decays = sorted(decays, key=lambda dec: dec.br, reverse=True)

    # Loop over decays of the daughter
    for decay in decays:
        if not decay.br:
            continue  # Skip decays with zero BRs
        daughters = []
        mom = mother.copy()
        mom.nodeWeight = decay.br
        for ptc in decay.daughters:
            ptcNode = ParticleNode(particle=ptc)
            daughters.append(ptcNode)

        decayTrees.append(Tree({mom: daughters}))

    return decayTrees


def addOneStepDecays(tree, sigmacut=None):
    """
    Given a tree, generates a list of new trees (Tree objects),
    where all the (unstable) nodes appearing at the end of the original tree
    have been decayed. Each entry in the list corresponds to a different combination
    of decays. If no decays were possible, return an empty list.

    :param tree: Tree (Tree object) for which to add the decays
    :param sigmacut: Cut on the tree weight (xsec*BR). Any tree with weights
                     smaller than sigmacut will be ignored.


    :return: List of trees with all possible 1-step decays added.
    """

    treeList = [tree]
    # Get all (current) final states which are the mothers
    # of the decays to be added:
    mothers = [n for n in tree.nodes() if tree.out_degree(n) == 0]
    for mom in mothers:
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

        # Get all decay trees for final state:
        decayTrees = getDecayTrees(mom)
        if not decayTrees:
            mom.isFinalState = True
            continue

        # Add all decay channels to all the trees
        newTrees = []
        for T in treeList:
            for decay in decayTrees:
                # The order below matters,
                # since we want to keep the mother from the decay tree (which holds the BR value)
                dec = decay.copy()
                newTree = dec.compose(T)
                if sigmacut is not None:
                    treeWeight = newTree.getTreeWeight()
                    if treeWeight is not None and treeWeight < sigmacut:
                        # Do not consider this decay or the next ones,
                        # since they are sorted accodring to BR and the subsequent
                        # decays will only contain smaller weights
                        break
                newTrees.append(newTree)

        if not newTrees:
            continue
        treeList = newTrees

    if len(treeList) == 1 and treeList[0] == tree:
        return []
    else:
        return treeList


def cascadeDecay(tree, sigmacut=None):
    """
    Given a tree, generates a list of new trees (Tree objects),
    where all the particles have cascade decayed to stable final states.

    :param tree: Tree (Tree object) for which to add the decays
    :param sigmacut: Cut on the tree weight (xsec*BR). Any tree with weights
                     smaller than sigmacut will be ignored.

    :return: List of trees with all possible decays added.
    """

    treeList = [tree]
    newTrees = True
    finalTrees = []
    while newTrees:
        newTrees = []
        for T in treeList:
            newT = addOneStepDecays(T, sigmacut)
            if not newT:
                finalNodes = [n for n in T.nodes if T.out_degree(n) == 0]
                # Make sure all the final states have decayed
                # (newT can be empty if there is no allowed decay above sigmacut)
                if any(fn.isFinalState is False for fn in finalNodes):
                    continue
                finalTrees.append(T)  # It was not possible to add any new decay to the tree
            else:
                newTrees += newT  # Add decayed trees to the next iteration

        if not newTrees:  # It was not possible to add any new decay
            break
        treeList = newTrees[:]

    return finalTrees
