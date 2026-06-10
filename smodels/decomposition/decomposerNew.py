#!/usr/bin/env python3

"""
.. module:: Decomposer
   :synopsis: Decomposition of SLHA events and creation of TopologyLists.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>

"""

import itertools
from typing import Iterator
from smodels.decomposition.theorySMS import TheorySMS
from smodels.decomposition.topologyDict import TopologyDict
from smodels.base.particleNode import ParticleNode
from smodels.base.physicsUnits import fb, GeV
from smodels.decomposition.exceptions import SModelSDecompositionError as SModelSError
from smodels.base.smodelsLogging import logger
from itertools import product
import numpy as np


def decomposeNew(model, sigmacut, massCompress, invisibleCompress, minmassgap, minmassgapISR):

    xSectionList = model.xsections
    sigmacutFB = sigmacut.asNumber(fb)  # sigmacut in fb (faster comparison)

    xSectionList.removeLowerOrder()
    # Order xsections by highest xsec value to improve performance
    xSectionList.sort()

    productionSMS = []
    smsTopDict = TopologyDict()
    for pdgs in xSectionList.getPIDpairs():
        weight = xSectionList.getXsecsFor(pdgs)
        maxWeight = weight.getMaxXsec().asNumber(fb)
        if maxWeight < sigmacutFB:
            continue
        pv = ParticleNode(model.getParticle(label='PV'))
        primaryMothers = [ParticleNode(model.getParticle(pdg=pdg)) for pdg in pdgs]
        newSMS = TheorySMS()
        newSMS.maxWeight = maxWeight
        newSMS.prodXSec = weight
        pvIndex = newSMS.add_node(pv)
        motherIndices = newSMS.add_nodes_from(primaryMothers)
        newSMS.add_edges_from(product([pvIndex], motherIndices))
        productionSMS.append(newSMS)

    maxXsec = max(sms.maxWeight for sms in productionSMS)
    minBR = sigmacutFB / maxXsec
    cache = {}
    for sms in productionSMS:
        for particleNode in sms.daughters(sms.rootIndex):
            if particleNode.particle is None:
                continue
            if particleNode.particle in cache:
                continue
            _, cache = build_subtree_cache(particleNode, memo=cache, minBR=minBR)

    for sms in productionSMS:
        all_subtrees = [
            cache.get(sms.indexToNode(daughterIndex).particle, [])
            for daughterIndex in sms.daughterIndices(sms.rootIndex)
        ]
        daughterIndices = list(sms.daughterIndices(sms.rootIndex))
        for primary_subtrees in itertools.product(*all_subtrees):
            totalBR = 1.0
            for subtree in primary_subtrees:
                totalBR *= subtree.decayBRs
            if sms.maxWeight * totalBR < sigmacutFB:
                continue

            smsDecayed = sms.copy()
            for idaughter, subtree in enumerate(primary_subtrees):
                old2newIndexMapping = {0: daughterIndices[idaughter]}
                for nodeIndex in subtree.nodeIndices:
                    # Skip subtree root: the primary mother is already present in smsDecayed.
                    if nodeIndex == subtree.rootIndex:
                        continue
                    node = subtree.indexToNode(nodeIndex)
                    newIndex = smsDecayed.add_node(node)
                    old2newIndexMapping[nodeIndex] = newIndex

                for edgeA, edgeB in subtree.edgeIndices:
                    smsDecayed.add_edge(old2newIndexMapping[edgeA], old2newIndexMapping[edgeB])

            smsDecayed.decayBRs = totalBR
            smsDecayed.maxWeight = sms.maxWeight * totalBR
            smsDecayed.setGlobalProperties()  # Set global properties for each tree
            smsDecayed.ancestors = [smsDecayed]  # Set ancestors (before compression)
            smsTopDict.addSMS(smsDecayed)



    if massCompress or invisibleCompress:
        smsTopDict.compress(massCompress, invisibleCompress, 
                            minmassgap, minmassgapISR)

    return smsTopDict

def build_subtree_cache(particleNode, memo=None, visiting=None, minBR=0.0):
    """Return all weighted TheorySMS subtrees rooted at `particle` and cache by particle."""

    if memo is None:
        memo = {}
    if visiting is None:
        visiting = set()
    if particleNode.particle in memo:
        return memo[particleNode.particle], memo

    if particleNode.particle in visiting:
        raise ValueError(
            f"Cycle detected at particle {particleNode.particle}; this DFS cache expects a DAG/tree-like decay graph."
        )

    logger.debug(f"Building subtree cache for particle {particleNode.particle}")
    visiting.add(particleNode.particle)
    subtrees = []
    if not hasattr(particleNode.particle, "decays") or not particleNode.particle.decays:
        decays = None
    else:
        decays = getDecayNodes(particleNode)
    
    if not decays:
        sms = TheorySMS()
        sms.add_node(particleNode)
        sms.decayBRs = 1.0
        subtrees.append(sms)
    else:
        for decay in decays:
            mom,daughters,br = decay
            if br < minBR:
                continue
            daughter_choices = []
            for daughter in daughters:
                daughter_subtrees, memo = build_subtree_cache(daughter, memo, visiting, minBR)
                daughter_choices.append(daughter_subtrees)
            
            for combo in product(*daughter_choices):
                sms = TheorySMS()
                root_index = sms.add_node(mom)

                all_weight = br*np.prod([daughter_sms.decayBRs for daughter_sms in combo])
                if all_weight < minBR:
                    continue
                
                for daughter_sms in combo:
                    index_map = {}
                    for node_index in daughter_sms.nodeIndices:
                        node = daughter_sms.indexToNode(node_index)
                        node_copy = node.copy() if hasattr(node, "copy") else node
                        new_index = sms.add_node(node_copy)
                        index_map[node_index] = new_index

                    for edge_a, edge_b in daughter_sms.edgeIndices:
                        sms.add_edge(index_map[edge_a], index_map[edge_b])

                    sms.add_edge(root_index, index_map[daughter_sms.rootIndex])
                    
                sms.decayBRs = all_weight
                subtrees.append(sms)

    visiting.remove(particleNode.particle)
    # Sort subtrees by decay BRs (descending) to improve performance of later pruning based on sigmacut.
    subtrees = sorted(subtrees, key=lambda sms: sms.decayBRs, reverse=True)
    memo[particleNode.particle] = tuple(subtrees)
    return memo[particleNode.particle], memo


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