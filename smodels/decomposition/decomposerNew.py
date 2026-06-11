#!/usr/bin/env python3

"""
.. module:: Decomposer
   :synopsis: Decomposition of SLHA events and creation of TopologyLists.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>

"""

import itertools
from typing import Iterator, List, Dict, Tuple, Union
from smodels.decomposition.theorySMS import TheorySMS
from smodels.decomposition.topologyDict import TopologyDict
from smodels.base.particleNode import ParticleNode
from smodels.base.physicsUnits import fb, GeV
from smodels.decomposition.exceptions import SModelSDecompositionError as SModelSError
from smodels.base.smodelsLogging import logger
from itertools import product
from collections import namedtuple
from smodels.base.particle import Particle, MultiParticle
import numpy as np

decaySubtreeTuple = namedtuple('decaySubtreeTuple', ['daughterIDs', 'br', 'maxBR'])
subtreeTuple = namedtuple('subtreeTuple', ['particleIDs', 'edges', 'decayBRs', 'canonName'])
decayTupleObj = namedtuple('decayTuple', ['mom', 'daughters', 'br'])
xsecTupleObj = namedtuple('xsecTuple', ['primaryMotherIDs', 'maxWeight', 'xsecList'])


def build_subtree_cache(particleID, decayDict, memo=None, visiting=None, minBR=0.0):
    """Build memoized subtree tuples for all descendants of particleID."""
    if memo is None:
        memo = {}
    if visiting is None:
        visiting = set()

    if particleID in memo:
        return memo

    if particleID in visiting:
        raise ValueError(f"Cycle detected at particle {particleID}; decay graph must be a DAG.")

    visiting.add(particleID)
    decays = decayDict.get(particleID, [])

    subtrees = []
    if not decays:
        parent_subtree = subtreeTuple(edges=[], particleIDs=[particleID,], decayBRs=1.0, canonName=10)
        subtrees.append(parent_subtree)
    else:
        for decay in decays:
            # Keep parity with decomposerNew.py behavior for direct decay pruning.
            if decay.br < minBR:
                continue

            child_ids = tuple(sorted(decay.daughters))
            daughter_choices = []
            for daughter_id in child_ids:
                memo = build_subtree_cache(daughter_id, decayDict, memo, visiting, minBR)
                daughter_subtree = memo[daughter_id]
                daughter_choices.append(daughter_subtree)

            for combo in product(*daughter_choices):
                all_BRs = decay.br*np.prod([daughter_subtree.decayBRs for daughter_subtree in combo])
                if all_BRs < minBR:
                    continue
                combo = sorted(combo, key=lambda s: s.canonName)  # Sort daughter subtrees by canonName for consistent ordering
                cName = '1'+"".join(f"{daughter_subtree.canonName}" for daughter_subtree in combo) + '0'
                parent_subtree = subtreeTuple(edges=[], particleIDs=[particleID,], decayBRs=all_BRs, canonName=int(cName))
                for daughter_subtree in combo:
                    index_map = {}
                    for idx,daughter_id in enumerate(daughter_subtree.particleIDs):
                        parent_subtree.particleIDs.append(daughter_id)
                        new_index = len(parent_subtree.particleIDs)-1
                        index_map[idx] = new_index                                                

                    for edge_a, edge_b in daughter_subtree.edges:
                        parent_subtree.edges.append((index_map[edge_a], index_map[edge_b]))

                    parent_subtree.edges.append((0, index_map[0]))

                subtrees.append(parent_subtree)
    subtrees.sort(key=lambda s: (s.canonName,s.particleIDs), reverse=False)
    visiting.remove(particleID)
    memo[particleID] = tuple(subtrees)
    return memo

def get_lightweight_decays(model) -> Dict[int, List[decayTupleObj]]:
    """
    Build a lightweight decay representation for all particles in the model, keyed by particle hash.
    """

    decaysDict = {}
    for p in model.BSMparticles + model.SMparticles:
        pid = hash(p)
        decays = getattr(p, 'decays', [])
        lightweight_decays = []
        for decay in decays:
            if decay is None:
                lightweight_decays.append(decayTupleObj(mom=pid, daughters=[], br=1.0))
            else:
                lightweight_decays.append(
                    decayTupleObj(
                        mom=pid,
                        daughters=[hash(d) for d in decay.daughters],
                        br=decay.br,
                    )
                )
        decaysDict[pid] = sorted(lightweight_decays, key=lambda d: d.br, reverse=True)
    
    return decaysDict

def get_lightweight_xsecs(model, sigmacutFB):


    xSectionList = model.xsections
    xSectionList.removeLowerOrder()
    xSectionList.sort()
    
    xsecTupleList = []
    for pdgs in xSectionList.getPIDpairs():
        xsecList = xSectionList.getXsecsFor(pdgs)
        maxWeight = xsecList.getMaxXsec().asNumber(fb)
        if maxWeight < sigmacutFB:
            continue
        primaryMotherIDs = tuple(hash(model.getParticle(pdg=pdg)) for pdg in pdgs)
        xsecTupleList.append(
            xsecTupleObj(primaryMotherIDs=primaryMotherIDs, maxWeight=maxWeight, xsecList=xsecList)
        )
    
    return xsecTupleList

def simplify_bsm_particles(model) -> Dict[int, Union[MultiParticle, Particle]]:
    """
    Simplify BSM particles by merging particles which can be considered as equal.
    These particles should be used to replaced the original particles in the SMS topologies
    and reduce the number of physically equivalent SMS generated during decomposition.
    """

    bsmList = []
    for bsm_particle in model.BSMparticles:
        if bsm_particle in bsmList:
            index = bsmList.index(bsm_particle)
            bsmList[index] = bsmList[index] + bsm_particle        
        else:
            bsmList.append(bsm_particle)

    # For the SM particles, directly used the defined particles/multiparticles without merging
    particleDict = {hash(p): p for p in model.SMparticles}
    for bsm_particle in bsmList:
        if type(bsm_particle) == Particle:
            particleDict[hash(bsm_particle)] = bsm_particle
        elif type(bsm_particle) == MultiParticle:
            for p in bsm_particle.particles:
                particleDict[hash(p)] = bsm_particle
        else:
            raise SModelSError(f"Unexpected particle type {type(bsm_particle)} in BSM particle list.")

    return particleDict

def decomposeNew(model, sigmacut, massCompress, invisibleCompress, minmassgap, minmassgapISR):
    
    
    # Lightweight decay representation keyed by particle hash.
    decaysDict = get_lightweight_decays(model)

    # Get BSM dict where equal BSM particles have been merged to be used when building topologies.
    particleDict = simplify_bsm_particles(model)
    
    
    sigmacutFB = sigmacut.asNumber(fb)
    xsecTupleList = get_lightweight_xsecs(model, sigmacutFB)

    smsTopDict = TopologyDict()
    if not xsecTupleList:
        return smsTopDict

    maxXsec = max(x.maxWeight for x in xsecTupleList)
    minBR = sigmacutFB / maxXsec if maxXsec > 0.0 else 0.0


    # Build subtree cache for all primary mothers appearing in production channeprintls.
    cache = {}
    for xsecTuple in xsecTupleList:
        for primaryMotherID in xsecTuple.primaryMotherIDs:
            if primaryMotherID in cache:
                continue
            cache = build_subtree_cache(primaryMotherID, decaysDict, memo=cache, minBR=minBR)

    pv = ParticleNode(model.getParticle(label='PV'))
    for xsecTuple in xsecTupleList:
        
        weight = xsecTuple.maxWeight
        primary_ids = xsecTuple.primaryMotherIDs
        all_subtrees = [cache[pid] for pid in primary_ids]
        
        for primary_subtrees in itertools.product(*all_subtrees):
            totalBR = 1.0
            for subtree in primary_subtrees:
                totalBR *= subtree.decayBRs
            if weight*totalBR < sigmacutFB:
                continue

            smsDecayed = TheorySMS()
            smsDecayed.maxWeight = xsecTuple.maxWeight
            smsDecayed.prodXSec = xsecTuple.xsecList
            pvIndex = smsDecayed.add_node(pv)
            primaryMothers = [ParticleNode(particleDict[momID]) for momID in primary_ids]
            motherIndices = smsDecayed.add_nodes_from(primaryMothers)
            smsDecayed.add_edges_from(product([pvIndex], motherIndices))
            for idaughter, subtree in enumerate(primary_subtrees):
                old2newIndexMapping = {0: motherIndices[idaughter]}
                for nodeIndex,daughter_id in enumerate(subtree.particleIDs):
                    # Skip subtree root: the primary mother is already present in smsDecayed.
                    if nodeIndex == 0:
                        continue
                    daughter = particleDict[daughter_id]
                    node = ParticleNode(daughter)
                    newIndex = smsDecayed.add_node(node)
                    old2newIndexMapping[nodeIndex] = newIndex
                for edgeA, edgeB in subtree.edges:
                    smsDecayed.add_edge(old2newIndexMapping[edgeA], old2newIndexMapping[edgeB])

            smsDecayed.decayBRs = totalBR
            smsDecayed.maxWeight = weight * totalBR 
            smsDecayed.weightList = smsDecayed.prodXSec*smsDecayed.decayBRs
            # Speed up here, since the  SMS are already sorted!
            smsDecayed.setGlobalProperties(sort=False)  # Set global properties for each tree
            smsDecayed.ancestors = [smsDecayed]  # Set ancestors (before compression)
            smsTopDict.addSMS(smsDecayed)


    if massCompress or invisibleCompress:
        smsTopDict.compress(massCompress, invisibleCompress, minmassgap, minmassgapISR)

    return smsTopDict