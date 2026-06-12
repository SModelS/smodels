#!/usr/bin/env python3

"""
.. module:: Decomposer
   :synopsis: Decomposition of an input model into SMS Topologies (TheorySMS objects).

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import itertools
import time
from typing import Dict, Iterable, List, Optional, Set, Tuple, Union
from smodels.decomposition.theorySMS import TheorySMS
from smodels.decomposition.topologyDict import TopologyDict
from smodels.base.particleNode import ParticleNode
from smodels.base.physicsUnits import fb, GeV
from smodels.base.model import Model
from smodels.decomposition.exceptions import SModelSDecompositionError as SModelSError
from smodels.base.smodelsLogging import logger
from itertools import product
from collections import namedtuple
from smodels.base.particle import Particle, MultiParticle
from smodels.base.physicsUnits import UnitXSec,UnitEnergy
import numpy as np

# Auxiliary tuples for a lightweight representation of decays, subtrees and cross-sections during the decomposition process.
subtreeTuple = namedtuple('subtreeTuple', ['particleIDs', 'edges', 'decayBRs', 'canonName'])
decayTupleObj = namedtuple('decayTuple', ['mom', 'daughters', 'br'])
xsecTupleObj = namedtuple('xsecTuple', ['primaryMotherIDs', 'maxWeight', 'xsecList'])


def lightweight_sortTrees(subtreeList: Iterable[subtreeTuple], 
                          particleOrderDict: Dict[int, int],) -> List[subtreeTuple]:
    """
    Sort a list of subtree tuples first by canonName, then by order of the particles appearing in it.
    It assumes the sub-subtrees are arleady sorted by the same criteria.
    """

    sorted_list = sorted(subtreeList, 
                         key=lambda s: (s.canonName,tuple(particleOrderDict.get(pid, 0) 
                                                          for pid in s.particleIDs)))
    
    return sorted_list

def lightweight_sortParticleIDs( particleIDs: List[int], 
                                particleOrderDict: Dict[int, int]) -> List[int]:
    """
    Sort a list of particle IDs based on their order in the model.
    """

    sorted_list = sorted(particleIDs, 
                         key=lambda pid: particleOrderDict.get(pid, 0))
    
    return sorted_list

def get_particle_order_dict(model: Model) -> Dict[int, int]:
    """
    Get a dictionary mapping particle hash to an integer representing the order of the particle in the model.
    This is used for sorting leaves and trees.
    """

    sorted_particles = sorted(model.SMparticles + model.BSMparticles)
    particle_order_dict = {}
    for i, p in enumerate(sorted_particles):
        particle_order_dict[hash(p)] = i
    return particle_order_dict

def get_lightweight_decays(model: Model, 
                           particleOrderDict: Dict[int, int]) -> Dict[int, List[decayTupleObj]]:
    """
    Build a lightweight decay representation for all particles in the model, keyed by particle hash.
    The daughters in a given decay are sorted by their order from particleOrderDict to ensure consistent ordering when building subtrees.
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
                daughters_id = [hash(d) for d in decay.daughters]
                sorted_daughters = lightweight_sortParticleIDs(daughters_id, particleOrderDict)
                lightweight_decays.append(
                    decayTupleObj(
                        mom=pid,
                        daughters=sorted_daughters,
                        br=decay.br,
                    )
                )
        decaysDict[pid] = sorted(lightweight_decays, key=lambda d: d.br, reverse=True)
    
    return decaysDict

def get_lightweight_xsecs(model: Model, sigmacutFB: float) -> List[xsecTupleObj]:
    """
    Build a lightweight cross-section representation for all particle pairs in the model.
    The primary mothers in a given production channel are sorted by their order from particleOrderDict to ensure consistent ordering when building subtrees.
    """


    xSectionList = model.xsections
    xSectionList.removeLowerOrder()
    xSectionList.sort()
    
    xsecTupleList = []
    for pdgs in xSectionList.getPIDpairs():
        xsecList = xSectionList.getXsecsFor(pdgs)
        maxWeight = xsecList.getMaxXsec().asNumber(fb)
        if maxWeight < sigmacutFB:
            continue
        primaryMotherIDs = [hash(model.getParticle(pdg=pdg)) for pdg in pdgs]
        xsecTupleList.append(
            xsecTupleObj(primaryMotherIDs=primaryMotherIDs, maxWeight=maxWeight, xsecList=xsecList)
        )
    
    return xsecTupleList

def get_lightweight_canonName(sorted_subtrees: List[subtreeTuple]) -> int:
    """
    Get a canonical name for a subtree based on the canonNames of its daughter subtrees.
    The canonName is constructed as '1' + concatenation of sorted daughter canonNames + '0'.
    The subtrees are assumed to be already sorted by canonName and particle ordering, so that the same physical subtree will always have the same canonName regardless of the order in which the daughters were combined.
    """

    cName = '1'+"".join(f"{subtree.canonName}" for subtree in sorted_subtrees) + '0'
    return int(cName)

def build_subtree_cacheFor(particleID: int, decayDict: Dict[int, List[decayTupleObj]], 
                        particleOrderDict: Dict[int, int],  
                        memo: Optional[Dict[int, Tuple[subtreeTuple, ...]]] = None,
                        visiting: Optional[Set[int]] = None, minBR: float = 0.0) -> Dict[int, Tuple[subtreeTuple,...]]:
    """
    Build memoized subtree tuples for all descendants of particleID.
    """

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
                memo = build_subtree_cacheFor(daughter_id, decayDict, particleOrderDict, memo, visiting, minBR)
                daughter_subtree = memo[daughter_id]
                daughter_choices.append(daughter_subtree)

            for combo in product(*daughter_choices):
                all_BRs = decay.br*np.prod([daughter_subtree.decayBRs for daughter_subtree in combo])
                if all_BRs < minBR:
                    continue
                combo = lightweight_sortTrees(combo,particleOrderDict)  # Sort daughter subtrees by canonName and particle ordering
                cName = get_lightweight_canonName(combo)
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
    subtrees = lightweight_sortTrees(subtrees, particleOrderDict)
    visiting.remove(particleID)
    memo[particleID] = tuple(subtrees)
    return memo

def simplify_bsm_particles(model: Model) -> Dict[int, Union[MultiParticle, Particle]]:
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

def decompose(model: Model, sigmacut: Union[float,int,UnitXSec] = 0*fb, 
                 massCompress: bool = True, invisibleCompress: bool = True, 
                 minmassgap: UnitEnergy = 0.0*GeV, minmassgapISR: UnitEnergy = 0.0*GeV) -> TopologyDict:
     
    t0 = time.time()
    t1= time.time()

    if massCompress and minmassgap.asNumber(GeV) < 0.:
        logger.error("Asked for compression without specifying minmassgap. Please set minmassgap.")
        raise SModelSError()

    if isinstance(sigmacut, UnitXSec):
        sigmacutFB = sigmacut.asNumber(fb)  # sigmacut in fb (faster comparison)
    else:
        sigmacutFB = sigmacut


    # Define particle ordering for building sorted subtrees and topologies.
    particleOrderDict = get_particle_order_dict(model)

    # Lightweight decay representation keyed by particle hash.
    decaysDict = get_lightweight_decays(model,particleOrderDict)

    # Get BSM dict where equal BSM particles have been merged 
    # (to be used when building topologies)
    particleDict = simplify_bsm_particles(model)
    
    # Get lightweight cross-section representation above sigmacutFB
    xsecTupleList = get_lightweight_xsecs(model, sigmacutFB)
    # Sort by maxWeight (maximum cross-section)
    xsecTupleList.sort(key=lambda x: x.maxWeight, reverse=True)  

    logger.debug(f"{len(xsecTupleList)} production cross-sections obtained in {time.time() - t1:.2f} s.")
    t1 = time.time()


    smsTopDict = TopologyDict()
    if not xsecTupleList:
        return smsTopDict

    # Define a minimum BR for all subtrees
    maxXsec = max(x.maxWeight for x in xsecTupleList)
    minBR = sigmacutFB / maxXsec if maxXsec > 0.0 else 0.0


    # Build subtree cache for all primary mothers appearing in production channeprintls.
    cache = {}
    for xsecTuple in sorted(xsecTupleList, key=lambda x: x.maxWeight, reverse=True):
        # Define a minimum BR for building the subtrees.
        # Since the cross-sections are ordered, minBR increases as we go down the list, allowing more aggressive pruning of the subtree cache for subtrees which
        # are only needed for production channels with smaller cross-sections.
        minBR = sigmacutFB/xsecTuple.maxWeight
        for primaryMotherID in xsecTuple.primaryMotherIDs:
            if primaryMotherID in cache:
                continue
            cache = build_subtree_cacheFor(primaryMotherID, decaysDict, particleOrderDict, 
                                        memo=cache, minBR=minBR)
            
    
    pv = ParticleNode(model.getParticle(label='PV'))
    cache = {}
    nCascadeTrees = 0
    for xsecTuple in xsecTupleList:
        minBR = sigmacutFB/xsecTuple.maxWeight
        # Build subtree cache for primary mothers of the current production channel
        # (if not already built for a previous channel with larger cross-section)
        for primaryMotherID in xsecTuple.primaryMotherIDs:
            cache = build_subtree_cacheFor(primaryMotherID, decaysDict, particleOrderDict, 
                                            memo=cache, minBR=minBR)
        
        weight = xsecTuple.maxWeight
        all_subtrees = [cache[pid] for pid in xsecTuple.primaryMotherIDs]
        
        for primary_subtrees in itertools.product(*all_subtrees):
            totalBR = 1.0
            for subtree in primary_subtrees:
                totalBR *= subtree.decayBRs
            if weight*totalBR < sigmacutFB:
                continue

            # Sort the primary subtrees by canonName and particle ordering to ensure the 
            # tree is sorted
            primary_subtrees_sorted = lightweight_sortTrees(primary_subtrees, particleOrderDict)
            # Reorder the mother IDs in the same order as the sorted subtrees
            primary_mothers_sorted = [particleDict[s.particleIDs[0]] for s in primary_subtrees_sorted]
                        
            smsDecayed = TheorySMS()
            smsDecayed.maxWeight = xsecTuple.maxWeight
            smsDecayed.prodXSec = xsecTuple.xsecList
            smsDecayed._canonName = get_lightweight_canonName(primary_subtrees_sorted)
            pvIndex = smsDecayed.add_node(pv)
            primaryMothers = [ParticleNode(mother) for mother in primary_mothers_sorted]
            motherIndices = smsDecayed.add_nodes_from(primaryMothers)
            smsDecayed.add_edges_from(product([pvIndex], motherIndices))
            for idaughter, subtree in enumerate(primary_subtrees_sorted):
                old2newIndexMapping = {0: motherIndices[idaughter]}
                for nodeIndex,daughter_id in enumerate(subtree.particleIDs):
                    # Skip subtree root: the primary mother is already present in smsDecayed.
                    if nodeIndex == 0:
                        continue
                    daughter = particleDict[daughter_id]
                    # print('mom=',primary_mothers_sorted[idaughter], 'daughter=', daughter)
                    node = ParticleNode(daughter)
                    newIndex = smsDecayed.add_node(node)
                    old2newIndexMapping[nodeIndex] = newIndex
                for edgeA, edgeB in subtree.edges:
                    smsDecayed.add_edge(old2newIndexMapping[edgeA], old2newIndexMapping[edgeB])

            smsDecayed.decayBRs = totalBR
            smsDecayed.maxWeight = weight * totalBR 
            smsDecayed.weightList = smsDecayed.prodXSec*smsDecayed.decayBRs
            smsDecayed._sorted = True
            smsDecayed.ancestors = [smsDecayed]  # Set ancestors (before compression)
            smsTopDict.addSMS(smsDecayed)
            nCascadeTrees += 1


    logger.debug(f"{nCascadeTrees} cascade topologies trees generated and added to TopoDict in {time.time() - t1:.2f} s.")
    t1 = time.time()
    
    if massCompress or invisibleCompress:
        smsTopDict.compress(massCompress, invisibleCompress, minmassgap, minmassgapISR)
    
    logger.debug(f"Compression done in {time.time() - t1:.2f} s.")
    t1 = time.time()

    # Sort the topology dictionary according to the canonical names
    smsTopDict.sort()
    # Set the SMS IDs
    smsTopDict.setSMSIds()


    logger.debug(f"decomposer done in {time.time() - t0:.2f} s.")
    

    return smsTopDict