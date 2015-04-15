#!/usr/bin/env python

"""
.. module:: theory.slhaDecomposer
   :synopsis: Decomposition of SLHA events and creation of TopologyLists.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import copy
import time
from smodels.theory import element, topology, crossSection
from smodels.theory.branch import Branch, decayBranches
from smodels.tools import modpyslha as pyslha
from smodels.tools.physicsUnits import fb, GeV
import smodels.particles
import logging
import sys

logger = logging.getLogger(__name__)

def decompose(slhafile, sigcut=.1 * fb, doCompress=False, doInvisible=False,
              minmassgap=-1.*GeV, useXSecs=None):
    """
    Perform SLHA-based decomposition.
    
    :param sigcut: minimum sigma*BR to be generated, by default sigcut = 0.1 fb
    :param doCompress: turn mass compression on/off
    :param doInvisible: turn invisible compression on/off
    :param minmassgap: maximum value (in GeV) for considering two R-odd particles
                       degenerate (only revelant for doCompress=True )
    :param useXSecs: optionally a dictionary with cross-sections for pair
                 production, by default reading the cross sections
                 from the SLHA file.
    :returns: list of topologies (TopologyList object)

    """
    t1 = time.time()

    if doCompress and minmassgap / GeV < 0.:
        logger.error("Asked for compression without specifying minmassgap. Please set minmassgap.")        
        sys.exit()

    if type(sigcut) == type(1.):
        sigcut = sigcut * fb

    from smodels.tools import modpyslha
    try:
        f=modpyslha.readSLHAFile ( slhafile )
    except modpyslha.ParseError,e:
        logger.error ( "The file %s cannot be parsed as an SLHA file: %s" % (slhafile, e) )
        sys.exit()

    # Get cross-section from file
    xSectionList = crossSection.getXsecFromSLHAFile(slhafile, useXSecs)
    # Get BRs and masses from file
    brDic, massDic = _getDictionariesFromSLHA(slhafile)
    # Only use the highest order cross-sections for each process
    xSectionList.removeLowerOrder()
    # Order xsections by PDGs to improve performance
    xSectionList.order()

    # Get maximum cross-sections (weights) for single particles (irrespective
    # of sqrtS)
    maxWeight = {}
    for pid in xSectionList.getPIDs():
        maxWeight[pid] = xSectionList.getXsecsFor(pid).getMaxXsec()    

    # Generate dictionary, where keys are the PIDs and values 
    # are the list of cross-sections for the PID pair (for performance)
    xSectionListDict = {}    
    for pids in xSectionList.getPIDpairs():
        xSectionListDict[pids] = xSectionList.getXsecsFor(pids)

    # Create 1-particle branches with all possible mothers
    branchList = []
    for pid in maxWeight:
        branchList.append(Branch())
        branchList[-1].momID = pid
        branchList[-1].daughterID = pid
        if not pid in massDic:
            logger.error ( "pid %d does not appear in masses dictionary %s in slhafile %s" % 
                    ( pid, massDic, slhafile ) )
        branchList[-1].masses = [massDic[pid]]
        branchList[-1].maxWeight = maxWeight[pid]

    # Generate final branches (after all R-odd particles have decayed)    
    finalBranchList = decayBranches(branchList, brDic, massDic, sigcut)
    # Generate dictionary, where keys are the PIDs and values are the list of branches for the PID (for performance)
    branchListDict = {}
    for branch in finalBranchList:
        if branch.momID in branchListDict:
            branchListDict[branch.momID].append(branch)
        else:
            branchListDict[branch.momID] = [branch]
    for pid in xSectionList.getPIDs():
        if not pid in branchListDict: branchListDict[pid] = []

    smsTopList = topology.TopologyList()
    # Combine pairs of branches into elements according to production
    # cross-section list
    for pids in xSectionList.getPIDpairs():
        weightList = xSectionListDict[pids]
        minBR = (sigcut/weightList.getMaxXsec()).asNumber()
        if minBR > 1.: continue
        for branch1 in branchListDict[pids[0]]:
            for branch2 in branchListDict[pids[1]]:
                finalBR = branch1.maxWeight * branch2.maxWeight / \
                            (maxWeight[pids[0]] * maxWeight[pids[1]])
                if type(finalBR) == type( 1. * fb):
                    finalBR = finalBR.asNumber()
                if finalBR < minBR: continue # Skip elements with xsec below sigcut

                newElement = element.Element([branch1, branch2])
                newElement.weight = weightList*finalBR
                allElements = [newElement]
                # Perform compression
                if doCompress or doInvisible:
                    allElements += newElement.compressElement(doCompress,
                                                                  doInvisible,
                                                                  minmassgap)

                for el in allElements:
                    top = topology.Topology(el)
                    smsTopList.addList([top])
    logger.debug("slhaDecomposer done in " + str(time.time() - t1) + " s.")
    return smsTopList


def _getDictionariesFromSLHA(slhafile):
    """
    Create mass and BR dictionaries from an SLHA file.
    Ignore decay blocks with R-parity violating or unknown decays

    """

    res = pyslha.readSLHAFile(slhafile)

    rOdd = smodels.particles.rOdd.keys()
    rEven = smodels.particles.rEven.keys()
    
    # Get mass and branching ratios for all particles
    brDic = {}
    for pid in res.decays.keys():
        if pid in rEven:
            logger.info("Ignoring %s decays",smodels.particles.rEven[pid])
            continue
        brs = []
        for decay in res.decays[pid].decays:
            nEven = nOdd = 0.
            for pidd in decay.ids:
                if pidd in rOdd: nOdd += 1
                elif pidd in rEven: nEven += 1
            if nOdd + nEven == len(decay.ids) and nOdd == 1:
                brs.append(decay)
            else:
                logger.info("Ignoring decay: %i -> [%s]",pid,decay.ids)

        brsConj = copy.deepcopy(brs)
        for br in brsConj:
            br.ids = [-x for x in br.ids]
        brDic[pid] = brs
        brDic[-pid] = brsConj
    # Get mass list for all particles
    massDic = dict(res.blocks['MASS'].items())
    for pid in massDic.keys()[:]:
        massDic[pid] *= GeV
        massDic[pid] = abs(massDic[pid])
        if not -pid in massDic: massDic[-pid] = massDic[pid]    
 
    return brDic, massDic

