"""
.. module:: theory.slhaDecomposer
   :synopsis: SLHA-based SMS decomposition
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

import sys
import os
import copy
import time
import element
import topology
from tools.physicsUnits import addunit, rmvunit
import pyslha
import crossSection
from branch import Branch, decayBranches
import logging

logger = logging.getLogger(__name__)

 
def decompose(slhafile, sigcut=0.1, doCompress=False, doInvisible=False, 
              minmassgap=-1, useXSecs=None):
    """
    Do SLHA-based decomposition.
    
    :param slhafile: file with mass spectrum and branching ratios and
    optionally with cross-sections
    :param Xsec: optionally a dictionary with cross-sections for pair
    production, by default reading the cross sections from the SLHA file.
    :param XsecsInfo: information about the cross-sections (sqrts, order and
    label). Only relevant for Xsec=None (reading from slha file). If defined as
    input or in crossSection.XSectionInfo restricts the cross-sections values
    in the SLHA file to the ones in XsecsInfo. If not defined, it will be
    generated from the SLHA file and stored in crossSection.XSectionInfo. Only
    generated if cross-sections are read from SLHA file and not previously
    created
    :param sigcut: minimum sigma*BR to be generated, by default sigcut = 0.1 fb
    :param doCompress: turn mass compressed topologies on/off
    :param doInvisible: turn invisibly compressed topologies on/off
    :param minmassgap: maximum value for considering two R-odd particles
    degenerate (only revelant for doCompress=True)        
    :returns: a TopologyList.
     
    """
    t1 = time.time()
    
    if doCompress and rmvunit(minmassgap,'GeV') == -1: 
        logger.error("Please set minmassgap.")
        return False

    if type(sigcut) == type(1.):
        sigcut = addunit(sigcut, 'fb')

    # get cross-section from file
    xSectionList = crossSection.getXsecFromSLHAFile(slhafile, useXSecs)
    # get BRs and masses from file
    brDic, massDic = getDictionariesFromSLHA(slhafile)
    
    # Get maximum cross-sections (weights) for single particles (irrespective
    # of sqrtS):
    maxWeight = {}
    for pid in xSectionList.getPIDs():
        maxWeight[pid] = xSectionList.getXsecsFor(pid).getMaxXsec()
           
    # Create 1-particle branches with all possible mothers
    branchList = []
    for pid in maxWeight:
        branchList.append(Branch())
        branchList[-1].momID = pid
        branchList[-1].daughterID = pid
        branchList[-1].masses = [massDic[pid]]
        branchList[-1].maxWeight = maxWeight[pid]

    # Generate final branches (after all R-odd particles have decayed)
    finalBranchList = decayBranches(branchList, brDic, massDic, sigcut)
    
    smsTopList = topology.TopologyList()    
    # Combine pairs of branches into elements according to production
    # cross-section list:    
    for pids in xSectionList.getPIDpairs():
        for branch1 in finalBranchList:
            for branch2 in finalBranchList:
                if branch1.momID == pids[0] and branch2.momID == pids[1]:
                    finalBR = branch1.maxWeight*branch2.maxWeight/(
                            maxWeight[pids[0]]*maxWeight[pids[1]])
                    if type(finalBR) == type(addunit(1.,'fb')):
                        finalBR = finalBR.asNumber()                 
                    weightList = xSectionList.getXsecsFor(pids)*finalBR
                    
                    newElement = element.Element([branch1, branch2])
                    newElement.weight = weightList
                    allElements = [newElement]
                    # Do compression:
                    if doCompress or doInvisible:
                        allElements += newElement.compressElement(doCompress,
                                                                  doInvisible,
                                                                  minmassgap)
                    
                    for el in allElements:
                        if el.weight.getMaxXsec() < sigcut:
                            # skip elements with xsec below sigcut 
                            continue                          
                        top = topology.Topology(el)            
                        smsTopList.addList([top])                       
    logger.debug("slhaDecomposer done in " + str(time.time()-t1) + " s.")

    return smsTopList        


def getDictionariesFromSLHA(slhafile):
    """
    Read a SLHA file and get the mass and BR dictionaries.
    
    """
    workdir = os.getcwd() 
    pyslhadir = workdir + "/pyslha-1.4.3"
    sys.path.append(pyslhadir)
    res = pyslha.readSLHAFile(slhafile)

    # Get mass and branching ratios for all particles:
    brDic = {}
    for pid in res.decays.keys():
        brs = copy.deepcopy(res.decays[pid].decays)
        brsConj = copy.deepcopy(brs)
        for br in brsConj:
            br.ids = [-x for x in br.ids]
        brDic[pid] = brs
        brDic[-pid] = brsConj
    # Get mass list for all particles
    massDic = {}
    for pid in res.decays.keys():
        if pid and res.decays[pid].mass != None:
            massDic[pid] = addunit(abs(res.decays[pid].mass), 'GeV')
            massDic[-pid] = addunit(abs(res.decays[pid].mass), 'GeV')

    return brDic, massDic
