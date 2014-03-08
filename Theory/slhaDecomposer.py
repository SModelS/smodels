"""
.. module:: Theory.slhaDecomposer
   :synopsis: SLHA-based SMS decomposition
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

import sys, os, copy
import element, topology
from Tools.PhysicsUnits import addunit, rmvunit
import pyslha2 as pyslha
import crossSection
from branch import Branch, decayBranches
import logging

logger = logging.getLogger(__name__)
 
def decompose(slhafile,sigcut=0.1,DoCompress=False,DoInvisible=False,minmassgap=-1,UseXSecs=None):
    """Do SLHA-based decomposition.
    
    .. todo:: currently using pyslha2 because we need this hack to handle SLHA
    files with XSECTION blocks.

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
    :param DoCompress: turn mass compressed topologies on/off
    :param DoInvisible: turn invisibly compressed topologies on/off
    :param minmassgap: maximum value for considering two R-odd particles
    degenerate (only revelant for DoCompress=True)        
    :returns: a TopologyList.
     
    """
    if DoCompress and rmvunit(minmassgap,'GeV') == -1: 
        print "[SLHAdecomposer] Please, set minmassgap"
        return False

    if type(sigcut) == type(1.):
        sigcut = addunit(sigcut, 'fb')

    # get cross-section from file
    XSectionList = crossSection.getXsecFromSLHAFile(slhafile,UseXSecs)
    # get BRs and masses from file
    BRdic, Massdic = getDictionariesFromSLHA(slhafile)
    
    # Get maximum cross-sections (weights) for single particles (irrespective of sqrtS):
    maxWeight = {}
    for pid in XSectionList.getPIDs():
        maxWeight[pid] = XSectionList.getXsecsFor(pid).getMaxXsec()
           
    # Create 1-particle branches with all possible mothers
    branchList = []
    for pid in maxWeight:
        branchList.append(Branch())
        branchList[-1].momID = pid
        branchList[-1].daughterID = pid
        branchList[-1].masses = [Massdic[pid]]
        branchList[-1].maxWeight = maxWeight[pid]

    # Generate final branches (after all R-odd particles have decayed)
    finalBranchList = decayBranches(branchList,BRdic,Massdic,sigcut)
    
    SMSTopList = topology.TopologyList()    
    # Combine pairs of branches into elements according to production cross-section list:    
    for pids in XSectionList.getPIDpairs():
        for branch1 in finalBranchList:
            for branch2 in finalBranchList:
                if branch1.momID == pids[0] and branch2.momID == pids[1]:
                    FinalBR = branch1.maxWeight*branch2.maxWeight/(maxWeight[pids[0]]*maxWeight[pids[1]])
                    if type(FinalBR) == type(addunit(1.,'fb')):
                        FinalBR = FinalBR.asNumber()                 
                    weightList = XSectionList.getXsecsFor(pids)*FinalBR
                    
                    newElement = element.Element([branch1,branch2])
                    newElement.weight = weightList
                    # Do compression:
                    if DoCompress or DoInvisible:
                        compElements = newElement.compressElement(DoCompress,DoInvisible,minmassgap)
                    allElements = [newElement] + compElements
                    
                    for el in allElements:
                        if el.weight.getMaxXsec() < sigcut:
                            continue    # skip elements with xsec below sigcut                       
                        Top = topology.Topology(el)            
                        SMSTopList.addList([Top])                       

    return SMSTopList        


def getDictionariesFromSLHA(slhafile):
    """
    Read a SLHA file and get the mass and BR dictionaries.
    
    """
    workdir = os.getcwd() 
    pyslhadir = workdir + "/pyslha-1.4.3"
    sys.path.append(pyslhadir)
    res = pyslha.readSLHAFile(slhafile)

    # Get mass and branching ratios for all particles:
    BRdic = {}
    for pid in res[1].keys():
        brs = copy.deepcopy(res[1][pid].decays)
        brs_conj = copy.deepcopy(brs)
        for br in brs_conj:
            br.ids = [-x for x in br.ids]
        BRdic[pid] = brs
        BRdic[-pid] = brs_conj
    # Get mass list for all particles
    Massdic = {}
    for pid in res[1].keys():
        if pid and res[1][pid].mass != None:
            Massdic[pid] = addunit(abs(res[1][pid].mass),'GeV')
            Massdic[-pid] = addunit(abs(res[1][pid].mass),'GeV')

    return BRdic,Massdic
