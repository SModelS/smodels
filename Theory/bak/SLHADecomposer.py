#!/usr/bin/env python

"""
.. module:: SLHADecomposer
        :synopsis: SLHA-based SMS decomposition
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

import sys, os, copy
import TopologyBuilder, SMSDataObjects, ParticleNames
from Tools.PhysicsUnits import addunit, rmvunit
import pyslha2 as pyslha
import CrossSection

        
def decompose(slhafile,sigcut=0.1,DoCompress=False,DoInvisible=False,minmassgap=-1):
    """Do SLHA-based decomposition.
            FIXME currently using pyslha2 because we need this hack to handle SLHA files with XSECTION blocks.

        :param slhafile: file with mass spectrum and branching ratios and optionally with cross-sections
        :param Xsec: optionally a dictionary with cross-sections for pair production, by default reading the cross sections from the SLHA file.
        :param XsecsInfo: information about the cross-sections (sqrts, order and label). Only relevant for Xsec=None (reading from slha file).\
             If defined as input or in CrossSection.XSectionInfo restricts the cross-sections values in the SLHA file to the ones in XsecsInfo. \
             If not defined, it will be generated from the SLHA file and stored in CrossSection.XSectionInfo.
Only generated if cross-sections are read from SLHA file and not previously created
        :param sigcut: minimum sigma*BR to be generated, by default sigcut = 0.1 fb
        :param DoCompress: turn mass compressed topologies on/off
        :param DoInvisible: turn invisibly compressed topologies on/off
        :param minmassgap: maximum value for considering two R-odd particles degenerate (only revelant for DoCompress=True)
        :returns: a TopologyList. 
    """

    if DoCompress and rmvunit(minmassgap,'GeV') == -1: 
        print "[SLHAdecomposer] Please, set minmassgap"
        return False

    if type(sigcut) == type(1.): sigcut = addunit(sigcut, 'fb')

#get cross-section from file
    XSectionList = CrossSection.getXsecFromFile(slhafile) 
#get BRs and masses from file
    BRdic, Massdic = getDictionariesFromSLHA(slhafile)

#Remove all cross-sections below sigcut:
    for xsec in XSectionList.XSections:
        if xsec.value < sigcut: XSectionList.delete(xsec)
#Get maximum cross-sections (weights) for single particles (irrespective of sqrtS):
    maxWeight = {}
    for pid in XSectionList.getPIDs(): maxWeight[pid] = XSectionList.getXsecsFor(pid).getMaxXsec()
           
#Construct 1-particle branches with all possible mothers
    branchList = []
    for pid in maxWeight:
        branchList.append(SMSDataObjects.Branch())
        branchList[-1].momID = pid
        branchList[-1].daughterID = pid
        branchList[-1].masses = [Massdic[pid]]
        branchList[-1].maxWeight = maxWeight[pid]

    allDecayed = False
#Now loop over branches and keep adding all decay possibilities till the end of the cascade decay
    while not allDecayed:        
        newBranchList = []
        allDecayed = True
        for branch in branchList:
            brs = BRdic[branch.daughterID]    #List of possible decays (BRs) for daughter in branch
            if len(brs) != 0: allDecayed = False  #False as long as there are unstable daughters
            newBranches = branch.addDecays(brs,Massdic)  #Add all possible decays to original branch
            newBranchList.extend(newBranches)

        branchList = []
        for newbranch in newBranchList:
            if newbranch.maxWeight > sigcut: branchList.append(copy.deepcopy(newbranch))  #Just keep the branches above sigcut


    SMSTopList = SMSDataObjects.TopologyList()
#Combine pairs of branches into Eelements according to production cross-section list:
    for pids in XSectionList.getPIDpairs():
        for branch1 in branchList:
            for branch2 in branchList:
                if branch1.momID == pids[0] and branch2.momID == pids[1]:
                    FinalBR = branch1.maxWeight*branch2.maxWeight/(maxWeight[pids[0]]*maxWeight[pids[1]])
                    if type(FinalBR) == type(addunit(1.,'fb')): FinalBR = FinalBR.asNumber()                 
                    weightList = XSectionList.getXsecsFor(pids)*FinalBR
                    if weightList.getMaxXsec() < sigcut: continue   #skip elements with xsec below sigcut

                    NewElement = SMSDataObjects.EElement()
                    NewElement.B = [copy.deepcopy(branch1),copy.deepcopy(branch2)]
                    NewElement.weight = weightList.getDictionary(groupBy="pids")[pids]         
                    EInfo = NewElement.getEinfo()
                    Top = SMSDataObjects.GTop()
                    Top.vertnumb = EInfo["vertnumb"]
                    Top.vertparts = EInfo["vertparts"]
                    Top.ElList.append(NewElement)
                    if not Top.checkConsistency(): continue
                    TopList = []
                    TopList.append(Top) 

                    # Do compression:                    
                    if DoCompress or DoInvisible:
                        FinalTopList = TopologyBuilder.compressTopology(TopList,DoCompress,DoInvisible,minmassgap)
                    else:
                        FinalTopList = TopList
                    # Add topologies to topology list:
                    SMSTopList.addList ( FinalTopList )

    return SMSTopList        


def getDictionariesFromSLHA(slhafile):
    """ Read a SLHA file and get the mass and BR dictionaries  """

    workdir = os.getcwd() 
    pyslhadir = workdir + "/pyslha-1.4.3"
    sys.path.append(pyslhadir)
    res = pyslha.readSLHAFile(slhafile)

#Get mass and branching ratios for all particles:
    BRdic = {}
    for pid in res[1].keys():
        brs = copy.deepcopy(res[1][pid].decays)
        brs_conj = copy.deepcopy(brs)
        for br in brs_conj:    br.ids = [-x for x in br.ids]
        BRdic[pid] = brs
        BRdic[-pid] = brs_conj
#Get mass list for all particles
    Massdic = {}
    for pid in res[1].keys():
        if pid and res[1][pid].mass != None:
            Massdic[pid] = addunit(abs(res[1][pid].mass),'GeV')
            Massdic[-pid] = addunit(abs(res[1][pid].mass),'GeV')

    return BRdic,Massdic

