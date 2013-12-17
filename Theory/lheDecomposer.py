#!/usr/bin/env python

"""
.. module:: LHEDecomposer
    :synopsis: smodels-decomposes LHE events, creating TopologyLists 

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


import LHEReader, topology, crossSection, element, pyslha2, branch
import logging
logger = logging.getLogger(__name__)

def decompose(lhefile,inputXsecs=None,nevts=None,DoCompress=False,DoInvisible=False,minmassgap=None):
    """ Do LHE-based decomposition. 

    :param lhefile: LHE file with e.g. pythia events
    :param inputXsecs: XSectionList object with cross-sections for the mothers appearing in the LHE file.\
      If None, use information from file
    :param nevts: (maximum) number of events used in the decomposition. If None, all events from \
      file are processed.
    :param doCompress: mass compression option (True/False)
    :param doInvisible: invisible compression option (True/False)
    :param minmassgap: minimum mass gap for mass compression (only used if doCompress=True)
    :returns: a TopologyList object 
    """
  
    reader = LHEReader.LHEReader(lhefile,nevts)
    SMSTopList=topology.TopologyList ( )
#get cross-section from file
    if not inputXsecs:  XSectionList = crossSection.getXsecFromLHEFile(lhefile)
    else: XSectionList = inputXsecs 

#Loop over events and decompose 
    for Event in reader:    
        momPDG = tuple(Event.getMom())  # Get mother PDGs
        eventweight = XSectionList.getXsecsFor(momPDG)   
# Get event element
        newElement = elementFromEvent(Event,eventweight)
#Do compression:
        if DoCompress or DoInvisible: compElements = newElement.compressElement(DoCompress,DoInvisible,minmassgap)
        allElements = [newElement] + compElements
        for el in allElements:
            Top = topology.Topology(el)            
            SMSTopList.addList([Top])                       

    return SMSTopList


def elementFromEvent(Event,weight):
    """ Creates a element from a LHE event and the corresponding event weight
    :param Event: LHE event
    :param weight: event weight. Must be a XSectionList object (usually with a single entry)
    :returns: element
    """

    if not Event.particles:
        logger.error('[lheDecomposer]: Empty event!')
        return None
   
#Get simple BR and Mass dictionaries for each branch
    branchBRs, Massdic = getDictionariesFromEvent(Event)   
   
#Creates branch list (list of mothers)
    momBranches = {}
    for ip,particle in enumerate(Event.particles):
        if 1 in particle.moms:    # 1 means particle came from initial state
            newBranch = branch.Branch()            
            newBranch.momID = particle.pdg
            newBranch.daughterID = particle.pdg
            newBranch.masses = [particle.mass]
            momBranches[ip] = newBranch

#Generate final branches (after all R-odd particles have decayed)
    finalBranchList = []
    for momPos in momBranches:
        BRdic = branchBRs[momPos]
        finalBranchList.append(branch.decayBranches([momBranches[momPos]],BRdic,Massdic))
    
    if len(finalBranchList) != 2:
        logger.error("[lheDecomposer]: "+str(len(finalBranchList))+" branches found in event. R-parity violation?")
        return False
#Finally create element from Event:
    newElement = element.Element(finalBranchList)
    
    return newElement    

def getDictionariesFromEvent(Event):
    """ Read an event and create simple mass and BR dictionaries for a single branch (mother)
    :param Event: LHE event
    :returns: one BR dictionary for each branch and a common Mass dictionary
    """

#Get mass and branching ratios for all particles:
    BRdic = {}
    Massdic = {}
    branchBRs = {}
    for ip,particle in enumerate(Event.particles):                
        Massdic[particle.pdg] = particle.mass
        parent =particle
        while not 1 in parent.moms:   #Stop when primary mother has been found
            if parent.moms[0] != parent.moms[1]:
                logger.error("[lheDecomposer]: More than one parent particle found!")
                return False
            momPos = parent.moms[0]  #Mother position
            parent = Event.particles[momPos]
        if not momPos in branchBRs: branchBRs[momPos] = {}
        BRdic = branchBRs[momPos]   #Store BRs for this specific branch
        BRdic[particle.pdg] = []
        daughterIDs = []
        for particle2 in Event.particles:
            if ip in particle2.moms: daughterIDs.append(particle2.pdg)
#Remove possible repetition of entries in event list:            
        if len(daughterIDs) == 1 and  daughterIDs[0] == particle.pdg: continue
        BRdic[particle.pdg] = [pyslha2.Decay(1.,len(daughterIDs),daughterIDs)]  #Store single decay

    return branchBRs,Massdic
    
    