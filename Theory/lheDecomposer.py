"""
.. module:: Theory.lheDecomposer
:synopsis: smodels-decomposes LHE events, creating TopologyLists 

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


import LHEReader, topology, crossSection, element, pyslha2, branch, ParticleNames
from Tools.PhysicsUnits import addunit
import logging, copy
logger = logging.getLogger(__name__)

def decompose(lhefile,inputXsecs=None,nevts=None,DoCompress=False,DoInvisible=False,minmassgap=None):
    """
    Do LHE-based decomposition. 

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
    # get cross-section from file (= event weight, assuming a common weight for all events)
    if not inputXsecs:
        XSectionList = crossSection.getXsecFromLHEFile(lhefile,addEvents=False)
    else:
        XSectionList = inputXsecs

    # Loop over events and decompose 
    for Event in reader:
        momPDG = tuple(sorted(Event.getMom()))  # Get mother PDGs
        eventweight = XSectionList.getXsecsFor(momPDG)        
        # Get event element        
        newElement = elementFromEvent(Event,eventweight)    
        # Do compression:        
        if DoCompress or DoInvisible:
            compElements = newElement.compressElement(DoCompress,DoInvisible,minmassgap)
        allElements = [newElement] + compElements
        for el in allElements:
            Top = topology.Topology(el)                      
            SMSTopList.addList([Top])                   

    return SMSTopList


def elementFromEvent(Event,weight):
    """
    Creates an element from a LHE event and the corresponding event weight.
    
    :param Event: LHE event
    :param weight: event weight. Must be a XSectionList object (usually with a single entry)
    :returns: element
    
    """
    if not Event.particles:
        logger.error('Empty event!')
        return None
          
    # Get simple BR and Mass dictionaries for each branch
    BRdic, Massdic = getDictionariesFromEvent(Event)   
   
    # Creates branch list (list of mothers branches)
    branchList = []
    for particle in Event.particles:
        if 1 in particle.moms:  # Means particle came from initial state (primary mother)
            branchList.append(branch.Branch())
            branchList[-1].momID = particle.pdg
            branchList[-1].daughterID = particle.pdg
            branchList[-1].masses = [Massdic[particle.pdg]]
            branchList[-1].maxWeight = weight.getMaxXsec()            
            
    # Generate final branches (after all R-odd particles have decayed)
    finalBranchList = branch.decayBranches(branchList,BRdic,Massdic,sigcut=addunit(0.,'fb'))     
    
    if len(finalBranchList) != 2:
        logger.error(str(len(finalBranchList))+" branches found in event. R-parity violation?")
        return False
    # Finally create element from Event:
    newElement = element.Element(finalBranchList)
    newElement.weight = copy.deepcopy(weight)

    return newElement


def getDictionariesFromEvent(Event):
    """
    Read an event and create simple mass and BR dictionaries for a single
    branch (mother).
    
    :param Event: LHE event
    :returns: one BR dictionary for each branch and a common Mass dictionary
    
    """
    # Get mass and branching ratios for all particles:
    BRdic = {}
    Massdic = {}
    for particle in Event.particles:
        if particle.pdg in ParticleNames.Reven:
            continue   #Ignore R-even particles        
        Massdic[particle.pdg] = addunit(particle.mass,'GeV')        
        BRdic[particle.pdg] = [pyslha2.Decay(0.,0,[],particle.pdg)]   #Empty  BRs
        
    for particle in Event.particles:
        if particle.status == -1:
            continue    # Ignore initial state particles
        if Event.particles[particle.moms[0]].status == -1:
            continue    # Ignore initial mothers
        if particle.moms[0] != particle.moms[1] and min(particle.moms) != 0:
            logger.error("More than one parent particle found!")
            return False        
        mom_pdg = Event.particles[max(particle.moms)-1].pdg
        if mom_pdg in ParticleNames.Reven:
            continue # Ignore R-even decays
        BRdic[mom_pdg][0].br = 1.   # BR = 1 always for an event
        BRdic[mom_pdg][0].nda += 1
        BRdic[mom_pdg][0].ids.append(particle.pdg)
        
    return BRdic,Massdic
    