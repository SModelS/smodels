"""
.. module:: theory.lheDecomposer
:synopsis: smodels-decomposes LHE events, creating TopologyLists 

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import LHEReader
import topology
import crossSection
import element
import pyslha2
import branch
import ParticleNames
from tools.PhysicsUnits import addunit
import logging
import copy

logger = logging.getLogger(__name__)


def decompose(lhefile, inputXsecs=None, nevts=None, doCompress=False, 
              doInvisible=False, minmassgap=None):
    """
    Do LHE-based decomposition. 

    :param lhefile: LHE file with e.g. pythia events
    :param inputXsecs: xSectionList object with cross-sections for the mothers
    appearing in the LHE file. If None, use information from file.
    :param nevts: (maximum) number of events used in the decomposition. If
    None, all events from file are processed.
    :param doCompress: mass compression option (True/False)
    :param doInvisible: invisible compression option (True/False)
    :param minmassgap: minimum mass gap for mass compression (only used if
    doCompress=True)
    :returns: a TopologyList object 
    
    """
    reader = LHEReader.LHEReader(lhefile, nevts)
    smsTopList = topology.TopologyList ( )
    # get cross-section from file (= event weight, assuming a common weight for
    # all events)
    if not inputXsecs:
        xSectionList = crossSection.getXsecFromLHEFile(lhefile, 
                                                       addEvents=False)
    else:
        xSectionList = inputXsecs

    # Loop over events and decompose 
    for event in reader:
        momPDG = tuple(sorted(event.getMom()))  # Get mother PDGs
        eventweight = xSectionList.getXsecsFor(momPDG)        
        # Get event element        
        newElement = elementFromEvent(event, eventweight)    
        # Do compression:        
        if doCompress or doInvisible:
            compElements = newElement.compressElement(doCompress, doInvisible,
                                                      minmassgap)
        allElements = [newElement] + compElements
        for el in allElements:
            top = topology.Topology(el)                      
            smsTopList.addList([top])                   

    return smsTopList


def elementFromEvent(event, weight):
    """
    Creates an element from a LHE event and the corresponding event weight.
    
    :param event: LHE event
    :param weight: event weight. Must be a XSectionList object (usually with a
    single entry)
    :returns: element
    
    """
    if not event.particles:
        logger.error('Empty event!')
        return None
          
    # Get simple BR and Mass dictionaries for each branch
    brDic, massdic = getDictionariesFromEvent(event)   
   
    # Creates branch list (list of mothers branches)
    branchList = []
    for particle in event.particles:
        # Means particle came from initial state (primary mother)
        if 1 in particle.moms:
            branchList.append(branch.Branch())
            branchList[-1].momID = particle.pdg
            branchList[-1].daughterID = particle.pdg
            branchList[-1].masses = [massdic[particle.pdg]]
            branchList[-1].maxWeight = weight.getMaxXsec()            
            
    # Generate final branches (after all R-odd particles have decayed)
    finalBranchList = branch.decayBranches(branchList, brDic, massdic,
                                           sigcut=addunit(0., 'fb'))     
    
    if len(finalBranchList) != 2:
        logger.error(str(len(finalBranchList))+" branches found in event. " +
                     "R-parity violation?")
        return False
    # Finally create element from event:
    newElement = element.Element(finalBranchList)
    newElement.weight = copy.deepcopy(weight)

    return newElement


def getDictionariesFromEvent(event):
    """
    Read an event and create simple mass and BR dictionaries for a single
    branch (mother).
    
    :param event: LHE event
    :returns: one BR dictionary for each branch and a common Mass dictionary
    
    """
    # Get mass and branching ratios for all particles:
    brDic = {}
    massdic = {}
    for particle in event.particles:
        if particle.pdg in ParticleNames.Reven:
            # Ignore R-even particles
            continue
        massdic[particle.pdg] = addunit(particle.mass, 'GeV')       
        # Empty BRs 
        brDic[particle.pdg] = [pyslha2.Decay(0., 0, [], particle.pdg)]
        
    for particle in event.particles:
        if particle.status == -1:
            # Ignore initial state particles
            continue
        if event.particles[particle.moms[0]].status == -1:
            # Ignore initial mothers
            continue
        if particle.moms[0] != particle.moms[1] and min(particle.moms) != 0:
            logger.error("More than one parent particle found!")
            return False        
        momPdg = event.particles[max(particle.moms)-1].pdg
        if momPdg in ParticleNames.Reven:
            # Ignore R-even decays
            continue
        # BR = 1 always for an event
        brDic[momPdg][0].br = 1.
        brDic[momPdg][0].nda += 1
        brDic[momPdg][0].ids.append(particle.pdg)
        
    return brDic, massdic
    