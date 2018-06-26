#!/usr/bin/env python

"""
.. module:: lheDecomposer
   :synopsis: Decomposition of LHE events and creation of TopologyLists.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.theory import lheReader, topology, crossSection, element
from smodels.theory import branch
from smodels.tools.physicsUnits import fb, GeV
from smodels.tools.smodelsLogging import logger
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
import pyslha
from smodels.particleDefinitions import SMpdgs,BSMpdgs
import copy


def decompose(lhefile, inputXsecs=None, nevts=None, doCompress=False,
              doInvisible=False, minmassgap=-1. * GeV ):
    """
    Perform LHE-based decomposition. 

    :param lhefile: LHE file with e.g. pythia events
    :param inputXsecs: xSectionList object with cross sections for the mothers
           appearing in the LHE file. If None, use information from file.
    :param nevts: (maximum) number of events used in the decomposition. If
                  None, all events from file are processed.
    :param doCompress: mass compression option (True/False)
    :param doInvisible: invisible compression option (True/False)
    :param minmassgap: minimum mass gap for mass compression (only used if
                       doCompress=True)
    :returns: list of topologies (TopologyList object) 
    
    """

    if doCompress and minmassgap < 0. * GeV:
        logger.error("Asked for compression without specifying minmassgap. Please set minmassgap.")
        raise SModelSError()

    reader = lheReader.LheReader(lhefile, nevts)
    smsTopList = topology.TopologyList()
    # Get cross section from file (= event weight, assuming a common weight for
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
        if not newElement:
            continue
        allElements = [newElement]
        # Perform compression
        if doCompress or doInvisible:
            allElements += newElement.compressElement(doCompress, doInvisible,
                                                      minmassgap)

        for el in allElements:
            el.sortBranches()            
            smsTopList.addElement(el)

    smsTopList._setElementIds()
    return smsTopList


def elementFromEvent(event, weight=None):
    """
    Creates an element from a LHE event and the corresponding event weight.
    
    :param event: LHE event
    :param weight: event weight. Must be a XSectionList object (usually with a
                   single entry) or None if not specified.
    :returns: element
    
    """
    if not event.particles:
        logger.error("Empty event")
        return None

    brDic, massDic = _getDictionariesFromEvent(event)

    # Create branch list
    finalBranchList = []
    for ip, particle in enumerate(event.particles):
        keys = SMpdgs+BSMpdgs
        if not particle.pdg in keys:
            logger.warning("Particle %i not defined in particle.py, events containing this particle will be ignored" %(particle.pdg))
            return None
        
        # Particle came from initial state (primary mother)
        if 1 in particle.moms:
            mombranch = branch.Branch()
            mombranch.PIDs = [[particle.pdg]]           
            if weight:
                mombranch.maxWeight = weight.getMaxXsec()
            else:
                mombranch.maxWeight = 0.*fb
            # Get simple BR and Mass dictionaries for the corresponding branch
            branchBR = brDic[ip]
            branchMass = massDic[ip]
            mombranch.masses = [branchMass[mombranch.PIDs[0][0]]]
            # Generate final branches (after all R-odd particles have decayed)
            finalBranchList += branch.decayBranches([mombranch], sigcut=0. * fb )
                                                     

    if len(finalBranchList) != 2:
        logger.error(str(len(finalBranchList)) + " branches found in event; "
                     "Possible R-parity violation")
        raise SModelSError()
    # Create element from event
    newElement = element.Element(finalBranchList)
    if weight:
        newElement.weight = copy.deepcopy(weight)

    return newElement


def _getDictionariesFromEvent(event):
    """
    Create mass and BR dictionaries for each branch in an event.
    
    :param event: LHE event
    :returns: BR and mass dictionaries for the branches in the event
    
    """

    particles = event.particles

    # Identify and label to which branch each particle belongs 
    #(the branch label is the position of the primary mother)
    branches = {}
    for ip, particle in enumerate(particles):
        if particle.status == -1:
            continue
        if particles[particle.moms[0]].status == -1:
            # If a primary mother, the branch index is its own position
            initMom = ip
        else:
            # If not a primary mother, check if particle has a single parent
            # (as it should)
            if particle.moms[0] != particle.moms[1] and \
                    min(particle.moms) != 0:
                logger.error("More than one parent particle found")
                raise SModelSError()
            initMom = max(particle.moms) - 1
            while particles[particles[initMom].moms[0]].status != -1:
                # Find primary mother (labels the branch)
                initMom = max(particles[initMom].moms) - 1
        branches[ip] = initMom

    # Get mass and BR dictionaries for all branches:
    massDic = {}
    brDic = {}
    for ibranch in branches.values():  #ibranch = position of primary mother
        massDic[ibranch] = {}
        brDic[ibranch] = {}
    for ip, particle in enumerate(particles):
        if particle.pdg in SMpdgs or particle.status == -1:
            # Ignore R-even particles and initial state particles
            continue
        ibranch = branches[ip]  # Get particle branch
        massDic[ibranch][particle.pdg] = round(particle.mass,1)* GeV
        # Create empty BRs
        brDic[ibranch][particle.pdg] = [pyslha.Decay(0., 0, [], particle.pdg)]

    # Get BRs from event
    for ip, particle in enumerate(particles):
        if particle.status == -1:
            # Ignore initial state particles
            continue
        if particles[particle.moms[0]].status == -1:
            # Ignore initial mothers
            continue
        ibranch = branches[ip]
        momPdg = particles[max(particle.moms) - 1].pdg
        if momPdg in SMpdgs:
            # Ignore R-even decays
            continue
        # BR = 1 always for an event
        brDic[ibranch][momPdg][0].br = 1.
        brDic[ibranch][momPdg][0].nda += 1
        brDic[ibranch][momPdg][0].ids.append(particle.pdg)

    return brDic, massDic
