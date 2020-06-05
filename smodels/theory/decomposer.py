#!/usr/bin/env python3

"""
.. module:: Decomposer
   :synopsis: Decomposition of SLHA events and creation of TopologyLists.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>

"""

import time
from smodels.theory import element, topology
from smodels.theory.branch import Branch, decayBranches
from smodels.tools.physicsUnits import fb, GeV
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger

def decompose(model, sigmacut= 0.1*fb, doCompress=True, doInvisible=True,
              minmassgap= 0*GeV):
    """
    Perform decomposition using the information stored in model.
    
    :param sigmacut: minimum sigma*BR to be generated, by default sigmacut = 0.1 fb
    :param doCompress: turn mass compression on/off
    :param doInvisible: turn invisible compression on/off
    :param minmassgap: maximum value (in GeV) for considering two R-odd particles
                       degenerate (only revelant for doCompress=True )
    :returns: list of topologies (TopologyList object)

    """
    t1 = time.time()
    
    xSectionList = model.xsections    
    pdgList = model.getValuesFor('pdg')

    if doCompress and minmassgap/GeV < 0.:
        logger.error("Asked for compression without specifying minmassgap. Please set minmassgap.")        
        raise SModelSError()

    if isinstance(sigmacut,(float,int)):
        sigmacut = sigmacut*fb
    sigmacut = sigmacut.asNumber(fb)

    xSectionList.removeLowerOrder()
    # Order xsections by PDGs to improve performance
    xSectionList.order()

    # Get maximum cross sections (weights) for single particles (irrespective
    # of sqrtS)
    maxWeight = {}
    for pid in xSectionList.getPIDs():
        maxWeight[pid] = xSectionList.getXsecsFor(pid).getMaxXsec().asNumber(fb)

    # Generate dictionary, where keys are the PIDs and values 
    # are the list of cross sections for the PID pair (for performance)
    xSectionListDict = {}    
    for pids in xSectionList.getPIDpairs():
        xSectionListDict[pids] = xSectionList.getXsecsFor(pids)

    # Create 1-particle branches with all possible mothers    
    branchList = []
    for pid in maxWeight:
        branchList.append(Branch())
        bsmParticle = model.getParticlesWith(pdg=pid)
        if not bsmParticle:
            raise SModelSError("Particle for pdg %i has not been defined.")
        if len(bsmParticle) != 1:
            raise SModelSError("Particle with pdg %i has multiple definitions.")
        branchList[-1].oddParticles = [bsmParticle[0]]
        if not pid in pdgList:
            logger.error("PDG %i has not been defined" %int(pid))
        branchList[-1].maxWeight = maxWeight[pid]        
        
    # Generate final branches (after all R-odd particles have decayed)
    finalBranchList = decayBranches(branchList, sigmacut)
    
    # Generate dictionary, where keys are the PIDs and values are the list of branches for the PID (for performance)
    branchListDict = {}
    for branch in finalBranchList:
        if branch.oddParticles[0].pdg in branchListDict:
            branchListDict[branch.oddParticles[0].pdg].append(branch)
        else:
            branchListDict[branch.oddParticles[0].pdg] = [branch]

    for pid in xSectionList.getPIDs():
        if not pid in branchListDict:
            branchListDict[pid] = []

    #Sort the branch lists by max weight to improve performance:
    for pid in branchListDict:
        branchListDict[pid] = sorted(branchListDict[pid], 
                                     key=lambda br: br.maxWeight, reverse=True)        

    smsTopList = topology.TopologyList()

    # Combine pairs of branches into elements according to production
    # cross section list
    for pids in xSectionList.getPIDpairs():
        weightList = xSectionListDict[pids]
        maxxsec = weightList.getMaxXsec().asNumber(fb)
        if maxxsec == 0.: ## protection
            continue
        minBR = sigmacut/maxxsec
        if minBR > 1.:
            continue
        for branch1 in branchListDict[pids[0]]:
            BR1 =  branch1.maxWeight/maxWeight[pids[0]]  #Branching ratio for first branch            
            if BR1 < minBR:
                break #Stop loop if BR1 is already too low
            for branch2 in branchListDict[pids[1]]:
                BR2 =  branch2.maxWeight/maxWeight[pids[1]]  #Branching ratio for second branch                
                if BR2 < minBR:
                    break #Stop loop if BR2 is already too low        
                                     
                finalBR = BR1*BR2
                if finalBR < minBR:
                    continue # Skip elements with xsec below sigmacut
                                       
                newElement = element.Element([branch1, branch2])
                newElement.weight = weightList*finalBR
                newElement.sortBranches()  #Make sure elements are sorted BEFORE adding them              
                smsTopList.addElement(newElement)                                                 
                                                    
    smsTopList.compressElements(doCompress, doInvisible, minmassgap)
    smsTopList._setElementIds()       
            
    logger.debug("decomposer done in %.2f s." % (time.time() -t1 ) )
  
    return smsTopList
