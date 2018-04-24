#!/usr/bin/env python

"""
.. module:: slhaDecomposer
   :synopsis: Decomposition of SLHA events and creation of TopologyLists.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>

"""

import copy
import time
import sys
import pyslha
from smodels.theory import element, topology, crossSection
from smodels.theory.branch import Branch, decayBranches
from smodels.tools.physicsUnits import fb, GeV
from smodels.particleDefinitions import BSMList, BSMpdgs
from smodels.theory.particleNames import getObjectFromPdg
from smodels.theory.updateParticles import addPromptAndDisplaced
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger

def decompose(slhafile, sigcut=.1 * fb, doCompress=False, doInvisible=False,
              minmassgap=-1.*GeV, useXSecs=None):
    """
    Perform SLHA-based decomposition.
    
    :param sigcut: minimum sigma*BR to be generated, by default sigcut = 0.1 fb
    :param doCompress: turn mass compression on/off
    :param doInvisible: turn invisible compression on/off
    :param minmassgap: maximum value (in GeV) for considering two R-odd particles
                       degenerate (only revelant for doCompress=True )
    :param useXSecs: optionally a dictionary with cross sections for pair
                 production, by default reading the cross sections
                 from the SLHA file.
    :returns: list of topologies (TopologyList object)

    """
    t1 = time.time()

    if doCompress and minmassgap / GeV < 0.:
        logger.error("Asked for compression without specifying minmassgap. Please set minmassgap.")        
        raise SModelSError()

    if type(sigcut) == type(1.):
        sigcut = sigcut * fb

    try:
        f=pyslha.readSLHAFile ( slhafile )
    except pyslha.ParseError as e:
        logger.error ( "The file %s cannot be parsed as an SLHA file: %s" % (slhafile, e) )
        raise SModelSError()

    # Get cross section from file
    xSectionList = crossSection.getXsecFromSLHAFile(slhafile, useXSecs) 
    # Only use the highest order cross sections for each process
    xSectionList.removeLowerOrder()
    # Order xsections by PDGs to improve performance
    xSectionList.order()

    # Get maximum cross sections (weights) for single particles (irrespective
    # of sqrtS)
    maxWeight = {}
    for pid in xSectionList.getPIDs():
        maxWeight[pid] = xSectionList.getXsecsFor(pid).getMaxXsec()   

    # Generate dictionary, where keys are the PIDs and values 
    # are the list of cross sections for the PID pair (for performance)
    xSectionListDict = {}    
    for pids in xSectionList.getPIDpairs():
        xSectionListDict[pids] = xSectionList.getXsecsFor(pids)

    # Create 1-particle branches with all possible mothers
    
    branchList = []
    for pid in maxWeight:
        branchList.append(Branch())
        branchList[-1].BSMparticles = [[getObjectFromPdg(pid)]] # [[pid]]
        if not pid in BSMpdgs:
            logger.error ( "pid %d does not appear in BSMList" % pid )
        branchList[-1].maxWeight = maxWeight[pid]        

    # Generate final branches (after all R-odd particles have decayed)         
    
    finalBranchList = decayBranches(branchList, sigcut)

    # Generate dictionary, where keys are the PIDs and values are the list of branches for the PID (for performance)
    branchListDict = {}
    for branch in finalBranchList:
        if len(branch.BSMparticles) != 1:
            logger.error("During decomposition the branches should \
                            not have multiple PID lists!")
            return False   
        if branch.BSMparticles[0][0].pdg in branchListDict:
            branchListDict[branch.BSMparticles[0][0].pdg].append(branch)
        else:
            branchListDict[branch.BSMparticles[0][0].pdg] = [branch]

    for pid in xSectionList.getPIDs():
        if not pid in branchListDict: branchListDict[pid] = []

    #Sort the branch lists by max weight to improve performance:
    for pid in branchListDict:
        branchListDict[pid] = sorted(branchListDict[pid], 
                                     key=lambda br: br.maxWeight, reverse=True)                                                                       
    
    smsTopList = topology.TopologyList()
    # Combine pairs of branches into elements according to production
    # cross section list
    
    for pids in xSectionList.getPIDpairs():
        weightList = xSectionListDict[pids]
        minBR = (sigcut/weightList.getMaxXsec()).asNumber()

        if minBR > 1.: continue
        for branch1 in branchListDict[pids[0]]:
            BR1 =  branch1.maxWeight/maxWeight[pids[0]]  #Branching ratio for first branch            
            
            # add BRs as a list consisting of BR for prompt and a displaced decay
            if branch1.particles: BRs1, BR1labels = addPromptAndDisplaced(branch1,BR1)      
            else: 
                BRs1 = [BR1]
                BR1labels = ['longlived']
            
            if max(BRs1) < minBR: break  #Stop loop if the largest possible BR1 is already too low            
            for branch2 in branchListDict[pids[1]]:           
            
                BR2 =  branch2.maxWeight/maxWeight[pids[1]]  #Branching ratio for second branch                
                        
                
                if branch2.particles: BRs2, BR2labels = addPromptAndDisplaced(branch2,BR2)
                else: 
                    BRs2 = [BR2]
                    BR2labels = ['longlived']
                
                if max(BRs2) < minBR: break  #Stop loop if BR2 is already too low
                
                finalBRs = [ BR1*BR2 for BR1 in BRs1 for BR2 in BRs2]
                finalLabels = [ BR1label + '+' +  BR2label for BR1label in BR1labels for BR2label in BR2labels ] 
        
                for finalBR in finalBRs:
                    if type(finalBR) == type(1.*fb):
                        finalBR = finalBR.asNumber()
                    
                if max(finalBRs) < minBR: continue # Skip elements with xsec below sigcut

                if len(branch1.BSMparticles) != 1 or len(branch2.BSMparticles) != 1:
                    logger.error("During decomposition the branches should \
                            not have multiple PID lists!")
                    return False    

                newElement = element.Element([branch1, branch2])   

                weights = []
                for finalBR in finalBRs:
                    weights.append(weightList*finalBR.asNumber())

                
                newElement.weight = weightList*finalBRs[0].asNumber() # why asNumber necessary??
                newElement.weightList = weights
                newElement.weightLabels = finalLabels

                newElement.sortBranches()  #Make sure elements are sorted BEFORE adding them                            
                smsTopList.addElement(newElement)                          

                                                    
    smsTopList.compressElements(doCompress, doInvisible, minmassgap)
    smsTopList._setElementIds()                 
                       
    logger.debug("slhaDecomposer done in %.2f s." % (time.time() -t1 ) )
    
    
    return smsTopList



