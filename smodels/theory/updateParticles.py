"""
.. module:: particles
   :synopsis: Updates masses, width and branches from SLHA file

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
"""

import pyslha
import copy
import itertools
from math import exp
from smodels.particleDefinitions import SMpdgs, BSMList, BSMpdgs
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import GeV, MeV, m, mm, fm
from smodels.theory.particleClass import particleInList
from smodels.SMparticleDefinitions import jetList, lList


def updateParticles(slhafile, BSMList):
    getMassWidthBranches(slhafile, BSMList)
    addLongLived(BSMList)
    

def getMassWidthBranches(slhafile, BSMList):
    """
    Update mass, total width and branches of BSM particles from input slha file. 
        
    :param BSMList: list of particle objects containing all BSM particles defined in particleDefinitions
        
    """
    res = pyslha.readSLHAFile(slhafile)

    rOdd = BSMpdgs
    rEven = SMpdgs

    writeIgnoreMessage ( res.decays.keys(), rEven, rOdd )

	    
    for pid in res.decays.keys():
        for particle in BSMList:
            if abs(particle.pdg) == pid:
                particle.mass = abs(res.blocks['MASS'][pid])*GeV
                particle.width = abs(res.decays[pid].totalwidth)*GeV

                if not pid in rOdd:
                    logger.warning("Particle %s is not in BSMList, its mass, width and branches cannot be included") %pid
                    continue
                brs = []
                for decay in res.decays[pid].decays:
                    nEven = nOdd = 0.
                    for pidd in decay.ids:
                        if pidd in rOdd: nOdd += 1
                        elif pidd in rEven: nEven += 1
                        else:
                            logger.warning("Particle %i not defined in particleDefinitions.py,decay %i -> [%s] will be ignored" %(pidd,pid,decay.ids))
                            break
                    if nOdd + nEven == len(decay.ids) and nOdd == 1:
                        brs.append(decay)
                    else:
                        logger.info("Ignoring decay: %i -> [%s]",pid,decay.ids)
                
                if particle.pdg == pid: particle.branches = brs
                elif (-1)*particle.pdg == pid:
                    brsConj = copy.deepcopy(brs)
                    for br in brsConj:
                        br.ids = [-x for x in br.ids]                              
                    particle.branches = brsConj


def writeIgnoreMessage ( keys, rEven, rOdd ):
    msg = ""
    for pid in keys:
        if not pid in list(rEven) + list(rOdd):
            logger.warning("Particle %i not defined in particleDefinitions.py, its decays will be ignored" %(pid))
            continue
        if pid in rEven:
            msg += "%i, " %pid # used to be string of particle name (getObjectFromPdg(pid).label causes circular dependence) -- pid not so pretty 
            continue         
    if len(msg)>0:
            logger.info ( "Ignoring %s decays" % msg[:-2] )
            
            
            
def addLongLived(BSMList):
    """
    Update BSM particles including whether a particle can be long-lived.
    
    Using the widths of the particles (see getMassWidthBranches), reweights the BR's with the fraction
    of decays within the detector and add the fraction of "long-lived decays".
    The fraction of decays within the detector and "long-lived decays" are defined as:
  
    F_long = exp(-width*l_outer/gb_outer)
    F_decay = 1 - F_long
    
    where l_outer is the outer radius of the detector
    and gb_outer is the estimate for the kinematical factor gamma*beta.
    We use l_outer=10.m and gb_outer= 0.6    
    """
    
    l_outer = 10.*m
    gb_outer = 0.6    
    
    hc = 197.327*MeV*fm  #hbar * c
    
    for particle in BSMList:       
        F_long = exp( -1*particle.width * l_outer /(gb_outer*hc) )
        """
        F_decay = (1. - F_long) 
        
        
        for branch in particle.branches:
            branch.br *= F_decay            
        """
        if F_long:
            stable = pyslha.Decay(br=1., nda=0 ,ids=[]) #F_long
            particle.branches.append(stable)


def addPromptAndDisplaced(branch):  
    """
    Add probabilities and branches, reweighted and relabeled by the probability functions to be longlived or MET or to decay promptly (within inner detector) or displaced (outside of inner detector but inside of outer detector) for each decay in the branch.
    :param branch: original branch
    :return: probabilities (depending on types of decay within branch), branches (with different labels depending on type of decay)
    """
    
    if not branch.particles: # no decays happened
        if branch.BSMparticles[0][0].isStable(): 
            branch.decayType = 'METonly'
            probabilities = [1.]
            branches = [branch]
            return probabilities, branches 
   
    F = []
    for particle in branch.BSMparticles[0]:
        if particle.isStable(): continue
        F_long, F_prompt, F_displaced = calculateProbabilities(particle)
        # allow for combinations of decays and a long lived particle only if the last BSM particle is the long lived one
        if F_long and particle == branch.BSMparticles[0][-1]: F.append([F_long])
        else: F.append([F_prompt,F_displaced])
    
    # call the whole branch: 
    # .) prompt: all decays within this branch are prompt
    # .) longlived: at least one particle is longlived (and the last one in the decay)
    # .) displaced: at least one decay is displaced and no longlived particle
    # discard others
         
    if len(F[-1]) == 1:  # last BSM particle is long lived
        branch.decayType = 'longlived'
        branches = [branch]            
        longValue = F[-1][0] 
        
        if len(F) > 1: 
        # decays before can have all combinations of prompt and displaced 
        #(only different from 1. if particle before has probability to be longlived)  
            F.pop(-1)
            Plist = [list(P) for P in itertools.product(*F)]
            F = []        
            for P in Plist:
                value = P[0]
                for i in range(len(P)-1):
                    value *= P[i+1]
                F.append(value)
            longValue *= sum(F)
          
        probabilities = [longValue]
    
    elif len(F[-1]) > 1: # last BSM particle (not MET) decayed    
        promptValue = F[0][0]            
        for i in range(len(F)-1):
            promptValue *= F[i+1][0]  
        
        displacedP = []                                                      
        Plist = [list(P) for P in itertools.product(*F)]
        Plist.pop(0)                
        for P in Plist:
            value = P[0]
            for i in range(len(P)-1):
                value *= P[i+1]
            displacedP.append(value)
        displacedValue = sum(displacedP)
        
        probabilities = [promptValue, displacedValue]
        branches = labelPromptDisplaced(branch)         

    else: logger.warning("No probabilities were calculated for the decay(s) in s%" %(branch))
        
    return probabilities, branches
    

def calculateProbabilities(particle):
    """
    The fraction of prompt and displaced decays are defined as:
    
    F_long = exp(-width*l_outer/gb_outer)
    F_prompt = 1 - exp(-width*l_inner/gb_inner)
    F_displaced = 1 - F_prompt - F_long
    
    where l_inner (l_outer) is the inner (outer) radius of the detector
    and gb_inner (gb_outer) is the estimate for the kinematical factor gamma*beta.
    We use l_inner=10.mm and gb_inner=10; l_outer=10.m and gb_outer=0.6.  
    :param particle: particle object for which probabilities should be calculated
    :return: probabilities for the particle not to decay (in the detector), to decay promptly or displaced.
    """
    
    l_inner = 10.*mm
    gb_inner = 10.    
    l_outer = 10.*m
    gb_outer = 0.6
    hc = 197.327*MeV*fm  #hbar * c
    
    F_long = exp( -1*particle.width * l_outer /(gb_outer*hc) )    
    F_prompt = 1. - exp( -1*particle.width * l_inner /(gb_inner*hc) )         
    F_displaced = 1. - F_prompt - F_long
    
    return F_long, F_prompt, F_displaced



def labelPromptDisplaced(branch):
    """
    Label displaced and prompt decays 
    :param branch: original branch
    :return: branches with correct labels  
    """
    promptBranch = branch.copy()
    promptBranch.decayType = 'prompt'
    
    displacedBranch = branch.copy()      
    if any(particleInList(particle,[jetList]) for particle in branch.getBranchParticles() ):
        displacedBranch.decayType = 'displacedJet'
    elif any(particleInList(particle,[lList]) for particle in branch.getBranchParticles() ):           
        displacedBranch.decayType = 'displacedLepton'
    else: displacedBranch.decayType = 'displaced(neither jet nor lepton)'
    
    branches = [promptBranch, displacedBranch]    
    return branches
        

