"""
.. module:: reweighting
   :synopsis: calculate probabilities and relabel branches  

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
"""

import copy
import itertools
from math import exp
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import MeV, GeV, m, mm, fm
from smodels.experiment.finalStateParticles import jetList, lList

        

def addPromptAndDisplaced(branch):  
    """
    Add probabilities and branches, reweighted and relabeled by the probability functions to be longlived or MET or to decay promptly (within inner detector) or displaced (outside of inner detector but inside of outer detector) for each decay in the branch.
    :param branch: original branch
    :return: probabilities (depending on types of decay within branch), branches (with different labels depending on type of decay)
    """
    
    if not branch.particles: # no decays happened
        if branch.BSMparticles[0].totalwidth == 0.*GeV: 
            branch._decayType = 'prompt'
            probabilities = [1.]
            branches = [branch]
            return probabilities, branches 
   
    F = []
    for particle in branch.BSMparticles:
        if particle.totalwidth == float('inf')*GeV:
            F_long, F_prompt, F_displaced = [0.,1.,0.]
        else:
            if particle.totalwidth == 0.*GeV: continue
            F_long, F_prompt, F_displaced = calculateProbabilities(particle)
        # allow for combinations of decays and a long lived particle only if the last BSM particle is the long lived one
        if F_long and particle == branch.BSMparticles[-1]: F.append([F_long])
        else: F.append([F_prompt,F_displaced])
    
    # call the whole branch: 
    # .) prompt: all decays within this branch are prompt
    # .) longlived: at least one particle is longlived (and the last one in the decay)
    # .) displaced: at least one decay is displaced and no longlived particle
    # discard others
         
    if len(F[-1]) == 1:  # last BSM particle is long lived
        branch._decayType = 'longlived'
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
    
    F_long = exp(-totalwidth*l_outer/gb_outer)
    F_prompt = 1 - exp(-totaltotalwidth*l_inner/gb_inner)
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
    
    F_long = exp( -1*particle.totalwidth * l_outer /(gb_outer*hc) )    
    F_prompt = 1. - exp( -1*particle.totalwidth * l_inner /(gb_inner*hc) )         
    F_displaced = 1. - F_prompt - F_long
    
    return F_long, F_prompt, F_displaced



def labelPromptDisplaced(branch):
    """
    Label displaced and prompt decays 
    :param branch: original branch
    :return: branches with correct labels  
    """
    promptBranch = branch.copy()
    promptBranch._decayType = 'prompt'
    
    displacedBranch = branch.copy()      
    if any(particle==jetList for particle in branch.particles ):
        displacedBranch._decayType = 'displacedJet'
    elif any(particle==lList for particle in branch.particles ):           
        displacedBranch._decayType = 'displacedLepton'
    else: displacedBranch._decayType = 'displaced(neither jet nor lepton)'
    
    branches = [promptBranch, displacedBranch]    
    return branches
        

