"""
.. module:: reweighting
   :synopsis: calculate probabilities and relabel branches  

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
"""

import itertools
from math import exp
from smodels.tools.physicsUnits import GeV
from smodels.experiment.databaseParticles import jetList, lList
    

def defaultEffReweight(element,minWeight=1e-10):
    """
    Computes the lifetime reweighting factor for the element efficiency
    based on the lifetimes of all intermediate particles and the last stable odd-particle
    appearing in the element.
    The fraction corresponds to the fraction of decays corresponding to prompt
    decays to all intermediate BSM particles and to a long-lived decay (outside the detector)
    to the final BSM state.

    :param element: Element object
    :param minWeight: Lower cut for the reweighting factor. Any value below this will be taken to be zero.

    :return: Reweight factor (float)
    """

    elWidths = [ [ x.asNumber(GeV) for x in L ] for L in element.totalwidth ]
    elFraction = 1.
    for branchWidths in elWidths:
        branchFraction = 1.
        for width in branchWidths[:-1]:
            branchFraction *= calculateProbabilities(width)['F_prompt']
        branchFraction *= calculateProbabilities(width=branchWidths[-1])['F_long']
        elFraction *= branchFraction

    if elFraction < minWeight:
        return 0.

    return elFraction

def defaultULReweight(element):
    """
    Computes the lifetime reweighting factor for the element upper limit
    based on the lifetimes of all intermediate particles and the last stable odd-particle
    appearing in the element.
    The fraction corresponds to the fraction of decays corresponding to prompt
    decays to all intermediate BSM particles and to a long-lived decay (outside the detector)
    to the final BSM state.

    :param element: Element object

    :return: Reweight factor (float)
    """

    effFactor = defaultEffReweight(element)
    if not effFactor:
        return None
    else:
        return 1./effFactor


def addPromptAndDisplaced(branch):  
    """
    Add probabilities and branches, reweighted and relabeled by the probability functions to be longlived or MET or to decay promptly (within inner detector) or displaced (outside of inner detector but inside of outer detector) for each decay in the branch.
    :param branch: original branch
    :return: probabilities (depending on types of decay within branch), branches (with different labels depending on type of decay)
    """
      
    F = []
    for particle in branch.oddParticles:
        if isinstance(particle.totalwidth, list):
            width = min(particle.totalwidth)
        else: width = particle.totalwidth
        if width == float('inf')*GeV:
            F_long, F_prompt, F_displaced = [0.,1.,0.]
        elif width == 0.*GeV:
            F_long, F_prompt, F_displaced = [1.,0.,0.]
        else:    
            probs = calculateProbabilities(width.asNumber(GeV))
            F_long, F_prompt, F_displaced = probs['F_long'],probs['F_prompt'],probs['F_displaced']
        # allow for combinations of decays and a long lived particle only if the last BSM particle is the long lived one
        if particle == branch.oddParticles[-1]: F.append([F_long])
        else: F.append([F_prompt,F_displaced])
    
    # call the whole branch: 
    # .) prompt: all decays within this branch are prompt and the last particle is the LSP
    # .) longlived: all decays within this branch are prompt and the last particle is not the LSP
    # .) displaced: at least one decay is displaced no matter what the last particle is
    # discard others         
    promptValue = F[0][0]            
    for i in range(len(F)-1):
        promptValue *= F[i+1][0]  
    
    
    displacedP = []  
    lastP = F[-1][0]
    F.pop(-1)                                                   
    Plist = [list(P) for P in itertools.product(*F)]
    Plist.pop(0)                
    for P in Plist:
        value = P[0]
        for i in range(len(P)-1):
            value *= P[i+1]
        displacedP.append(value)
    displacedValue = sum(displacedP)
    displacedValue *= lastP
    
    probabilities = [promptValue, displacedValue]
    branches = labelPromptDisplaced(branch)         
        
    return probabilities, branches
    

def calculateProbabilities(width, l_inner = .001,
                           gb_inner=1.3,l_outer=10.,gb_outer=1.43):
    """
    The fraction of prompt and displaced decays are defined as:
    
    F_long = exp(-totalwidth*l_outer/gb_outer)
    F_prompt = 1 - exp(-totaltotalwidth*l_inner/gb_inner)
    F_displaced = 1 - F_prompt - F_long
    
    :param l_inner: is the inner radius of the detector, given in meters
    :param l_outer: is the outer radius of the detector, given in meters
    :param gb_inner: is the estimate for the kinematical factor gamma*beta for prompt decays.
    :param gb_outer: is the estimate for the kinematical factor gamma*beta for decays outside the detector.  
    :param width: particle width for which probabilities should be calculated (in GeV)
    
    :return: Dictionary with the probabilities for the particle not to decay (in the detector), to decay promptly or displaced.
    """
    
    # hc = 197.327*GeV*fm  #hbar * c
    hc = 197.327*1e-18  # hbar * c * m / GeV, expecting the width to be in GeV, and lengths in m

    F_long = exp(-width*l_outer/(gb_outer*hc))
    F_prompt = 1. - exp(-width*l_inner/(gb_inner*hc))
    F_displaced = 1. - F_prompt - F_long

    return {'F_long' : F_long, 'F_prompt' : F_prompt, 'F_displaced' : F_displaced}



def labelPromptDisplaced(branch):
    """
    Label displaced and prompt decays 
    :param branch: original branch
    :return: branches with correct labels  
    """
    promptBranch = branch.copy()
    if branch.oddParticles[-1].isMET():
        if not branch.evenParticles:
            promptBranch._decayType = 'METonly'
        else: promptBranch._decayType = 'prompt'

    else: promptBranch._decayType = 'longlived'
    
    displacedBranch = branch.copy()       
    if any(particle==jetList for vertex in branch.evenParticles for particle in vertex):
        displacedBranch._decayType = 'displacedJet'
    elif any(particle==lList for vertex in branch.evenParticles for particle in vertex):           
        displacedBranch._decayType = 'displacedLepton'
    else: displacedBranch._decayType = 'displaced(neither jet nor lepton)'
    
    branches = [promptBranch, displacedBranch]    
    return branches
        

