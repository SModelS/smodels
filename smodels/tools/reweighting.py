"""
.. module:: reweighting
   :synopsis: calculate probabilities and relabel branches  

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
"""

from math import exp
from smodels.tools.physicsUnits import GeV
from smodels.theory.element import Element
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger
    

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


def reweightFactorFor(element,resType='prompt'):
    """
    Computer the reweighting factor for the element according to the experimental result type.
    Currently only two result types are supported: 'prompt' and 'displaced'.
    If resultType = 'prompt', returns the reweighting factor for all decays in the element
    to be prompt and the last odd particle to be stable.
    If resultType = 'displaced', returns the reweighting factor for ANY decay in the element
    to be displaced and the last odd particle to be stable.

    :param element: Element object
    :param resType: Type of result to compute the reweight factor for (either 'prompt' or 'displaced')
    :return: probabilities (depending on types of decay within branch), branches (with different labels depending on type of decay)
    """

    if not isinstance(element,Element):
        logger.error('element should be an Element object and not %s' %type(element))
        raise SModelSError()
    rType = resType.lower()
    if not rType in ['prompt', 'displaced']:
        logger.error('resultType should be prompt or displaced and not %s' %str(resType))
        raise SModelSError()
    
    Fprompt = 1. #Keep track of the probability of all decays being prompt
    Fdisp = 1. #Keep track of the probability of AT LEAST one decay being displaced
    for branch in element.branches:
        for particle in branch.oddParticles[:-1]:
            if isinstance(particle.totalwidth, list):
                width = min(particle.totalwidth)
            else:
                width = particle.totalwidth
            if width == float('inf')*GeV:
                Fprompt *= 1.
                Fdisp *= 1.
            elif width == 0.*GeV:
                Fprompt *= 0.
                Fdisp *= 1.
            else:
                probs = calculateProbabilities(width.asNumber(GeV))
                Fprompt *= probs['F_prompt']
                Fdisp *= (probs['F_prompt'] + probs['F_long'])
    Fdisp = 1. - Fdisp #Prob(any decay being displaced) = 1 - Prob(all decays being either prompt or long-lived)
    #Include Flong factor for final BSM states:
    for branch in element.branches:
        width = branch.oddParticles[-1].totalwidth
        if width == float('inf')*GeV:
            Flong = 0.
        elif width == 0*GeV:
            Flong = 1.
        else:
            Flong = calculateProbabilities(width.asNumber(GeV))['F_long']
        Fdisp *= Flong
        Fprompt *= Flong
        
    if rType == 'prompt':
        return Fprompt
    if rType == 'displaced':
        return Fdisp

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
