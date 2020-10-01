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

Leff_inner_default = 0.000769
Leff_outer_default = 7.0

def defaultEffReweight(element,Leff_inner=None,Leff_outer=None,
                                minWeight=1e-10):
    """
    Computes the lifetime reweighting factor for the element efficiency
    based on the lifetimes of all intermediate particles and the last stable odd-particle
    appearing in the element.
    The fraction corresponds to the fraction of decays corresponding to prompt
    decays to all intermediate BSM particles and to a long-lived decay (outside the detector)
    to the final BSM state.

    :param element: Element object or nested list of widths
    :param minWeight: Lower cut for the reweighting factor. Any value below this will be taken to be zero.
    :param Leff_inner: is the effective inner radius of the detector, given in meters. If None,
                        use default value.
    :param Leff_outer: is the effective outer radius of the detector, given in meters. If None,
                        use default value.


    :return: Reweight factor (float)
    """

    if Leff_inner is None:
        Leff_inner = Leff_inner_default
    if Leff_outer is None:
        Leff_outer = Leff_outer_default


    if isinstance(element,list):
        elWidths = [[w.asNumber(GeV) for w in br] for br in element]
    else:
        elWidths = [ [ x.asNumber(GeV) for x in L ] for L in element.totalwidth ]
    elFraction = 1.
    for branchWidths in elWidths:
        branchFraction = 1.
        for width in branchWidths[:-1]:
            branchFraction *= calculateProbabilities(width,Leff_inner,Leff_outer)['F_prompt']
        branchFraction *= calculateProbabilities(branchWidths[-1],
                                                  Leff_inner,Leff_outer)['F_long']
        elFraction *= branchFraction

    if elFraction < minWeight:
        return 0.

    return elFraction

def defaultULReweight(element,Leff_inner=None,Leff_outer=None):
    """
    Computes the lifetime reweighting factor for the element upper limit
    based on the lifetimes of all intermediate particles and the last stable odd-particle
    appearing in the element.
    The fraction corresponds to the fraction of decays corresponding to prompt
    decays to all intermediate BSM particles and to a long-lived decay (outside the detector)
    to the final BSM state.

    :param element: Element object
    :param Leff_inner: is the effective inner radius of the detector, given in meters. If None,
                        use default value.
    :param Leff_outer: is the effective outer radius of the detector, given in meters. If None,
                        use default value.


    :return: Reweight factor (float)
    """

    if Leff_inner is None:
        Leff_inner = Leff_inner_default
    if Leff_outer is None:
        Leff_outer = Leff_outer_default


    effFactor = defaultEffReweight(element,Leff_inner,Leff_outer)
    if not effFactor:
        return None
    else:
        return 1./effFactor


def reweightFactorFor(element,resType='prompt',
                        Leff_inner=None,Leff_outer=None):
    """
    Computer the reweighting factor for the element according to the experimental result type.
    Currently only two result types are supported: 'prompt' and 'displaced'.
    If resultType = 'prompt', returns the reweighting factor for all decays in the element
    to be prompt and the last odd particle to be stable.
    If resultType = 'displaced', returns the reweighting factor for ANY decay in the element
    to be displaced and no long-lived decays and the last odd particle to be stable.
    Not that the fraction of "long-lived (meta-stable) decays" is usually included
    in topologies where the meta-stable particle appears in the final state. Hence
    it should not be included in the prompt or displaced fractions.

    :param element: Element object
    :param resType: Type of result to compute the reweight factor for (either 'prompt' or 'displaced')
    :param Leff_inner: is the effective inner radius of the detector, given in meters. If None,
                        use default value.
    :param Leff_outer: is the effective outer radius of the detector, given in meters. If None,
                        use default value.

    :return: probabilities (depending on types of decay within branch), branches (with different labels depending on type of decay)
    """

    if Leff_inner is None:
        Leff_inner = Leff_inner_default
    if Leff_outer is None:
        Leff_outer = Leff_outer_default


    if not isinstance(element,Element):
        logger.error('element should be an Element object and not %s' %type(element))
        raise SModelSError()
    rType = resType.lower()
    if not rType in ['prompt', 'displaced']:
        logger.error('resultType should be prompt or displaced and not %s' %str(resType))
        raise SModelSError()

    FpromptTotal = 1. #Keep track of the probability of all decays being prompt
    FlongTotal = 1. #Keep track of the probability for at least one decay being long-lived
    for branch in element.branches:
        for particle in branch.oddParticles[:-1]:
            if isinstance(particle.totalwidth, list):
                width = min(particle.totalwidth)
            else:
                width = particle.totalwidth
            if width == float('inf')*GeV:
                Fprompt = 1.
                Fdisp = 0.
            elif width == 0.*GeV:
                Fprompt = 0.
                Fdisp = 0.
            else:
                probs = calculateProbabilities(width.asNumber(GeV),
                                                Leff_inner,Leff_outer)
                Fprompt = probs['F_prompt']
                Fdisp = probs['F_displaced']
            FlongTotal *= (Fprompt+Fdisp)
            FpromptTotal *= Fprompt
    #Prob(at least one decay being long-lived) = 1 - Prob(all decays being either prompt or displaced)
    FlongTotal = 1. - FlongTotal
    #Prob(at least one displaced decay and no long-lived decays) = 1 -Prob(at least one long-lived) - Prob(all prompt)
    FdispTotal = 1. - FlongTotal - FpromptTotal
    #Include Flong factor for final BSM states:
    for branch in element.branches:
        width = branch.oddParticles[-1].totalwidth
        if width == float('inf')*GeV:
            Flong = 0.
        elif width == 0*GeV:
            Flong = 1.
        else:
            Flong = calculateProbabilities(width.asNumber(GeV),
                                            Leff_inner,Leff_outer)['F_long']
        FdispTotal *= Flong
        FpromptTotal *= Flong

    if rType == 'prompt':
        return FpromptTotal
    if rType == 'displaced':
        return FdispTotal

def calculateProbabilities(width, Leff_inner, Leff_outer):
    """
    The fraction of prompt and displaced decays are defined as:

    F_long = exp(-totalwidth*l_outer/gb_outer)
    F_prompt = 1 - exp(-totaltotalwidth*l_inner/gb_inner)
    F_displaced = 1 - F_prompt - F_long

    :param Leff_inner: is the effective inner radius of the detector, given in meters
    :param Leff_outer: is the effective outer radius of the detector, given in meters
    :param width: particle width for which probabilities should be calculated (in GeV)

    :return: Dictionary with the probabilities for the particle not to decay (in the detector), to decay promptly or displaced.
    """

    # hc = 197.327*GeV*fm  #hbar * c
    hc = 197.327*1e-18  # hbar * c * m / GeV, expecting the width to be in GeV, and lengths in m

    F_long = exp(-width*Leff_outer/hc)
    F_prompt = 1. - exp(-width*Leff_inner/hc)
    F_displaced = 1. - F_prompt - F_long

    return {'F_long' : F_long, 'F_prompt' : F_prompt, 'F_displaced' : F_displaced}
