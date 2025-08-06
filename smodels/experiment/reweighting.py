"""
.. module:: reweighting
   :synopsis: calculate probabilities and relabel branches

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
"""

from math import exp
from smodels.base.physicsUnits import GeV
from smodels.base.genericSMS import GenericSMS
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from smodels.base.smodelsLogging import logger

Leff_inner_default = 0.000769
Leff_outer_default = 7.0


def getWidthsFromSMS(sms):
    """
    Extracts all the widths of unstable particles in the SMS and the widths
    of BSM particles appearing as final states (undecayed).

    :param sms: SMS object

    :return: List of unstable widths and list of stable widths
    """

    if not isinstance(sms, GenericSMS):
        msgError = "Input to getWidthsFromSMS should be an SMS object"
        logger.error(msgError)
        raise SModelSError(msgError)

    # Get the unstable and stable widths:
    unstableWidths = []
    stableWidths = []
    for nodeIndex in sms.nodeIndices:
        if nodeIndex == sms.rootIndex:
            continue
        # Convert from index to node object
        node = sms.indexToNode(nodeIndex)
        if node.isInclusive:
            continue

        if sms.out_degree(nodeIndex) != 0:
            width = node.totalwidth
            if isinstance(width, list):  # If width is a list (MultiParticle obj), use smallest
                width = min(width)
            unstableWidths.append(width)
        else:
            if node.isSM:
                continue
            width = node.totalwidth
            if isinstance(width, list):  # If width is a list (MultiParticle obj), use largest
                width = max(width)
            stableWidths.append(width)

    return unstableWidths, stableWidths


def defaultEffReweight(sms=None, unstableWidths=None, stableWidths=None,
                       Leff_inner=None, Leff_outer=None,
                       minWeight=1e-10):
    """
    Computes the lifetime reweighting factor for the SMS efficiency
    based on the widths of the BSM particles.
    The widths can be defined through the unstableWidths and stableWidths arguments
    or they will be extracted from the SMS.
    The reweighting factor corresponds to the fraction of prompt decays (for the unstableWidths)
    and the fractions of detector-stable decays (for the stableWidths).

    :param sms: SMS object
    :param unstableWidths: List of widths for particles appearing as prompt decays
    :param stableWidths: List of widths for particles appearing as stable
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

    if unstableWidths is None or stableWidths is None:
        if sms is None or not isinstance(sms, GenericSMS):
            logger.error("If unstableWidths and stableWidths are not defined, the SMS should be given.")
            raise SModelSError()
        # Extract SMS widths:
        unstableWidths, stableWidths = getWidthsFromSMS(sms)

    unstableWidths = [w.asNumber(GeV) for w in unstableWidths[:]]
    stableWidths = [w.asNumber(GeV) for w in stableWidths[:]]

    elFraction = 1.
    for width in unstableWidths:
        elFraction *= calculateProbabilities(width, Leff_inner, Leff_outer)['F_prompt']

    for width in stableWidths:
        elFraction *= calculateProbabilities(width,
                                             Leff_inner, Leff_outer)['F_long']

    if elFraction < minWeight:
        return 0.

    return elFraction


def defaultULReweight(sms=None, unstableWidths=None, stableWidths=None,
                      Leff_inner=None, Leff_outer=None):
    """
    Computes the lifetime reweighting factor for the SMS upper limit
    based on the lifetimes of all intermediate particles and the last stable odd-particle
    appearing in the SMS.
    The fraction corresponds to the fraction of decays corresponding to prompt
    decays to all intermediate BSM particles and to a long-lived decay (outside the detector)
    to the final BSM state.

    :param sms: SMS object
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

    effFactor = defaultEffReweight(sms=sms, unstableWidths=unstableWidths,
                                   stableWidths=stableWidths,
                                   Leff_inner=Leff_inner, Leff_outer=Leff_outer)
    if not effFactor:
        return None
    else:
        return 1./effFactor


def reweightFactorFor(sms=None, resType='prompt',
                      unstableWidths=None, stableWidths=None,
                      Leff_inner=None, Leff_outer=None):
    """
    Computer the reweighting factor for the SMS according to the experimental result type.
    Currently only two result types are supported: 'prompt' and 'displaced'.
    If resultType = 'prompt', returns the reweighting factor for all decays in the SMS
    to be prompt and the last odd particle to be stable.
    If resultType = 'displaced', returns the reweighting factor for ANY decay in the SMS
    to be displaced and no long-lived decays and the last odd particle to be stable.
    Not that the fraction of "long-lived (meta-stable) decays" is usually included
    in topologies where the meta-stable particle appears in the final state. Hence
    it should not be included in the prompt or displaced fractions.

    :param sms: SMS object
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

    if unstableWidths is None or stableWidths is None:
        if sms is None or not isinstance(sms, GenericSMS):
            logger.error("If unstableWidths and stableWidths are not defined, the SMS should be given.")
            raise SModelSError()
        # Extract SMS widths:
        unstableWidths, stableWidths = getWidthsFromSMS(sms)

    rType = resType.lower()
    if rType not in ['prompt', 'displaced']:
        logger.error(f'resultType should be prompt or displaced and not {str(resType)}')
        raise SModelSError()

    FpromptTotal = 1.  # Keep track of the probability of all decays being prompt
    FlongTotal = 1.  # Keep track of the probability for at least one decay being long-lived
    for width in unstableWidths:
        if width == float('inf')*GeV:
            Fprompt = 1.
            Fdisp = 0.
        elif width == 0.*GeV:
            Fprompt = 0.
            Fdisp = 0.
        else:
            probs = calculateProbabilities(width.asNumber(GeV),
                                           Leff_inner, Leff_outer)
            Fprompt = probs['F_prompt']
            Fdisp = probs['F_displaced']
        FlongTotal *= (Fprompt+Fdisp)
        FpromptTotal *= Fprompt

    # Prob(at least one decay being long-lived) = 1 - Prob(all decays being either prompt or displaced)
    FlongTotal = 1. - FlongTotal
    # Prob(at least one displaced decay and no long-lived decays) = 1 -Prob(at least one long-lived) - Prob(all prompt)
    FdispTotal = 1. - FlongTotal - FpromptTotal
    # Include Flong factor for final BSM states:
    for width in stableWidths:
        if width == float('inf')*GeV:
            Flong = 0.
        elif width == 0*GeV:
            Flong = 1.
        else:
            Flong = calculateProbabilities(width.asNumber(GeV),
                                           Leff_inner, Leff_outer)['F_long']
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

    return {'F_long': F_long, 'F_prompt': F_prompt, 'F_displaced': F_displaced}
