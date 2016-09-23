#!/usr/bin/env python

"""
.. module:: statistics
   :synopsis: Code that computes CLs, p values, etc.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
        
"""

import logging
from smodels.tools.physicsUnits import fb
from smodels.tools.caching import _memoize
from scipy import stats, optimize, integrate,special
from numpy import sqrt, random, exp, log
import math

@_memoize
def upperLimit ( Nobs, Nexp, sigmaexp, lumi, alpha=.05, toys=50000 ):
    """ computes the 95% CL upper limit on the production cross section """
    ret = _upperLimitMadAnalysis ( Nobs, Nexp, sigmaexp, 1.-alpha, toys ) / lumi
    return ret

def _upperLimitMadAnalysis ( nev, xbg, sbg, cl=.95, numberoftoys=50000, upto = 5.0, 
                             return_nan=False ):
    """ upper limit obtained via mad analysis 5 code 
    :param nev: number of observed events
    :param sac: relative uncertainty in acceptance
    :param sbg: uncertainty in background
    :param  cl: desired CL
    :param numberoftoys: how many toy experiments do we make?
    :param upto: search for 95% upto this number times nev

    """
    
    def f( sig ):
        return CLs ( nev, xbg, sbg, sig, numberoftoys ) - cl
    try:
        return optimize.brentq ( f, 0, upto * max(nev,xbg,sbg) )
    except (ValueError,RuntimeError),e:
        if not return_nan:
            return _upperLimitMadAnalysis ( nev, xbg, sbg, cl, 5*numberoftoys, 
                                            5.0*upto, upto>30. )
        else:
            return float("nan")

def CLs(NumObserved, ExpectedBG, BGError, SigHypothesis, NumToyExperiments):
    """ this method has been taken from MadAnalysis5, see
        Official website: <https://launchpad.net/madanalysis5>

        Thanks to the MadAnalysis for granting us the permission to use
        the code here! """
        
    ## testing whether scipy is there
    try:
        import scipy.stats
    except ImportError:
        logging.warning('scipy is not installed... the CLs module cannot be used.')
        logging.warning('Please install scipy.')
        return False
    # generate a set of expected-number-of-background-events, one for each toy
    # experiment, distributed according to a Gaussian with the specified mean
    # and uncertainty
    ExpectedBGs = scipy.stats.norm.rvs(loc=ExpectedBG, scale=BGError, size=NumToyExperiments)

    # Ignore values in the tail of the Gaussian extending to negative numbers
    ExpectedBGs = [value for value in ExpectedBGs if value > 0]

    # For each toy experiment, get the actual number of background events by
    # taking one value from a Poisson distribution created using the expected
    # number of events.
    ToyBGs = scipy.stats.poisson.rvs(ExpectedBGs)
    ToyBGs = map(float, ToyBGs)

    # The probability for the background alone to fluctutate as LOW as
    # observed = the fraction of the toy experiments with backgrounds as low as
    # observed = p_b.
    # NB (1 - this p_b) corresponds to what is usually called p_b for CLs.
    p_b = scipy.stats.percentileofscore(ToyBGs, NumObserved, kind='weak')*.01

    # Toy MC for background+signal
    ExpectedBGandS = [expectedbg + SigHypothesis for expectedbg in ExpectedBGs]
    ToyBplusS = scipy.stats.poisson.rvs(ExpectedBGandS)
    ToyBplusS = map(float, ToyBplusS)

    # Calculate the fraction of these that are >= the number observed,
    # giving p_(S+B). Divide by (1 - p_b) a la the CLs prescription.
    p_SplusB = scipy.stats.percentileofscore(ToyBplusS, NumObserved, kind='weak')*.01

    if p_SplusB>p_b:
        return 0.
    else:
        return 1.-(p_SplusB / p_b) # 1 - CLs

if __name__ == "__main__":
    print upperLimit ( 4, 3.6, 0.1, 20. / fb )
