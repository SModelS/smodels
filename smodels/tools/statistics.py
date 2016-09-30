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
def upperLimit ( Nobs, Nexp, sigmaexp, lumi, alpha=.05, toys=200000 ):
    """ computes the 95% CL upper limit on the production cross section """
    #ret = _upperLimitMadAnalysis ( Nobs, Nexp, sigmaexp, 1.-alpha, toys ) / lumi
    #print "ret=",ret
    computer = UpperLimitComputer ( toys, lumi, 1.-alpha )
    ret = computer.compute ( Nobs, Nexp, sigmaexp, 5. )
    #print "ret=",ret
    return ret

class UpperLimitComputer:
    def __init__ ( self, numberoftoys, lumi, cl=.95):
        """
        :param numberoftoys: how many toy experiments do we make?
        :param lumi: integrated luminosity
        :param cl: desired CL
        """
        self.origNToys=numberoftoys
        self.currentNToys = self.origNToys
        self.lumi = lumi
        self.cl = cl

    def f( self, sig ):
        # self.currentNToys=int ( self.currentNToys * 1.10 )
        cls = CLs ( self.nev, self.xbg, self.sbg, sig, self.currentNToys ) - self.cl
        # print "[ULC] running",sig," with",self.currentNToys,"cls=",cls
        return cls

    def compute ( self, nev, xbg, sbg, upto=5.0, return_nan=False ):
        """ upper limit obtained via mad analysis 5 code 
        :param nev: number of observed events
        :param sac: relative uncertainty in acceptance
        :param xbg: expected bg
        :param sbg: uncertainty in background
        :param upto: defines the interval to be probed
        """
        self.nev = nev
        self.xbg = xbg
        self.sbg = sbg
        self.currentNToys = self.origNToys
        
        try:
            ## up = upto ## 5. ##  upto * max(nev,xbg,sbg)
            dn = max( 0, nev - xbg )
            up = upto * max( dn + math.sqrt(nev), dn + math.sqrt(xbg),sbg)
            #print "checking between 0 and ",up ## ,self.f(0),self.f(up)
            return optimize.brentq ( self.f, 0, up, rtol=1e-3 ) / self.lumi
        except (ValueError,RuntimeError),e:
            #print "exception: >>",type(e),e
            if not return_nan:
                # print "compute again, upto=",upto
                # self.origNToys = 5* self.origNToys
                return self.compute ( nev, xbg, sbg, 5.0*upto, upto>200. )
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
    ExpectedBGs = scipy.stats.norm.rvs( loc=ExpectedBG, scale=BGError, 
                                        size=NumToyExperiments )

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
