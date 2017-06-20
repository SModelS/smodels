#!/usr/bin/env python

"""
.. module:: statistics
   :synopsis: Code that computes CLs, p values, etc.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
from smodels.tools.physicsUnits import fb
from smodels.tools.smodelsLogging import logger
from smodels.tools.caching import _memoize
from scipy import stats, optimize, integrate, special
from numpy import sqrt, exp, log, sign
import numpy
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
        cls = CLs ( self.nev, self.xbg, self.sbg, sig, self.currentNToys ) - self.cl
        return cls

    def compute ( self, nev, xbg, sbg, upto=5.0, return_nan=False ):
        if type (nev) not in [ list, tuple ]:
            return self.compute1d ( nev, xbg, sbg, upto, return_nan )
        """ upper limit obtained via mad analysis 5 code
        :param nev: number of observed events
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
        except (ValueError,RuntimeError) as e:
            #print "exception: >>",type(e),e
            if not return_nan:
                # print "compute again, upto=",upto
                # self.origNToys = 5* self.origNToys
                return self.compute ( nev, xbg, sbg, 5.0*upto, upto>200. )
            else:
                return float("nan")


    def compute1d ( self, nev, xbg, sbg, upto=5.0, return_nan=False ):
        """ upper limit obtained via mad analysis 5 code
        :param nev: number of observed events
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
        except (ValueError,RuntimeError) as e:
            #print "exception: >>",type(e),e
            if not return_nan:
                # print "compute again, upto=",upto
                # self.origNToys = 5* self.origNToys
                return self.compute ( nev, xbg, sbg, 5.0*upto, upto>200. )
            else:
                return float("nan")

def CLs(NumObserved, ExpectedBG, BGError, SigHypothesis, NumToyExperiments):
    """ CLs, with vectors and covariances """
    ## testing whether scipy is there
    if type(NumObserved) not in [ list, tuple ]:
        return CLs1D ( NumObserved, ExpectedBG, BGError, SigHypothesis, NumToyExperiments )
    try:
        import scipy.stats
    except ImportError:
        logger.warning('scipy is not installed... the CLs module cannot be used.')
        logger.warning('Please install scipy.')
        return False
    # generate a set of expected-number-of-background-events, one for each toy
    # experiment, distributed according to a Gaussian with the specified mean
    # and uncertainty
    ExpectedBGs = numpy.random.multivariate_normal( mean=ExpectedBG,
            cov=BGError, size=NumToyExperiments )

    ## discard all negatives
    ExpectedBGs = [value for value in ExpectedBGs if (value > 0).all() ]

    ## complain if too many negatives
    p = float(len(ExpectedBGs )) / NumToyExperiments
    if p < 0.9:
        logger.warning ( "only %d%s of points are all-positive" % (100*p,"%"))

    ToyBGs = list ( map ( scipy.stats.poisson.rvs, ExpectedBGs ) )

    # The probability for the background alone to fluctutate as LOW as
    # observed = the fraction of the toy experiments with backgrounds as low as
    # observed = p_b.
    # NB (1 - this p_b) corresponds to what is usually called p_b for CLs.
    p_b = scipy.stats.percentileofscore(ToyBGs, NumObserved, kind='weak')*.01

    # Toy MC for background+signal
    ExpectedBGandS = [expectedbg + SigHypothesis for expectedbg in ExpectedBGs]
    ToyBplusS = list ( map ( scipy.stats.poisson.rvs, ExpectedBGandS) )

    # Calculate the fraction of these that are >= the number observed,
    # giving p_(S+B). Divide by (1 - p_b) a la the CLs prescription.
    p_SplusB = scipy.stats.percentileofscore(ToyBplusS, NumObserved, kind='weak')*.01

    if p_SplusB>p_b:
        return 0.
    else:
        return 1.-(p_SplusB / p_b) # 1 - CLs


def CLs1D(NumObserved, ExpectedBG, BGError, SigHypothesis, NumToyExperiments):
    """ this method has been taken from MadAnalysis5, see
        Official website: <https://launchpad.net/madanalysis5>
        It is for the 1d case only .

        Thanks to the MadAnalysis for granting us the permission to use
        the code here! """

    ## testing whether scipy is there
    try:
        import scipy.stats
    except ImportError:
        logger.warning('scipy is not installed... the CLs module cannot be used.')
        logger.warning('Please install scipy.')
        return False
    # generate a set of expected-number-of-background-events, one for each toy
    # experiment, distributed according to a Gaussian with the specified mean
    # and uncertainty
    # numpy.random.multivariate_normal ( [1.,5.], [[1.0,0.0],[0.0,1.0]], 10 )
    ExpectedBGs = scipy.stats.norm.rvs( loc=ExpectedBG, scale=BGError,
                                        size=NumToyExperiments )
    #ExpectedBGs = numpy.random.normal( loc=ExpectedBG,
    #        scale=BGError, size=NumToyExperiments )

    # All negative coordinates are drawn again
    ExpectedBGs = [value for value in ExpectedBGs if value > 0]

    # For each toy experiment, get the actual number of background events by
    # taking one value from a Poisson distribution created using the expected
    # number of events.
    # print ( "ExpectedBGs=",ExpectedBGs)
    ToyBGs = map ( scipy.stats.poisson.rvs, ExpectedBGs )
    ToyBGs = list ( map(float, ToyBGs) )

    # The probability for the background alone to fluctutate as LOW as
    # observed = the fraction of the toy experiments with backgrounds as low as
    # observed = p_b.
    # NB (1 - this p_b) corresponds to what is usually called p_b for CLs.
    p_b = scipy.stats.percentileofscore(ToyBGs, NumObserved, kind='weak')*.01

    # Toy MC for background+signal
    ExpectedBGandS = [expectedbg + SigHypothesis for expectedbg in ExpectedBGs]
    ToyBplusS = map ( scipy.stats.poisson.rvs, ExpectedBGandS)
    ToyBplusS = list ( map(float, ToyBplusS) )

    # Calculate the fraction of these that are >= the number observed,
    # giving p_(S+B). Divide by (1 - p_b) a la the CLs prescription.
    p_SplusB = scipy.stats.percentileofscore(ToyBplusS, NumObserved, kind='weak')*.01

    if p_SplusB>p_b:
        return 0.
    else:
        return 1.-(p_SplusB / p_b) # 1 - CLs


def likelihood(nsig, nobs, nb, deltab, deltas):
        """
        Return the likelihood to observe nobs events given the
        predicted background nb, error on this background (deltab),
        expected number of signal events nsig and the error on the signal (deltas).

        :param nsig: predicted signal (float)
        :param nobs: number of observed events (float)
        :param nb: predicted background (float)
        :param deltab: uncertainty on background (float)
        :param deltas: uncertainty on signal (float)

        :return: likelihood to observe nobs events (float)

        """

        #Set signal error to 20%, if not defined
        if deltas is None:
            deltas = 0.2*nsig

        #     Why not a simple gamma function for the factorial:
        #     -----------------------------------------------------
        #     The scipy.stats.poisson.pmf probability mass function
        #     for the Poisson distribution only works for discrete
        #     numbers. The gamma distribution is used to create a
        #     continuous Poisson distribution.
        #
        #     Why not a simple gamma function for the factorial:
        #     -----------------------------------------------------
        #     The gamma function does not yield results for integers
        #     larger than 170. Since the expression for the Poisson
        #     probability mass function as a whole should not be huge,
        #     the exponent of the log of this expression is calculated
        #     instead to avoid using large numbers.


        #Define integrand (gaussian_(bg+signal)*poisson(nobs)):
        def prob(x,nsig, nobs, nb, deltab, deltas):
            poisson = exp(nobs*log(x) - x - math.lgamma(nobs + 1))
            gaussian = stats.norm.pdf(x,loc=nb+nsig,scale=sqrt(deltab**2 + deltas**2))

            return poisson*gaussian

        #Compute maximum value for the integrand:
        sigma2 = deltab**2 + deltas**2
        xm = nb + nsig - sigma2
        #If nb + nsig = sigma2, shift the values slightly:
        if xm == 0.:
            xm = 0.001
        xmax = xm*(1.+sign(xm)*sqrt(1. + 4.*nobs*sigma2/xm**2))/2.

        #Define initial integration range:
        nrange = 5.
        a = max(0.,xmax-nrange*sqrt(sigma2))
        b = xmax+nrange*sqrt(sigma2)
        like = integrate.quad(prob,a,b,(nsig, nobs, nb, deltab, deltas),
                                      epsabs=0.,epsrel=1e-3)[0]

        #Increase integration range until integral converges
        err = 1.
        while err > 0.01:
            like_old = like
            nrange = nrange*2
            a = max(0.,xmax-nrange*sqrt(sigma2))
            b = xmax+nrange*sqrt(sigma2)
            like = integrate.quad(prob,a,b,(nsig, nobs, nb, deltab, deltas),
                                      epsabs=0.,epsrel=1e-3)[0]
            err = abs(like_old-like)/like

        #Renormalize the likelihood to account for the cut at x = 0.
        #The integral of the gaussian from 0 to infinity gives:
        #(1/2)*(1 + Erf(mu/sqrt(2*sigma2))), so we need to divide by it
        #(for mu - sigma >> 0, the normalization gives 1.)
        norm = (1./2.)*(1. + special.erf((nb+nsig)/sqrt(2.*sigma2)))
        like = like/norm

        return like


def chi2(nsig, nobs, nb, deltab, deltas=None):
        """
        Computes the chi2 for a given number of observed events nobs
        given the predicted background nb, error on this background deltab,
        expected number of signal events nsig and, if given, the error on signal (deltas).
        If deltas is not given, assume an error of 20% on the signal.

        :param nsig: predicted signal (float)
        :param nobs: number of observed events (float)
        :param nb: predicted background (float)
        :param deltab: uncertainty in background (float)
        :param deltas: uncertainty in signal acceptance (float)

        :return: chi2 (float)

        """

        #Set signal error to 20%, if not defined
        if deltas is None:
            deltas = 0.2*nsig

        # Compute the likelhood for the null hypothesis (signal hypothesis) H0:
        llhd = likelihood(nsig, nobs, nb, deltab, deltas)

        #Percentual signal error:
        deltas_pct = deltas/(1.0*nsig)

        # Compute the maximum likelihood H1, which sits at nsig = nobs - nb
        # (keeping the same % error on signal):
        maxllhd = likelihood(nobs-nb, nobs, nb, deltab, deltas_pct*(nobs-nb))

        # Return infinite likelihood if it is zero
        # This can happen in case e.g. nb >> nobs
        if llhd == 0.:
            return float('inf')

        # Return the test statistic -2log(H0/H1)
        return -2*log(llhd/maxllhd)



if __name__ == "__main__":
    """
    f=open("bla.txt","w")
    for i in range(100):
        cls1 = CLs ( 10, 7., 2., 11.0, 20000 )
        cls2 = CLs ( [10], [7.], [[4.0]], [11.0], 20000 )
        print ( "CLs =", cls1 )
        print ( "CLsC=", cls2 )
        f.write ( "%.3f %.3f\n" % ( cls1, cls2 ) )
    f.close()
    """
    #print ( "CLsC=", CLsCov ( [25,15], [7.,15.], [[1.415,0.],[0.,10.]], [11.0,14.00], 100 ) )
    # import doctest
    # doctest.testmod()
    print ( upperLimit ( [4], [3.6], [[0.1**2]], 20. / fb, .05, 1000 ) )
    print ( upperLimit ( 4, 3.6, 0.1, 20. / fb, .05, 1000 ) )
