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
def upperLimit ( Nobs, Nexp, sigmaexp, lumi, alpha=.05, toys=10000 ):
    """ computes the 95% CL upper limit on the production cross section """

#     ret = _getUL(Nobs, Nexp, sigmaexp,alpha)/lumi
#     ret = _bayesianUpperLimit(Nobs,0.00001,Nexp,sigmaexp,1.-alpha)/lumi
    ret = _upperLimitMadAnalysis ( Nobs, Nexp, sigmaexp, 1.-alpha, toys ) / lumi

    return ret

def _computeCLInterval( Nobs, Nexp, lumi, alpha=.05 ):
    """ Get experimental limit for the signal cross-section*efficiency in the analysis signal region.
                    
    :returns: (1-alpha) C.L. experimental upper limit for the signal cross-section in the signal
              region        
    """
                   
    Nmax = 0.5*stats.chi2.isf(alpha,2*(Nobs+1)) - Nexp  #Upper limit on number of signal events
    maxSignalEvents = Nmax  #DOES NOT INCLUDE SYSTEMATIC UNCERTAINTIES
            
    return maxSignalEvents/lumi


def _upperLimitMadAnalysis ( nev, xbg, sbg, cl=.95, numberoftoys=10000, upto = 5.0, 
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
        ##print "[statistics._upperLimitMadAnalysis] sig=",sig
        ##print "[statistics._upperLimitMadAnalysis] f=",exclusion_CLs.CLs ( nev, xbg, sbg, sig, numberoftoys )
        return CLs ( nev, xbg, sbg, sig, numberoftoys ) - cl
    try:
        ##print "upto=",upto,"max(nv,xbg,sbg)=",max(nev,xbg,sbg)
        return optimize.brentq ( f, 0, upto * max(nev,xbg,sbg) )
    except (ValueError,RuntimeError),e:
        if not return_nan:
            return _upperLimitMadAnalysis ( nev, xbg, sbg, cl, 5*numberoftoys, 5.0*upto, upto>30. )
        else:
            return float("nan")

def _getPValue(Nsig,Nobs,Nbg,NbgErr):
    """
    Computes the p-value using the signal cross-section (signalxsec) and the systematic
    error in the BG (bgsysError = systematic error/expected BG).
    Assumes a Gaussian distribution for the BG systematical error.
    """

   
    #Signal + BG prediction:
    Ntot = Nbg+Nsig
    #Total systematical error in BG+signal:
    NErr = NbgErr
            
    #Normalization        
    n = (1./2.)*(1. + special.erf(Ntot/(sqrt(2.)*NErr)))
            
    #P-value integrand
    def pint(x): 
        pInt = stats.poisson.cdf(Nobs,x)   #poisson.cdf with mean x (=total number of predicted events distributed according to gaussian)
        pInt *= stats.norm.pdf(x,loc=Ntot,scale=NErr)  #systematical error weight
        return pInt
    
    #P-value integral
    p = n*integrate.quad(pint,max(0.,Ntot-10.*NbgErr),Ntot+10.*NbgErr)[0]
            
    return p

def _getUL(Nobs,Nbg,NbgErr,alpha=0.05):
    """
    Computes the 95% upper limit on the signal*efficiency cross-section,
    given number of observed events (Nobs), number of expected BG events (Nbg)
    and the systematical error on the background (NbgErr)
    
    :param Nobs: number of observed events (integer)
    :param Nbg: number of expected BG events (float)
    :param NbgErr:  systematical error on the background (float)
    :param alpha: 1-C.L. i.e. for 95% C.L. alpha = 0.05
    :returns: the 95% confidence level on signal*efficiency*luminosity (float)
    """
   
    n0 = 0.
    n1 = abs(Nobs - Nbg) + 4*sqrt(NbgErr**2 + Nbg)
    while _getPValue(n1,Nobs,Nbg,NbgErr) > alpha:
        n1 += 4.*sqrt(NbgErr**2 + Nbg)

    def pValm(x):
        return alpha - _getPValue(x,Nobs,Nbg,NbgErr)

    nmax = optimize.brentq(pValm,n0,n1)
    return nmax
    
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

## print "cls=",CLs(100 ,100 ,0.0001,2.0,10000 )


def _bayesianUpperLimit ( nev, sac, xbg, sbg, cl=.95, prec=None, smax=None ):
    """ Computes the upper limit.

        :param nev: number of observed events
        :param sac: relative uncertainty in acceptance
        :param xbg: expected background
        :param sbg: uncertainty in background
        :param    cl: desired CL
        :param prec: integration step size
        :param smax: upper limit of integration
    """
    if smax == None: 
        smax = 10. * ( nev + 1 )
    if prec == None:
        prec=xbg*.001
        if prec > 200:
                prec=xbg*0.00001
    if xbg < 0.:
        print "[BayesianUpperLimit] error: cannot deal with negative expected background"
        return float('nan')

    if sac < 0. or sac > 1.:
        print "[BayesianUpperLimit] error: signal acceptance must be between 0 and 1"
        return float('nan')

    if nev < 0:
        print "[BayesianUpperLimit] error: negative number of observed events."
        return float('nan')

    if sbg < 0.:
        print "[BayesianUpperLimit] error: negative sigma on bkg expectation!"
        return float('nan')

    if ( cl < 0. or cl > 1. ):
        print "[BayesianUpperLimit] error, confidence limit must be between 0 and 1"
        return float('nan')

    if ( smax < 2. * nev or smax > 100. * nev ):
        print "[BayesianUpperLimit] warning: strange choice for smax, smax=",\
                    smax, ", nev=", nev
    if ( prec < 0. ):
        print "[BayesianUpperLimit] error: negative precision."
        return float('nan')

    if ( prec > .3 * nev and nev>0 ):
        print "[BayesianUpperLimit] error: precision too low, prec=",prec," nev=", nev
        return float('nan')

    bsum=0.
    xev=prec/2.
    blist=[0]*10000
    xlist=[0]*10000

    nlist =0

    while True:
        xlike = _blike ( nev,sac,xbg,sbg,xev)
        if ( math.isinf(xlike)) or math.isnan(xlike):
            if ( nlist==0 ):
                print "[BayesianUpperLimit] first likelihood is nan/inf! return nan!"
                return float('nan')

            print
            print "(D=" << nev << ", s=" << xev << ", l=" << xlike
        if nlist>=len(blist):
                print "[BayesianUpperLimit] need to extend array"
                for i in range(10000):
                        blist.append(0)
                        xlist.append(0)
        xlist[nlist]=xev
        blist[nlist]=xlike
        bsum+=xlike
        if False: 
            print "(D=" << nev << ", s=" << xev << ", l=" << xlike << ")\033[1A"
        xev+=prec
        if ( blist[nlist]/blist[0] < 1e-6 ): break
        if ( xev > smax ): break
        nlist+=1
    if False: print

    icl=0
    #    find place just below threshold
    bint=0.
    bcl=0.
#    print "blist=",blist
    for i in range(nlist):
        if ( ( bint < cl * bsum ) and ( (bint+blist[i])>cl*bsum) ):
#        if blist[i]>(1-cl) and blist[i+1]<(1-cl): ## fixme
            icl=i
            bcl=bint
        bint+=blist[i]

    # interpolate linearly
    plim=xlist[icl]+( xlist[icl+1]-xlist[icl] ) * (cl*bsum-bcl)/blist[icl+1]
    return plim

def _blike ( nev, sac, xbg, sbg, xev ):
        """ return likelihood to observe nev events given expected background xbg, 
                error on background sbg, expected number of signal events xev, 
                error on signal acceptance sac """
        xint=0.

        ## perform double Gaussian integral by Monte Carlo

        nmax=2000
        for i in range(nmax):
                # pick expected background and signal from Gaussian
                bg, sig = -1., -1.
                while bg < 0. or sig < 0.:
                        a,b = random.normal(), random.normal()
                        bg = xbg + a * sbg
                        sig = xev * ( 1. + b * sac )

                # total expected
                ex= bg + sig

                # value of integrand, poisson likelihood value
                # xxx = e (-ex) * ex^nev / nev! 
                xxx = exp(    nev * log ( ex ) - ex - math.lgamma(nev+1 ) )
#                 if math.isinf (xxx) or math.isnan ( xxx ):
#                         print "[blike] xxx=", xxx, " yyex=", yyex, " nev=", nev, ", nev!=", math.gamma(nev+1)
                xint +=xxx
        # print "[_blike] returns",xint/float(nmax)
        return xint/float(nmax)
