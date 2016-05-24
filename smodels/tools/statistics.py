"""
.. module:: tools.statistics
   :synopsis: Code that computes CLs, p values, etc.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
        
"""

from smodels.tools.physicsUnits import fb
from smodels.tools.caching import _memoize

@_memoize
def upperLimit ( Nobs, Nexp, sigmaexp, lumi, alpha=.05 ):
    """ computes the 95% CL upper limit on the production cross section """

#     ret = _getUL(Nobs, Nexp, sigmaexp,alpha)/lumi
#     ret = _bayesianUpperLimit(Nobs,0.00001,Nexp,sigmaexp,1.-alpha)/lumi
    ret = _upperLimitMadAnalysis ( Nobs, Nexp, sigmaexp, 1.-alpha ) / lumi

    return ret

def _computeCLInterval( Nobs, Nexp, lumi, alpha=.05 ):
    """ Get experimental limit for the signal cross-section*efficiency in the analysis signal region.
                    
    :returns: (1-alpha) C.L. experimental upper limit for the signal cross-section in the signal
              region        
    """
    from scipy import stats
                   
    Nmax = 0.5*stats.chi2.isf(alpha,2*(Nobs+1)) - Nexp  #Upper limit on number of signal events
    maxSignalEvents = Nmax  #DOES NOT INCLUDE SYSTEMATIC UNCERTAINTIES
            
    return maxSignalEvents/lumi

def _bayesianUpperLimit ( nev, sac, xbg, sbg, cl=.95, prec=None, smax=None ):
    """ conway's bayesian method 
    :param nev: number of observed events
    :param sac: relative uncertainty in acceptance
    :param xbg: expected background
    :param sbg: uncertainty in background
    :param  cl: desired CL
    :param smax: upper limit of integration
    :param prec: integration step size """
    from smodels.tools import BayesianUpperLimit
    return BayesianUpperLimit.upperLimit ( nev, sac, xbg, sbg, cl, prec, smax )

def _upperLimitMadAnalysis ( nev, xbg, sbg, cl=.95, numberoftoys=10000, upto = 5.0, return_nan=False ):
    """ upper limit obtained via mad analysis 5 code 
    :param nev: number of observed events
    :param sac: relative uncertainty in acceptance
    :param sbg: uncertainty in background
    :param  cl: desired CL
    :param numberoftoys: how many toy experiments do we make?
    :param upto: search for 95% upto this number times nev
    """
    import exclusion_CLs
    import scipy.optimize
    def f( sig ):
        ##print "[statistics._upperLimitMadAnalysis] sig=",sig
        ##print "[statistics._upperLimitMadAnalysis] f=",exclusion_CLs.CLs ( nev, xbg, sbg, sig, numberoftoys )
        return exclusion_CLs.CLs ( nev, xbg, sbg, sig, numberoftoys ) - cl
    try:
        ##print "upto=",upto,"max(nv,xbg,sbg)=",max(nev,xbg,sbg)
        return scipy.optimize.brentq ( f, 0, upto * max(nev,xbg,sbg) )
    except Exception,e:
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
    from numpy import sqrt
    from scipy import integrate,special, stats
   
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
    from numpy import sqrt
    from scipy import optimize
    
    n0 = 0.
    n1 = abs(Nobs - Nbg) + 4*sqrt(NbgErr**2 + Nbg)
    while _getPValue(n1,Nobs,Nbg,NbgErr) > alpha:
        n1 += 4.*sqrt(NbgErr**2 + Nbg)

    def pValm(x):
        return alpha - _getPValue(x,Nobs,Nbg,NbgErr)

    nmax = optimize.brentq(pValm,n0,n1)
    return nmax
    
