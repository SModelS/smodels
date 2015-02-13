"""
.. module:: tools.statistics
   :synopsis: Code that computes CLs, p values, etc.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

from numpy import sqrt,inf
from scipy import stats,special,integrate,optimize
import numpy as np
from scipy import stats,special,integrate

def computeCLInterval( Nobs, Nexp, lumi, alpha=.05 ):
    """ Get experimental limit for the signal cross-section*efficiency in the analysis signal region.
                    
    :returns: (1-alpha) C.L. experimental upper limit for the signal cross-section in the signal
              region        
    """
                   
    Nmax = 0.5*stats.chi2.isf(alpha,2*(Nobs+1)) - Nexp  #Upper limit on number of signal events
    maxSignalEvents = Nmax  #DOES NOT INCLUDE SYSTEMATIC UNCERTAINTIES
            
    return maxSignalEvents/lumi

def bayesianUpperLimit ( nev, sac, xbg, sbg, cl=.95, prec=-1., smax=0. ):
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

def getPValue(Nobs,Nbg,NbgErr,Nsig):
    """
    Computes the p-value using the signal cross-section (Nsig) and the systematic
    error in the BG (bgsysError = systematic error/expected BG).
    Assumes a Gaussian distribution for the BG systematical error.
    
    :param Nsig: number of signal events (float)
    :param Nobs: number of observed events (integer)
    :param Nbg: number of expected BG events (float)
    :param NbgErr:  systematical error on the background (float)
    :returns: the p-value (float)
    """
    #Signal + BG prediction:
    Ntot = Nbg+Nsig
    #Normalization        
    n = (1./2.)*(1. + special.erf(Ntot/(sqrt(2.)*NbgErr)))
    #P-value integrand
    def pint(x):            
        pInt = stats.poisson.cdf(Nobs,x)   #poisson.cdf with mean x (=total number of predicted events distributed according to gaussian)
        pInt *= stats.norm.pdf(x,loc=Ntot,scale=NbgErr)  #systematical error weight
        return pInt
    #P-value integral
    p = n*integrate.quad(pint,0,inf)[0]        
    return p

def getUL(Nobs,Nbg,NbgErr):
    """
    Computes the 95% upper limit on the signal*efficiency cross-section,
    given number of observed events (Nobs), number of expected BG events (Nbg)
    and the systematical error on the background (NbgErr)
    
    :param Nobs: number of observed events (integer)
    :param Nbg: number of expected BG events (float)
    :param NbgErr:  systematical error on the background (float)
    :returns: the 95% confidence level on signal*efficiency*luminosity (float)
    """
    
    def pValm(x):
        return 0.0227 - getPValue(x,Nobs,Nbg,NbgErr) 
    
    nmax = optimize.brentq(pValm,0.,float(Nobs))
    return nmax
    
