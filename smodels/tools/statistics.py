"""
.. module:: tools.statistics
   :synopsis: Code that computes CLs, p values, etc.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

from numpy import sqrt,inf
from scipy import stats,special,integrate

def computeCLInterval( Nobs, Nexp, lumi, alpha=.05 ):
    """ Get experimental limit for the signal cross-section*efficiency in the analysis signal region.
                    
    :returns: 1-alpha C.L. experimental upper limit for the signal cross-section in the signal
              region        
    """
                   
    Nmax = 0.5*stats.chi2.isf(alpha,2*(Nobs+1)) - Nexp  #Upper limit on number of signal events
    maxSignalEvents = Nmax  #DOES NOT INCLUDE SYSTEMATIC UNCERTAINTIES
            
    return maxSignalEvents/lumi

def getPValue(Nobs,Nbg,NbgErr,Nsig):
    """
    Computes the p-value using the signal cross-section (Nsig) and the systematic
    error in the BG (bgsysError = systematic error/expected BG).
    Assumes a Gaussian distribution for the BG systematical error.
    
    :param signalxsec: signal cross-section*efficiency for the signal region with units (Unum object)
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
