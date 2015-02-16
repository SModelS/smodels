#!/usr/bin/env python


#from smodels.tools.statistics import getUL,getPValue,computeCLInterval
import sys
from numpy import sqrt,inf
from scipy import stats,special,integrate,optimize
from smodels.tools.statistics import getUL

def getPValue(Nsig,Nobs,Nbg,NbgErr):
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
        p = n*integrate.quad(pint,Ntot-10.*NbgErr,Ntot+10.*NbgErr)[0]
        
        return p

def pValm(x):
    return 0.0227 - getPValue(x,Nobs,Nbg,NbgErr)

def agetUL(Nobs,Nbg,NbgErr):
    """
    Computes the 95% upper limit on the signal*efficiency cross-section,
    given number of observed events (Nobs), number of expected BG events (Nbg)
    and the systematical error on the background (NbgErr)
    
    :param Nobs: number of observed events (integer)
    :param Nbg: number of expected BG events (float)
    :param NbgErr:  systematical error on the background (float)
    :returns: the 95% confidence level on signal*efficiency*luminosity (float)
    """
    
    n0 = max((Nobs - Nbg) - 4*sqrt(NbgErr**2 + Nbg),0.)
    n1 = (Nobs - Nbg) + 4*sqrt(NbgErr**2 + Nbg)

    nmax = optimize.brentq(pValm,n0,n1)
    return nmax
    



Nobs = 4
Nbg = 6.
NbgErr = 4.1

#print computeCLInterval(Nobs, Nbg, 1.)

#sys.exit()

x = agetUL(Nobs,Nbg,NbgErr)
print x
x = getUL(Nobs,Nbg,NbgErr)
print x
#print getPValue(10.,Nobs,Nbg,NbgErr)


