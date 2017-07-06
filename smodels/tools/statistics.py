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
from numpy import sqrt, exp, log, sign, array, matrix
import numpy
import math
import sys

@_memoize
def upperLimit ( Nobs, Nexp, sigmaexp, lumi, alpha=.05, toys=200000 ):
    """ computes the 95% CL upper limit on the production cross section """
    #ret = _upperLimitMadAnalysis ( Nobs, Nexp, sigmaexp, 1.-alpha, toys ) / lumi
    #print "ret=",ret
    computer = UpperLimitComputer ( toys, lumi, 1.-alpha )
    ret = computer.ulSigmaTimesEpsilon ( Nobs, Nexp, sigmaexp, 5. )
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

    def root_func ( self, sig ):
        """ The function whose root is to be found (CLs - 0.95 ) """
        return self.CLs ( sig ) - self.cl

    """
    def root_funcMV( self, mu ):
        # print ( "MV try sig=",sig )
        sig = [ mu*x for x in self.sig ]
        cls = self.CLsMV ( sig ) - self.cl
        return cls
    """

    def ulSigma ( self, nev, xbg, cov, eff ):
        """ upper limit obtained from combined efficiencies
        :param nev: number of observed events per dataset
        :param xbg: expected bg per dataset
        :param cov: uncertainty in background, as a covariance matrix
        :param eff: dataset effective signal efficiencies
        :returns: upper limit on *production* xsec (efficiencies unfolded)
        """
        #if type(sig[0] ) == type(fb):
        #    sig = [ float(x.asUnit(fb) * self.lumi) for x in sig ]
        effs = numpy.array ( eff )
        # sigs = numpy.array ( sig )
        llhds={}
        upto = 5 * max ( nev + xbg ) / min(eff)
        dx = upto/20.
        start = dx/2.
        computer = LikelihoodComputer ( nev, xbg, cov )
        while True:
            for sig in numpy.arange ( start, upto, dx ):
                csig = sig * effs
                #lold = LLHD ( nev, xbg, cov, csig, 10*self.origNToys )
                #print ( "l=",lold,"ntoys=",self.origNToys )
                l = computer.likelihood ( csig )
                # print ( "l2=",l )
                llhds[float(sig)]=l
            norm = sum ( llhds.values() )
            last = llhds[sig]/norm
            if last < .0007:
                break ## ok, we sampled well enough?
            if last < .00001:
                ## dubious! we may have sampled to scarcely!
                logger.error ( "when integrating pdf, last bin is suspiciously small %f" % last )

            start = upto + dx/2.
            upto = 2*upto
            logger.error ( "last=%f. need to extend to %f" % ( last, upto ) )
        for k,v in llhds.items():
            llhds[k]=v/norm
        keys = llhds.keys()
        keys.sort()
        cdf = 0.
        for k in keys:
            v = llhds[k]
            cdf += v
            # print ( "k=",k,"v=",v,"cdf=",cdf )
            if cdf > .95: # perform a simple linear interpolation over cdf
                f = ( cdf - .95 ) / v
                # ret = ( k - f * dx )
                ret = ( k + dx * ( .5 - f ) )
                # print ( "ret=",ret )
                # return (k + dx * ( 1 - f )) / self.lumi
                return ret / self.lumi

    def CLs( self, SigHypothesis ):
        """ this method has been taken from MadAnalysis5, see
            Official website: <https://launchpad.net/madanalysis5>
            It is for the 1d case only .

            Thanks to the MadAnalysis for granting us the permission to use
            the code here! """
        NumObserved = self.nev
        ExpectedBG = self.xbg
        BGError = self.sbg
        NumToyExperiments = self.currentNToys

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


    def ulSigmaTimesEpsilon ( self, nev, xbg, sbg, upto=5.0, return_nan=False ):
        """ upper limit obtained via mad analysis 5 code, given for sigma
        :param nev: number of observed events
        :param xbg: expected bg
        :param sbg: uncertainty in background
        :param upto: defines the interval to be probed
        :returns: upper limit on xsec * epsilon (i.e. efficiencies *not*
                  unfolded)
        """
        self.nev = nev
        self.xbg = xbg
        self.sbg = sbg
        self.currentNToys = self.origNToys

        try:
            ## up = upto ## 5. ##  upto * max(nev,xbg,sbg)
            dn = max( 0, nev - xbg )
            up = upto * max( dn + math.sqrt(nev), dn + math.sqrt(xbg),sbg)
            return optimize.brentq ( self.root_func, 0, up, rtol=1e-3 ) / self.lumi
        except (ValueError,RuntimeError) as e:
            #print "exception: >>",type(e),e
            if not return_nan:
                # print "compute again, upto=",upto
                # self.origNToys = 5* self.origNToys
                return self.ulSigmaTimesEpsilon ( nev, xbg, sbg, 5.0*upto, upto>200. )
            else:
                return float("nan")
                
    def CLsMV( self, SigHypothesis ):
        """ CLs, with vectors and covariances """
        ## testing whether scipy is there
        NumObserved = self.nev
        ExpectedBG = self.xbg
        BGError = self.sbg
        NumToyExperiments = self.currentNToys
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
        ExpectedBGs = [ value for value in ExpectedBGs if (value > 0).all() ]

        ## complain if too many negatives
        p = float(len(ExpectedBGs )) / NumToyExperiments
        if p < 0.9:
            logger.warning ( "only %d%s of points are all-positive" % (100*p,"%"))

        ToyBGs = list ( map ( scipy.stats.poisson.rvs, ExpectedBGs ) )

        # The probability for the background alone to fluctutate as LOW as
        # observed = the fraction of the toy experiments with backgrounds as low as
        # observed = p_b.
        # NB (1 - this p_b) corresponds to what is usually called p_b for CLs.
        #print ( "ToyBGs=", ToyBGs[:3] )
        #print ( "NumObs=", NumObserved )
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

class LikelihoodComputer:
    def __init__ ( self, nobs, nb, covb ):
        """
        :param nobs: numbers of observed events (float or 1d array)
        :param nb: predicted backgrounds (float or 1d array)
        :param covb: covariance matrix of backgrounds (float or 2d array)
        """
        self.nobs = self.convert ( nobs )
        self.nb = self.convert ( nb )
        self.covb = self.convert ( covb )

    def convert ( self, x ):
        """ turn everything into numpy arrays """
        if type ( x ) in ( list, tuple ):
            return array ( x )
        return x

    #Define integrand (gaussian_(bg+signal)*poisson(nobs)):
    # def prob(x0, x1 )
    def probMV( self, *xar ):
        x = numpy.array ( xar )
        poisson = numpy.exp(self.nobs*numpy.log(x) - x - special.gammaln(self.nobs + 1))
        gaussian = stats.multivariate_normal.pdf(x,mean=self.nb+self.nsig,cov=(self.covb+numpy.diag(self.deltas**2)))
        ret = gaussian * ( reduce(lambda x, y: x*y, poisson) )
        return ret

    #Define integrand (gaussian_(bg+signal)*poisson(nobs)):
    def prob( self, x ):
        poisson = exp(self.nobs*log(x) - x - math.lgamma(self.nobs + 1))
        gaussian = stats.norm.pdf( x, loc=self.nb+self.nsig,\
                                   scale=sqrt(self.covb+self.deltas**2))
        return poisson*gaussian

    def findMax ( self ):
            #Compute maximum value for the integrand:
            sigma2 = self.covb + numpy.diag ( self.deltas**2 )
            dsigma2 = numpy.diag ( sigma2 )
            #print ( "sigma2=", sigma2 )
            #print ( "dsigma2=", dsigma2 )
            ## for now deal with variances only
            ntot = self.nb + self.nsig
            xm = self.nb + self.nsig - dsigma2
            #If nb + nsig = sigma2, shift the values slightly:
            for ctr,i in enumerate(xm):
                if i == 0.:
                    xm[ctr]=1e-3
            """
            print ( "self.nobs=",self.nobs )
            print ( "self.nb=",self.nb )
            print ( "nsig=", self.nsig )
            print ( "ntot=", ntot )
            print ( "sigma2=", sigma2 )
            """
            weight = ( numpy.matrix (sigma2 ) )**(-1) ## weight matrix
            q = self.nobs / numpy.diag ( weight ) ## q_i= nobs_i * w_ii^-1
            p = ntot - 1. / numpy.diag ( weight )
            xmax1 = p/2. * ( 1 + sqrt ( 1. + 4*q / p**2 ) ) ## no cov iteration
            #print ( "xmax1 w/ cov,p(xmax)=", xmax1, self.probMV ( *xmax1 ) )
            ndims = len(p)
            #print ( "ndims=", ndims )
            for i in range(ndims):
                for j in range(ndims):
                    if i==j: 
                        continue ## treat covariance terms
                    p[i]+=(ntot[j]-xmax1[j])*weight[i,j] / weight[i,i]
            xmax = p/2. * ( 1 + sqrt ( 1. + 4*q / p**2 ) )
            #print ( "xmax2 w/ cov,p(xmax2)=", xmax, self.probMV ( *xmax ) )
            return xmax

    def findMaxNoCov ( self ):
            """ find maximum in gauss*poisson function, disregarding
                off-diagonal elements in covariance matrix """
            sigma2 = self.covb + numpy.diag ( self.deltas**2 )
            dsigma2 = numpy.diag ( sigma2 )
            xm = self.nb + self.nsig - dsigma2
            xmax = xm*(1.+sign(xm)*sqrt(1. + 4.*self.nobs*dsigma2/xm**2))/2.
            # print ( "xmax no cov,p(xmax)=", xmax, self.probMV ( *xmax ) )
            return xmax

    def _mvLikelihood( self, nsig, deltas ):
            #     Why not a simple poisson function for the factorial
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

            #Compute maximum value for the integrand:
            sigma2 = self.covb + numpy.diag ( deltas**2 )
            dsigma2 = numpy.diag ( sigma2 )
            #print ( "sigma2=", sigma2 )
            #print ( "dsigma2=", dsigma2 )
            ## for now deal with variances only
            ntot = self.nb + nsig
            xm = self.nb + nsig - dsigma2
            #If nb + nsig = sigma2, shift the values slightly:
            for ctr,i in enumerate(xm):
                if i == 0.:
                    xm[ctr]=1e-3
            self.nsig, self.deltas = nsig, deltas ## store for integration
            xmax = self.findMaxNoCov ()
            # xmax = xm*(1.+sign(xm)*sqrt(1. + 4.*self.nobs*dsigma2/xm**2))/2.
            #print ( "xmax=", xmax )
            #print ( "xmax2=", self.findMax ( ) )
            #print ( "xmax,p(xmax)=", xmax, self.probMV ( xmax ) )

            #Define initial integration range:
            nrange = 5.
            #print ( "xmax=",xmax)
            #import IPython
            #IPython.embed()
            low = xmax-nrange*sqrt(dsigma2)
            low[low<0.] = 0. ## set all negatives to zero!
            #print ( "low=",low )
            # a = max(0.,xmax-nrange*sqrt(sigma2))
            up = xmax+nrange*sqrt(dsigma2)
            #print ( "up=",up )
            #print ( "low,up,mv",low[0],up[0],self.probMV(xmax) )
            #a = numpy.array ( [1.]*len(nobs) ) # FIXME wrong
            #b = numpy.array ( [2.]*len(nobs) ) # FIXME wrong
            self.nsig, self.deltas = nsig, deltas ## store for integration
            like = integrate.nquad( self.probMV, zip(low,up) )[0] ## fixme so wrong
            #                              epsabs=0.,epsrel=1e-3)[0]
            #print ( "like=", like )

            norm = (1./2.)*(1. + special.erf((self.nb+nsig)/sqrt(2.*dsigma2)))#[0][0]
            #norm=1.
            #print ( "norm=", norm )
            #Increase integration range until integral converges
            err = 1.
            while err > 0.01: # wrong
                like_old = like
                nrange = nrange*2
                low = xmax-nrange*sqrt(dsigma2)
                low[low<0.] = 0. ## set all negatives to zero!
                up = xmax+nrange*sqrt(dsigma2)
                #print ( "integrating from",low,"to",up )
                #a = numpy.array ( [0.]*len(nobs) ) ## FIXME wrong
                #b = numpy.array ( [4.]*len(nobs) ) ## FIXME wrong
                like = integrate.nquad(self.probMV,zip(low,up))[0]
 #                                      epsabs=0.,epsrel=1e-3)[0] FIXME so wrong
                if like == 0.:
                    err = 1.
                    continue
                err = abs(like_old-like)/like
                #print ( "like_old,like,err=",like_old/norm,like/norm,err )

            #print ( "mvLikelihood, like=", like )
            #Renormalize the likelihood to account for the cut at x = 0.
            #The integral of the gaussian from 0 to infinity gives:
            #(1/2)*(1 + Erf(mu/sqrt(2*sigma2))), so we need to divide by it
            #(for mu - sigma >> 0, the normalization gives 1.)
            norm = (1./2.)*(1. + special.erf((self.nb+nsig)/sqrt(2.*dsigma2)))#[0][0]
            # norm = 1.
            #print ( "like=",like)
            like = like/ reduce(lambda x, y: x*y, norm ) 
            #print ( "like=",like)
            return like

    def likelihood ( self, nsig, deltas = None ):
        """ compute likelihood for nsig, marginalized the nuisances """
        nsig = self.convert ( nsig )
        if type(deltas) == type(None):
            deltas = 0.2*nsig
        if type( nsig ) in [ int, float, numpy.float64 ]:
            return self._likelihood1d( nsig, deltas )
        return self._mvLikelihood( nsig, deltas )

    def _likelihood1d( self, nsig, deltas ):
            #     Why not a simple poisson function for the factorial
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

            #Compute maximum value for the integrand:
            sigma2 = self.covb + deltas**2
            xm = self.nb + nsig - sigma2
            #If nb + nsig = sigma2, shift the values slightly:
            if xm == 0.:
                xm = 0.001
            xmax = xm*(1.+sign(xm)*sqrt(1. + 4.*self.nobs*sigma2/xm**2))/2.

            #Define initial integration range:
            nrange = 5.
            a = max(0.,xmax-nrange*sqrt(sigma2))
            b = xmax+nrange*sqrt(sigma2)
            #print ( "1d,a,b,",a,b )
            self.nsig, self.deltas = nsig, deltas ## store for integral
            like = integrate.quad(self.prob,a,b,epsabs=0.,epsrel=1e-3)[0]
            #print ( "1d,like",like )

            #Increase integration range until integral converges
            err = 1.
            while err > 0.01:
                like_old = like
                nrange = nrange*2
                a = max(0.,xmax-nrange*sqrt(sigma2))
                b = xmax+nrange*sqrt(sigma2)
                like = integrate.quad(self.prob,a,b, epsabs=0.,epsrel=1e-3)[0]
                err = abs(like_old-like)/like

            #Renormalize the likelihood to account for the cut at x = 0.
            #The integral of the gaussian from 0 to infinity gives:
            #(1/2)*(1 + Erf(mu/sqrt(2*sigma2))), so we need to divide by it
            #(for mu - sigma >> 0, the normalization gives 1.)
            norm = (1./2.)*(1. + special.erf((self.nb+nsig)/sqrt(2.*sigma2)))
            # print ( "norm=", norm )
            like = like/norm

            return like

    def chi2( self, nsig, deltas=None ):
            """
            Computes the chi2 for a given number of observed events nobs
            given the predicted background nb, error on this background deltab,
            expected number of signal events nsig and, if given, the error on signal (deltas).
            :return: chi2 (float)

            """
            nsig = self.convert ( nsig )
            if deltas == None:
                deltas = 0.2 * nsig
            # Compute the likelhood for the null hypothesis (signal hypothesis) H0:
            llhd = self.likelihood( nsig, deltas )

            #print ( "deltas=",deltas )
            #print ( "nsig=",nsig )
            #Percentual signal error:
            deltas_pct = deltas / nsig ## float(nsig)

            # Compute the maximum likelihood H1, which sits at nsig = nobs - nb
            # (keeping the same % error on signal):
            dn = self.nobs-self.nb
            deltas = deltas_pct*( dn )
            maxllhd = self.likelihood( dn, deltas )

            # Return infinite likelihood if it is zero
            # This can happen in case e.g. nb >> nobs
            if llhd == 0.:
                return float('inf')

            chi2=-2*log(llhd/maxllhd)

            # Return the test statistic -2log(H0/H1)
            return chi2

if __name__ == "__main__":
    computer = UpperLimitComputer ( 1000, 1. / fb, .95 )
    nsig_,nobs_,nb_,deltab_,deltas_,covb_=1,4,3.6,.1,None,.08**2
    print ( computer.ulSigma ( [nobs_,nobs_], [nb_,nb_], [[deltab_**2,covb_],[covb_,deltab_**2]], [.1,.1] ) )
    print ( computer.ulSigma ( [nobs_,nobs_], [nb_,nb_], [[deltab_**2,covb_],[covb_,deltab_**2]], [.01,.01] ) )
    # print ( computer.ulSigma ( [4,4,4], [3.6,3.6,3.6], [[0.1**2,0,0],[0.,0.1**2,0],[0.,0,0.1**2]], [0.02,.02,.02] ) )

    # print ( "LLHD", LLHD ( [nobs_], [nb_], [[deltab_**2]], [ nsig_ ], 1000 ) )
    # print ( "LLHD", LLHD ( [nobs_,nobs_], [nb_,nb_], [[deltab_**2,0.**2],[0.**2,deltab_**2]], [ nobs_, nobs_ ], 100 ) )
    computer = LikelihoodComputer ( nobs_, nb_, deltab_**2 )
    print ( "1d, computer:", computer.likelihood( nsig_, deltas_ )  )
    #print ( "1d, chi2:",computer.chi2 ( nsig_ ) )
    computer = LikelihoodComputer ( [nobs_], [nb_], numpy.diag([deltab_**2]) )
    print ( "mv 1d, computer:",computer.likelihood( [nsig_], deltas_) )
    print ( "mv 1d, chi2:",computer.chi2 ( [nsig_] ) )
    #sys.exit()
    cov = numpy.diag ([deltab_**2,deltab_**2])
    cov[0,1] = .01
    cov[1,0] = .01
    computer = LikelihoodComputer ( [nobs_,nobs_], [nb_,nb_], cov )
    l = computer.likelihood ( array([nsig_,nsig_]), deltas_  )
    print ( "mv 2d, computer:", l )
    print ( "mv 2d, chi2:", computer.chi2 ( [nsig_,nsig_] ) )
    # print ( likelihoodMV(array([nsig]), array([nobs]), array([nb]), numpy.diag([deltab**2]), deltas) )
    # print ( likelihoodMV(array([nsig,nsig]), array([nobs,nobs]), array([nb,nb]), numpy.diag([deltab**2,deltab**2]), deltas) )

    dummy_nobs = [ 1964, 877, 354, 182, 82, 36, 15, 11 ]
    dummy_nbg = [ 2006.4, 836.4, 350., 147.1, 62., 26.2, 11.1, 4.7 ]
    dummy_si = [ 47., 29.4, 21.1, 14.3, 9.4, 7.1, 4.7, 4.3 ]
    dummy_cov = [ [ 16787.2, -1691.3, -4520.3, -3599.9, -2286.4, -1316.5, -719.8, -381.1 ] ,
                  [ -1691.3, 603.1, 754.6, 513.3, 294., 154.9, 78.1, 38.3 ],
                  [ -4520.3, 754.6, 1454., 1110.9, 691.1, 392.3, 212.1, 111.2 ],
                  [ -3599.9, 513.3, 1110.9, 871.2, 551.8, 318.1, 174.3, 92.5 ],
                  [ -2286.4, 294., 691.1, 551.8, 353.9, 206.2, 114.1, 61. ],
                  [ -1316.5, 154.9, 392.3, 318.1, 206.2, 121.3, 67.6, 36.4 ],
                  [ -719.8, 78.1, 212.1, 174.3, 114.1, 67.6, 38.0, 20.6 ],
                  [ -381.1, 38.3, 111.2, 92.5, 61.0, 36.4, 20.6, 11.2 ],
    ]
    # print ( LLHD ( dummy_nobs, dummy_nbg, dummy_cov,dummy_si, 100 ) )
    # print ( computer.ulSigma ( dummy_nobs, dummy_nbg, dummy_cov,dummy_si ) )
