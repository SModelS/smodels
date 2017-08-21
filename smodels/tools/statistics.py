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
from functools import reduce
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

    def ulSigma ( self, nev, xbg, cov, eff ):
        """ upper limit obtained from combined efficiencies
        :param nev: number of observed events per dataset
        :param xbg: expected bg per dataset
        :param cov: uncertainty in background, as a covariance matrix
        :param eff: dataset effective signal efficiencies
        :returns: upper limit on *production* xsec (efficiencies unfolded)
        """
        computer = LikelihoodComputer ( nev, xbg, cov )
        effs = numpy.array ( eff )
        llhds={}
        mu_hat = computer.findMuHat ( eff )
        upto = mu_hat + 5.*computer.getSigmaMu ( eff, 1., mu_hat )
        first_upto = upto
        n_bins = 100.
        dx = upto / n_bins ## FIXME
        start = dx/2.
        # start = 0.
        #print ( "compute ulSigma",nev,xbg,cov,eff )
        #print ( "start,upto=",start,upto )
        while True:
            lst = [] ## also store in list
            for sig in numpy.arange ( start, upto, dx ):
                csig = sig * effs
                # l = computer.likelihood ( csig )
                # FIXME marginalize or profile?
                l = computer.profileLikelihood ( csig )
                llhds[float(sig)]=l
                lst.append ( l )

            ##print ( "likelihoods",llhds )
            norm = sum ( llhds.values() )
            #print ( "norm",norm )
            last = lst[-1]/norm
            # print ( "maximum at bin #", lst.index ( max ( lst ) ) )
            if lst.index ( max ( lst ) ) == 0 and lst[0]/lst[1] > 1.1:
                upto = .24 * upto
                logger.error ( "maximum is too close to the left (first two bins are %g, %g). lets go for smaller range: %f" % ( lst[0], lst[1], upto ) )
                dx = upto / n_bins
                start = dx/2.
                llhds={}
                continue ## and again
            if last < 1e-8:
                ## dubious! we may have sampled too scarcely!
                upto = .23 * upto
                logger.error ( "when integrating pdf, last bin is suspiciously small: %g. Take 23 pc the range: %s." % ( last, upto ) )
                dx = upto / n_bins ## FIXME
                llhds={}
                continue
            else:
                if last < .0005: ## last bin contributes > 1e-7 and < .0005?
                    break ## ok, we sampled well enough.
            if last > 1e-2:
                ## the last bin contributes so much, better start from scratch
                logger.error ( "very large last bin: %g. extend by factor of 3.5." % last )
                llhds={}
                upto = 3.5 * upto
                dx = upto / n_bins
                continue
            start = upto + dx/2.
            upto = 2*upto
            logger.error ( "last bin contributes %f to integral. need to extend to %f" % ( last, upto ) )
        for k,v in llhds.items():
            llhds[k]=v/norm
        ## now find the 95% quantile by interpolation
        ret = self.interpolate ( llhds, dx )
        self.plot ( llhds, dx, first_upto, upto, ret, computer, effs )
        # computer.plotLTheta ( )
        return ret

    def plot ( self, llhds, dx, upto, final_upto, cl95, computer, effs ):
        """ function to plot likelihoods and the 95% CL, for
            debugging purposes.
        :param llhds: the likelihood values, as a dictionary
        :param dx: the delta x between the likelihood values. FIXME redundant.
        :param upto: first guess for how far we need to integrate.
        :param final_upto: last guess for how far we need to integrate.
        :param cl95: final 95% CL upper limit.
        :param xmax: maximum likelihood computer
        """
        expected = ( computer.nobs == computer.nb ).all()
        # logger.error ( "expected=%s", expected )
        import ROOT
        xvals = list ( llhds.keys() )
        xvals.sort()
        yvals = []
        for x in xvals:
            yvals.append ( llhds[x] )
        lumi = self.lumi.asNumber(1/fb)
        logger.error ( "Plotting likelihoods: int.upto: %s, int.upto1: %s, cl95: %s" % \
                       ( upto/lumi, final_upto/lumi, cl95 ) )
        t=ROOT.TGraph ( len ( xvals ))
        l=ROOT.TLegend( .6,.7,.98,.89)
        for ctr,(x,y) in enumerate ( zip ( xvals, yvals ) ):
            xv = x / lumi
            t.SetPoint ( ctr, xv, y )
            # logger.error ( "x,y=%s,%s" % ( x,y ) )
        t.Draw("AC*")
        # logger.error ( "integral %f" % t.Integral() )
        t.GetXaxis().SetTitle ( "cross section [fb]" )
        title = " likelihood for signal cross section (%d SRs)" % len(effs)
        if expected:
            title = "expected" + title
        else:
            title = "observed" + title
        t.SetTitle ( "" )
        tt = ROOT.TText ( .1, .92, title )
        tt.SetTextSize(.04)
        tt.SetNDC()
        tt.Draw()
        t3=ROOT.TLine ( cl95.asNumber(fb), 0., cl95.asNumber(fb), max(yvals) )
        t3.SetLineStyle(3)
        t3.SetLineColor(ROOT.kRed)
        t3.SetLineWidth(2)
        t3.Draw("SAME")
        l.AddEntry(t3,"95% limit","L" )
        t4=ROOT.TLine ( upto/lumi, 0., upto/lumi, max(yvals) )
        t4.SetLineColor(ROOT.kBlue)
        t4.SetLineStyle(4)
        t4.SetLineWidth(2)
        t4.Draw("SAME")
        l.AddEntry(t4,"first guess for integration upper limit","L" )
        t5=ROOT.TLine ( final_upto/lumi, 0., final_upto/lumi, max(yvals) )
        t5.SetLineColor(ROOT.kGreen )
        t5.SetLineStyle(5)
        t5.SetLineWidth(2)
        t5.Draw("SAME")
        l.AddEntry(t5,"final upper int limit","L" )
        sigma_max = computer.findMuHat( effs, lumi )
        t6=ROOT.TLine ( sigma_max, 0., sigma_max, max(yvals) )
        t6.SetLineColor(ROOT.kCyan )
        t6.SetLineStyle(6)
        t6.SetLineWidth(2)
        t6.Draw("SAME")
        l.AddEntry(t6,"max value","L" )
        maxp1 = sigma_max + 1.96 * computer.getSigmaMu ( effs, lumi, sigma_max )
        t7=ROOT.TLine ( maxp1, 0., maxp1, max(yvals) )
        t7.SetLineColor(ROOT.kMagenta )
        t7.SetLineStyle(7)
        t7.SetLineWidth(2)
        t7.Draw("SAME")
        l.AddEntry(t7,"max + 1.96*sigma","L" )
        # logger.error ( "95 pc ul=%s" % cl95 )
        l.Draw()
        fname = "plot%d.png" % expected
        ROOT.c1.Print ( fname )

    def interpolate ( self, llhds, dx ):
        """ interpolate likelihoods to have the best possible
            estimate for the 95% CL
        :param llhds: dictionary of (normalized) likelihoods
        :param dx: step size in sigma_sig
        :returns: 95% CL
        """
        keys = list ( llhds.keys() )
        keys.sort()
        cdf = 0.
        quant = .95 ## 95% quantile
        for ctr,k in enumerate(keys):
            v = llhds[k]
            cdf += v
            if cdf > quant: # perform a simple linear interpolation over cdf
                if ctr < 0.1 * len(keys):
                    logger.error ( "for the 95% CL we needed fewer than 10% of points. bad." )
                    sys.exit()
                if ctr > 0.9 * len(keys):
                    logger.error ( "for the 95% CL we needed more than 90% of points. bad." )
                    sys.exit()
                f = ( cdf - quant ) / v
                # ret = ( k - f * dx )
                ret = ( k + dx * ( .5 - f ) )
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

    def dLdMu ( self, mu, effs, theta_hat ):
        """ dL/dmu * ( 1/L ), if L is the likelihood. The function
            whose root gives us muhat, i.e. the mu that maximizes
            the likelihood. """
        denominator = mu*effs + self.nb + theta_hat
        ret = self.nobs * effs / denominator - effs
        if type (ret ) in [ numpy.array, numpy.ndarray, list ]:
            ret = sum ( ret )
        return ret

    def plotMuHatRootFinding ( self, effs, theta_hat, lumi, mu_hat ):
        """ plot the function whose root gives us mu_hat: 1/L dL/dmu.
        """
        mu_c = ( self.nobs - self.nb - theta_hat ) / effs
        mu_ini = sum ( mu_c ) / len(mu_c)

        rnge = numpy.arange ( -1.*mu_ini, 4.*mu_ini, .5*mu_ini )
        import ROOT
        r=ROOT.TGraph ( len(rnge) )
        r.SetTitle ( "1/L dL/d#mu" ) ##  #frac{1}{L}#frac{dL}#frac{d#mu}" )
        mu_c = ( self.nobs - self.nb - theta_hat ) / effs
        logger.error ( "mu_c = %s " % (mu_c/lumi) )
        for ctr,i in enumerate(rnge):
            dldmu = self.dLdMu ( i, effs, theta_hat )
            # logger.error ( " mu=%s, p=%s" % ( i, dldmu ) )
            r.SetPoint ( ctr, i / lumi, dldmu )
        r.Draw("AC*")
        r.GetHistogram().SetXTitle ( "#mu [fb]" )
        r2=ROOT.TGraph ( 1 )
        r2.SetPoint(0,mu_hat/lumi, 0. )
        r2.SetMarkerStyle(29)
        r2.SetMarkerColor(ROOT.kRed)
        r2.SetMarkerSize(2)
        r2.Draw("PSAME" )
        ROOT.c1.Print ( "dldmu.pdf" )

    def findMuHat ( self, effs, lumi=1. ):
        """
        find the most likely signal strength mu
        FIXME what about the profiled nuisances?
        :param lumi: return yield (lumi=1.) or cross section ("real" lumi)
        :returns: mu_hat, either as signal yield (lumi=1.), or as cross section.
        """
        if ( self.nb == self.nobs ).all():
            return 0. / lumi
        if type ( effs ) in [ list, numpy.ndarray ]:
            effs = numpy.array ( effs )
        self.nsig = self.nobs - self.nb
        self.deltas = numpy.array ( [0.] * len(self.nobs ) )
        theta_hat = self.findThetaHat()
        logger.error ( "   theta hat=%s" % theta_hat[:11] )
        logger.error ( "        nobs=%s" % self.nobs[:11] )
        logger.error ( "          nb=%s" % self.nb[:11] )
        logger.error ( "        effs=%s" % effs[:11] )
        mu_c = ( self.nobs - self.nb - theta_hat ) / effs
        ## find mu_hat by finding the root of 1/L dL/dmu. We know 
        ## that the zero has to be between min(mu_c) and max(mu_c).
        lower,upper = min(mu_c),max(mu_c)
        try:
            mu_hat = optimize.brentq ( self.dLdMu, lower, upper, args=(effs, theta_hat ) )
        except ValueError as e:
            logger.error ( "could not find root of dLdMu! Values are [%f=%f,%f=%f]. This should not happen, but for now we proceed with larger range." % ( lower,self.dLdMu(lower,effs,theta_hat),upper,self.dLdMu(upper,effs,theta_hat)) )
            mu_hat = optimize.brentq ( self.dLdMu, 0., 5*upper, args=(effs, theta_hat ) )
        # self.plotMuHatRootFinding( effs, theta_hat, lumi, mu_hat ) 
        logger.error ( "   mu_hat  =%s" % (mu_hat/lumi) )
        check = sum ( self.nobs * effs / ( mu_hat * effs + self.nb ) - effs )
        logger.error ( "   check(mu_hat) =%s" % check )
        return mu_hat / lumi

    def getSigmaMu ( self, effs, lumi, mu_hat ):
        """
        get a rough estimate for the variance of mu around mu_max
        FIXME need to do this more thorougly.
        """
        s_effs = numpy.array ( effs )
        # ret = math.sqrt ( sum ( ( mu_hat * lumi * s_effs + self.nb )**2 / (self.nobs * s_effs**2 ) ) ) / lumi
        # logger.info ( "cov = %s" % (ret) )
        if type ( effs ) in [ list, numpy.ndarray ]:
            s_effs = sum ( effs )
        sgm_mu = math.sqrt ( sum (self.nobs ) + sum ( numpy.diag ( self.covb ) ) ) / s_effs / lumi
        logger.info ( "sgm_mu = %s" % sgm_mu )
        return sgm_mu


    def convert ( self, x ):
        """ turn everything into numpy arrays """
        if type ( x ) in ( list, tuple ):
            return array ( x )
        return x

    #Define integrand (gaussian_(bg+signal)*poisson(nobs)):
    # def prob(x0, x1 )
    def probMV( self, *xar ):
        x = numpy.array ( xar )
        #if sum ( x<0. ) > 0:
        #    logger.error ( "negative values for x: %s" % x )
        #    sys.exit()
        #poisson = numpy.exp(self.nobs*numpy.log(x) - x - special.gammaln(self.nobs + 1))
        #gaussian = stats.multivariate_normal.pdf(x,mean=self.nb+self.nsig,cov=(self.covb+numpy.diag(self.deltas**2)))
        ntot = self.nb + self.nsig
        xtot = x+ntot
        for ctr,i in enumerate ( xtot ):
            if i==0.:
                logger.info ( "encountered a zero in probMV at position %d. Will set to small value." % ctr )
                xtot[ctr]=1e-20
        poisson = numpy.exp(self.nobs*numpy.log(xtot) - x -ntot - special.gammaln(self.nobs + 1))
        gaussian = stats.multivariate_normal.pdf(x,mean=[0.]*len(x),cov=(self.covb+numpy.diag(self.deltas**2)))
        ret = gaussian * ( reduce(lambda x, y: x*y, poisson) )
        return ret

    #Define integrand (gaussian_(bg+signal)*poisson(nobs)):
    def prob( self, x ):
        #poisson = exp(self.nobs*log(x) - x - math.lgamma(self.nobs + 1))
        #gaussian = stats.norm.pdf( x, loc=self.nb+self.nsig,\
        #                           scale=sqrt(self.covb+self.deltas**2))
        ntot = self.nb + self.nsig
        poisson = exp(self.nobs*log(x+ntot) - x - ntot - math.lgamma(self.nobs + 1))
        gaussian = stats.norm.pdf( x, loc=0.,\
                                   scale=sqrt(self.covb+self.deltas**2))
        return poisson*gaussian

    def getThetaHat ( self, nobs, nb, nsig, covb, deltas ):
            """ Compute nuisance parameter theta that
            maximizes our likelihood (poisson*gauss).  """
            sigma2 = covb + numpy.diag ( deltas**2 )
            dsigma2 = numpy.diag ( sigma2 )
            ## for now deal with variances only
            ntot = nb + nsig
            cov = numpy.matrix ( sigma2 )
            weight = ( numpy.matrix (sigma2 ) )**(-1) ## weight matrix
            diag_cov = numpy.diag(cov)
            q = nobs * diag_cov ## q_i= nobs_i * w_ii^-1
            p = ntot - diag_cov ## caveat, this is usually called '-p', not 'p'
            # we start with assuming all nuisances and covariances to be zero
            xmax = p/2. * ( 1 + sign(p) * sqrt ( 1. + 4*q / p**2 ) )
            # logger.error ( "xmax 1  = %s" % xmax )
            ndims = len(p)
            for itr in range(41): ## correct a maximum of three times
                for i in range(ndims):
                    for j in range(ndims):
                        if i==j:
                            continue ## treat covariance terms
                        p[i]=ntot[i] - diag_cov[i] + (ntot[j]-xmax[j])*weight[i,j] / weight[i,i]
                        # p[i]+=(ntot[j]-xmax[j])*weight[i,j] / weight[i,i]
                new_xmax = p/2. * ( 1 + sign(p)* sqrt ( 1. + 4*q / p**2 ) )
                for ctr,i in enumerate(xmax):
                    if i==0.: ## watch out for zeroes.
                        xmax[ctr]=1e-20
                # logger.error ( "xmax=%s" % xmax )
                change = sum ( abs ( new_xmax - xmax ) / xmax ) / len(xmax)
                # logger.error ( "change=%s" % change )
                if math.isnan ( change ):
                    sys.exit()
                if change < 1e-3:
                    ## smaller 0.1% change? stop!
                    xmax = new_xmax
                    break
                if itr == 5:
                    logger.error ( "profiling did not converge in %d steps: %s vs %s. change=%f " % (itr, xmax[:3], new_xmax[:3],change ) )
                    xmax = new_xmax
                    break
                xmax = new_xmax
                # logger.error ( "xmax %d = %s" % (itr+2, xmax ) )
            return xmax - nb - nsig ## FIXME wrong!!!
            # return xmax

    def findThetaHat ( self ):
            """ Compute nuisance parameter theta that maximizes our likelihood
                (poisson*gauss).
            """
            return self.getThetaHat ( self.nobs, self.nb, self.nsig, self.covb, \
                                      self.deltas )

    def findThetaHatNoCov ( self ):
            """ find nuisance parameter theta that maximizes our likelihood
                (gauss*poisson) function, disregarding
                off-diagonal elements in covariance matrix """
            sigma2 = self.covb + numpy.diag ( self.deltas**2 )
            dsigma2 = numpy.diag ( sigma2 )
            xm = self.nb + self.nsig - dsigma2
            # print ( "[findMaxNoCov] xm=", xm )
            xmax = xm*(1.+sign(xm)*sqrt(1. + 4.*self.nobs*dsigma2/xm**2))/2.
            # print ( "xmax no cov,p(xmax)=", xmax, self.probMV ( *xmax ) )
            return xmax

    def _mvProfileLikelihood( self, nsig, deltas ):
            # compute the profiled (not normalized) likelihood of observing
            # nsig signal events
            self.nsig, self.deltas = nsig, deltas ## store for integration
            theta_hat = self.findThetaHat ()
            return self.probMV ( *theta_hat )

    def plotLTheta ( self, nsig=None ):
        """ plot the likelihood, but as a function of the nuisance theta! """
        if ( self.nobs == self.nb ).all():
            # only plot we this isnt the 'expected' case.
            return 
        import ROOT
        oldsig = self.nsig
        theta_hat = self.findThetaHat ()
        v_tm = self.probMV ( *theta_hat )
        logger.error ( "point at %s %s" % ( theta_hat, v_tm ) )
        rng = list ( numpy.arange ( -1e10, 1e10, 1e8 ) )
        rng = list ( numpy.arange ( -1.2e9, 1.2e9, 0.8e7 ) )
        t=ROOT.TGraph(len(rng))
        t.SetTitle ( "L(#theta/#hat{#theta}), %d SRs" % len(theta_hat) )
        NLLs = []
        llhds=[]
        for ctr,i in enumerate(rng):
            v=i*theta_hat
            L = self.probMV ( *v )
            llhds.append ( L )
            t.SetPoint( ctr, i, L )
            if L > 0.:
                nll = - math.log ( L )
                NLLs.append ( nll )
                t.SetPoint( ctr, i, L )
            else:
                t.SetPoint( ctr, i, 450. ) ## float("nan") )
        mNLL = min(NLLs)
        for ctr,L in enumerate(llhds):
            if L == 0.:
                pass
                # t.SetPoint ( ctr, rng[ctr], 1.1*max(NLLs) )
        # t.SetPoint(0,theta_hat, - math.log ( self.probMV ( *theta_hat ) ) )
        t.Draw("AC*")
        t.GetHistogram().SetXTitle ( "#theta/#hat{#theta}" )
        # line = ROOT.TLine ( 1., min(NLLs), 1., max(NLLs) )
        line = ROOT.TLine ( 1., min(llhds), 1., max(llhds) )
        line.SetLineColor ( ROOT.kRed )
        line.SetLineStyle(2)
        line.Draw()
        ROOT.c1.Print("Ltheta.png" )
        self.nsig = oldsig

    def _mvLikelihood( self, nsig, deltas ):
            # compute the marginalized likelihood of observing nsig signal
            # events


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
            # xmax = self.findMaxNoCov ()
            theta_hat = self.findThetaHat ()
            #print ( "found xmax at", xmax, "l=",self.probMV(*xmax) )

            #Define initial integration range:
            nrange = 5.
            low = theta_hat-nrange*sqrt(dsigma2)
            low[low<0.] = 0. ## set all negatives to zero!
            up = theta_hat+nrange*sqrt(dsigma2)
            ## first integral can be coarse, we will anyhow do another round
            opts = { "epsabs":1e-2, "epsrel":1e-1 }
            like = integrate.nquad( self.probMV, zip(low,up),opts=opts )[0]
            norm = (1./2.)*(1. + special.erf((self.nb+nsig)/sqrt(2.*dsigma2)))#[0][0]
            norms = reduce(lambda x, y: x*y, norm )
            err = 1.
            #print ( "found xmax2 at", xmax2, "l2=",self.probMV(*xmax2) )
            # sys.exit()
            # print ( "like=", self.probMV(*xmax)/ norms )
            while err > 0.01: # wrong
                opts = { "epsabs":1e-4, "epsrel":1e-2 }
                like_old = like
                nrange = nrange*2
                low = theta_hat-nrange*sqrt(dsigma2)
                low[low<0.] = 0. ## set all negatives to zero!
                up = theta_hat+nrange*sqrt(dsigma2)
                res = integrate.nquad(self.probMV,zip(low,up),opts=opts)
                like = res[0]
                if like == 0.: ## too large!
                    return 0.
                err = abs(like_old-like)/like
                ## print ( "res=",res[1]/res[0],"err=",err,"nrange=",nrange )

            #Renormalize the likelihood to account for the cut at x = 0.
            #The integral of the gaussian from 0 to infinity gives:
            #(1/2)*(1 + Erf(mu/sqrt(2*sigma2))), so we need to divide by it
            #(for mu - sigma >> 0, the normalization gives 1.)
            like = like/ norms ## reduce(lambda x, y: x*y, norm )
            # print ( "like now=", like )
            return like

    def profileLikelihood ( self, nsig, deltas = None ):
        """ compute the profiled likelihood for nsig.
            Warning: not normalized. """
        nsig = self.convert ( nsig )
        if type(deltas) == type(None):
            deltas = 1e-5*nsig ## FIXME for backwards compatibility
        if type( nsig ) in [ int, float, numpy.float64 ]:
            ## FIXME write profiled likelihood for 1d case.
            return self._likelihood1d( nsig, deltas )
        return self._mvProfileLikelihood( nsig, deltas )

    def likelihood ( self, nsig, deltas = None ):
        """ compute likelihood for nsig, marginalized the nuisances """
        nsig = self.convert ( nsig )
        if type(deltas) == type(None):
            deltas = 1e-5*nsig ## FIXME for backwards compatibility
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
                if like == 0.:
                    return like
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
    computer = UpperLimitComputer ( 50000, 20.5 / fb, .95 )
    nsig_,nobs_,nb_,deltab_,eff_=1,15,17.5,3.2,0.00454755
    #nsig_,nobs_,nb_,deltab_,eff_=1,15,17.5,1.0,1.
    ul = computer.ulSigmaTimesEpsilon ( nobs_, nb_, deltab_ )
    print ( "    ulSigmaTimesEpsilon=",ul)
    print ( "                 ul/eff=",ul/eff_)
    uls = computer.ulSigma ( [nobs_], [nb_], [[deltab_**2]], [eff_] )
    print ( "                ulSigma=",uls)
    sys.exit()
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
