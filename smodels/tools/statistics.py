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
from smodels.tools.likelihoods import LikelihoodComputer
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
    debug_mode = False

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
        mu_hat = computer.findMuHat ( eff )
        effs = numpy.array ( eff )
        theta_hat = computer.findThetaHat ( mu_hat * effs )
        sigma_mu = computer.getSigmaMu ( effs, 1., mu_hat, theta_hat )

        llhds={}
        upto = mu_hat + 4.0 * sigma_mu
        n_bins = 50.
        dx = upto / n_bins ## FIXME
        start = dx/2.
        # start = 0.
        #print ( "compute ulSigma",nev,xbg,cov,eff )
        #print ( "start,upto=",start,upto )
        errs = []

        while True:
            lst = [] ## also store in list
            # for sig in [ 947.567985654831431930, 1579.279976091385833570, 2210.991966527939894149 ]:
            for sig in numpy.arange ( start, upto, dx ):
                csig = sig * effs
                # l = computer.likelihood ( csig )
                # FIXME marginalize or profile?
                l,err = computer.profileLikelihood ( csig )
                # logger.error ( "llhd %.18f -> %s [%s]" % ( sig, l, err ) )
                llhds[float(sig)]=l
                lst.append ( l )
                errs.append ( err )
            # sys.exit()

            ##print ( "likelihoods",llhds )
            norm = sum ( llhds.values() )
            if norm == 0.:
                logger.error ( "Zero norm!" )
            last = lst[-1]/norm
            # print ( "maximum at bin #", lst.index ( max ( lst ) ) )
            if lst.index ( max ( lst ) ) == 0 and lst[0]/lst[1] > 1.1:
                upto = .24 * upto
                logger.error ( "maximum is too close to the left (first two bins are %g, %g). lets go for smaller range: %f" % ( lst[0], lst[1], upto ) )
                dx = upto / n_bins
                start = dx/2.
                llhds={}
                continue ## and again
            if last < 1e-20:
                ## dubious! we may have sampled too scarcely!
                upto = .23 * upto
                logger.error ( "when integrating pdf, last bin is suspiciously small: %g. Take 23 pc the range: %s." % ( last, upto ) )
                # logger.error ( "here are the last 10 bins  %s" % (llhds) )
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
        self.plot ( llhds, dx, upto, ret, computer, effs, norm, errs )
        computer.plotLTheta ( nsig=mu_hat*effs )
        computer.printProfilingStats()
        return ret

    def plot ( self, llhds, dx, upto, cl95, computer, effs, norm, errs ):
        """ function to plot likelihoods and the 95% CL, for
            debugging purposes.
        :param llhds: the likelihood values, as a dictionary
        :param dx: the delta x between the likelihood values. FIXME redundant.
        :param upto: guess for how far we need to integrate.
        :param cl95: final 95% CL upper limit.
        :param xmax: maximum likelihood computer
        """
        if not self.debug_mode: return
        expected = ( computer.nobs == computer.nb ).all()
        import ROOT, time
        xvals = list ( llhds.keys() )
        xvals.sort()
        yvals = []
        for x in xvals:
            yvals.append ( llhds[x] )
        lumi = self.lumi.asNumber(1/fb)
        s_expected = "observed"
        if expected:
            s_expected = "expected"
        logger.warning ( "Plotting %s likelihoods: int.upto: %s, cl95: %s" % \
                       ( s_expected, upto/lumi, cl95 ) )
        t=ROOT.TGraph ( len ( xvals ))
        n_errs = numpy.array ( errs )
        t_err=ROOT.TGraph ( len ( n_errs[n_errs!=0] ) )
        t_err.SetMarkerColor ( ROOT.kRed )
        t_err.SetMarkerStyle ( 21 )
        l=ROOT.TLegend( .6,.7,.98,.89)
        errctr=0
        for ctr,(x,y) in enumerate ( zip ( xvals, yvals ) ):
            xv = x / lumi
            t.SetPoint ( ctr, xv, y )
            if errs[ctr]!=0:
                t_err.SetPoint ( errctr, xv, y )
                errctr+=1
            # logger.error ( "x,y=%s,%s" % ( x,y ) )
        t.Draw("AC*")
        if errctr>0:
            t_err.Draw("*SAME")
        # logger.error ( "integral %f" % t.Integral() )
        t.GetXaxis().SetTitle ( "cross section [fb]" )
        title = "%s likelihood for signal cross section (%d SRs)" % (s_expected, len(effs) )
        t.SetTitle ( "" )
        tt = ROOT.TText ( .1, .92, title )
        tt.SetTextSize(.04)
        tt.SetNDC()
        tt.Draw()
        t3=ROOT.TLine ( cl95.asNumber(fb), 0., cl95.asNumber(fb), max(yvals) )
        t3.SetLineStyle(1)
        t3.SetLineColor(ROOT.kGreen)
        t3.SetLineWidth(3)
        t3.Draw("SAME")
        l.AddEntry(t3,"95% limit","L" )
        t4=ROOT.TLine ( upto/lumi, 0., upto/lumi, max(yvals) )
        t4.SetLineColor(ROOT.kBlue)
        t4.SetLineStyle(4)
        t4.SetLineWidth(2)
        t4.Draw("SAME")
        l.AddEntry(t4,"guess for integration upper limit","L" )
        mu_hat = computer.findMuHat( effs, lumi )
        t6=ROOT.TLine ( mu_hat, 0., mu_hat, max(yvals) )
        t6.SetLineColor(ROOT.kCyan )
        t6.SetLineStyle(6)
        t6.SetLineWidth(2)
        t6.Draw("SAME")
        l.AddEntry(t6,"#hat{#mu}","L" )
        theta_hat = computer.findThetaHat ( mu_hat * effs )
        sigma_mu = computer.getSigmaMu ( effs, lumi, mu_hat, theta_hat )
        maxp1 = mu_hat + 1.96 * sigma_mu
        ts = ROOT.TText ( .01, .01, time.asctime() )
        ts.SetNDC()
        ts.SetTextSize(.03)
        ts.Draw()
        t7=ROOT.TLine ( maxp1, 0., maxp1, max(yvals) )
        t7.SetLineColor(ROOT.kMagenta )
        t7.SetLineStyle(7)
        t7.SetLineWidth(2)
        t7.Draw("SAME")
        l.AddEntry(t7,"#hat{#mu} + 1.96*#sigma","L" )
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
