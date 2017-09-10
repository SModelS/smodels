#!/usr/bin/env python

"""
.. module:: likelihoods
   :synopsis: Holds the likelihood computer

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
from smodels.tools.physicsUnits import fb
from smodels.tools.smodelsLogging import logger
from scipy import stats, optimize, integrate, special
from numpy import sqrt, exp, log, sign, array, matrix
from functools import reduce
import time
import numpy
import math
import sys

class LikelihoodComputer:

    """ the default value for delta_s, the assumed relative error
    on signal yields. """
    deltas_default = 0.2
    debug_mode = False

    def __init__ ( self, nobs, nb, covb ):
        """
        :param nobs: numbers of observed events (float or 1d array)
        :param nb: predicted backgrounds (float or 1d array)
        :param covb: covariance matrix of backgrounds (float or 2d array)
        """
        self.nobs = self.convert ( nobs )
        self.nb = self.convert ( nb )
        self.covb = self.convert ( covb )
        self.timer = { "fmin1": 0., "fmin2": 0., "theta_hat_ini": 0., "brentq": 0., 
                       "profile": 0., "find_mu_hat": 0. }

    def printProfilingStats ( self ):
        return
        print ()
        print ( "Profiling likelihood, %d SRs:" % len(self.nobs) )
        print ( "=================================" )
        for k,v in self.timer.items():
            print ( "%s: %f" % ( k, v ) )
        print ()

    def dLdMu ( self, mu, effs, theta_hat ):
        """ d ln L/dmu, if L is the likelihood. The function
            whose root gives us muhat, i.e. the mu that maximizes
            the likelihood. """
        denominator = mu*effs + self.nb + theta_hat 
        for ctr,d in enumerate ( denominator ):
            if d==0.:
                if (self.nobs[ctr]*effs[ctr]) == 0.:
                    logger.debug ( "zero denominator, but numerator also zero, so we set denom to 1." )
                    denominator[ctr]=1.
                else:
                    logger.error ( "we have a zero value in the denominator at pos %d, with a non-zero numerator. dont know how to handle." % ctr )
        ret = self.nobs * effs / denominator - effs
        if type (ret ) in [ numpy.array, numpy.ndarray, list ]:
            ret = sum ( ret )
        return ret

    def plotMuHatRootFinding ( self, effs, theta_hat, lumi, mu_hat ):
        """ plot the function whose root gives us mu_hat: 1/L dL/dmu.
        """
        if not self.debug_mode: return
        if ( self.nobs == self.nb ).all():
            return ## not for expected 
        mu_c = numpy.abs ( self.nobs - self.nb - theta_hat ) / effs
        # mu_ini = sum ( mu_c ) / len(mu_c)
        mu_ini = 2. * max ( mu_c )

        rnge = numpy.arange ( -0.*mu_ini, 2*mu_ini, .01*mu_ini )
        import ROOT
        r=ROOT.TGraph ( len(rnge) )
        r.SetTitle ( "- #partial/#partial#mu ln(L)" ) ##  #frac{1}{L}#frac{dL}#frac{d#mu}" )
        mu_c = ( self.nobs - self.nb - theta_hat ) / effs
        # logger.error ( "mu_c = %s " % (mu_c/lumi) )
        for ctr,i in enumerate(rnge):
            dldmu = self.dLdMu ( i, effs, theta_hat )
            # logger.error ( " mu=%s, p=%s" % ( i, dldmu ) )
            r.SetPoint ( ctr, i / lumi, dldmu )
        r.Draw("AC*")
        r.GetHistogram().SetXTitle ( "#mu [fb]" )
        if mu_hat != None:
            r2=ROOT.TGraph ( 1 )
            r2.SetPoint(0,mu_hat/lumi, 0. )
            r2.SetMarkerStyle(29)
            r2.SetMarkerColor(ROOT.kRed)
            r2.SetMarkerSize(2)
            r2.Draw("PSAME" )
        zero = ROOT.TLine(min(rnge),0.,max(rnge),0.)
        zero.SetLineColor ( ROOT.kGreen )
        zero.SetLineStyle(5)
        zero.Draw()
        ROOT.c1.Print ( "dldmu.pdf" )

    def findMuHat ( self, effs, lumi=1. ):
        """
        find the most likely signal strength mu
        :param lumi: return yield (lumi=1.) or cross section ("real" lumi)
        :returns: mu_hat, either as signal yield (lumi=1.), or as cross section.
        """
        if ( self.nb == self.nobs ).all():
            return 0. / lumi
        self.timer["find_mu_hat"]-=time.time()
        if type ( effs ) in [ list, numpy.ndarray ]:
            effs = numpy.array ( effs )
        effs[effs==0.]=1e-20
        if sum ( effs<0. ):
            logger.error ( "Negative efficiencies!" )
            sys.exit()
        ## we need a very rough initial guess for mu(hat), to come
        ## up with a first theta
        self.nsig = numpy.array ( [0.] * len(self.nobs ) )
        self.deltas = numpy.array ( [0.] * len(self.nobs) )
        ## we start with theta_hat being all zeroes
        theta_hat = numpy.array ( [0.] * len(self.nobs ) )
        mu_hat_old, mu_hat = 0., 1.
        ctr=0
        while abs ( mu_hat - mu_hat_old )/ mu_hat > 1e-2 and ctr < 20:
            ctr+=1
            mu_hat_old = mu_hat
            #logger.info ( "theta hat[%d]=%s" % (ctr,list( theta_hat[:11] ) ) )
            #logger.info ( "   mu hat[%d]=%s" % (ctr, mu_hat ) )
            mu_c = numpy.abs ( self.nobs - self.nb - theta_hat ) / effs
            ## find mu_hat by finding the root of 1/L dL/dmu. We know 
            ## that the zero has to be between min(mu_c) and max(mu_c).
            lower,upper = 0.*max(mu_c),3.*max(mu_c)
            lower_v = self.dLdMu ( lower, effs, theta_hat )
            upper_v = self.dLdMu ( upper, effs, theta_hat )
            total_sign = numpy.sign ( lower_v * upper_v )
            if total_sign > -.5:
                if upper_v < lower_v < 0.:
                    ## seems like we really want to go for mu_hat = 0.
                    return 0. / lumi
                logger.debug ( "weird. cant find a zero in the Brent bracket for finding mu(hat). Let me try with a very small value." )
                lower = 1e-4*max(mu_c)
                lower_v = self.dLdMu ( lower, effs, theta_hat )
                total_sign = numpy.sign ( lower_v * upper_v )
                if total_sign > -.5:
                    logger.error ( "cant find zero in Brentq bracket. l,u=%s,%s" % ( lower, upper ) )
                    self.plotMuHatRootFinding( effs, theta_hat, lumi, None ) 
                    sys.exit()
                # return 0. / lumi
            self.timer["brentq"]-=time.time()
            mu_hat = optimize.brentq ( self.dLdMu, lower, upper, args=(effs, theta_hat ) )
            self.timer["brentq"]+=time.time()
            theta_hat,err = self.findThetaHat( mu_hat * effs )
            ctr+=1
            
        self.timer["find_mu_hat"]+=time.time()
        return mu_hat / lumi

    def getSigmaMu ( self, effs, lumi, mu_hat, theta_hat ):
        """
        get a rough estimate for the variance of mu around mu_max
        FIXME need to do this more thorougly.
        """
        # print ( "lumi=",lumi )
        #print ( "mu_hat=%s" % mu_hat )
        #print ( "theta_hat=%s" % str(theta_hat) )
        #print ( "eff=%s" % effs )
        #print ( "nb=%s" % self.nb )
        #print ( "mu_hat * effs=%s" % (mu_hat * effs + self.nb + theta_hat[0] ) )
        #num = mu_hat * effs + self.nb + theta_hat[0]
        # fisher = math.sqrt ( sum ( num**2 / ( self.nobs * effs**2 ) ) )
        #fisher = math.sqrt ( sum ( num ) )
        #logger.error ( "sigma_mu=%s" % ret )
        #print ( "fisher=%s, %s" % (fisher, fisher/lumi) )
        #return fisher # / lumi
        # logger.info ( "cov = %s" % (ret) )
        if type ( effs ) in [ list, numpy.ndarray ]:
            s_effs = sum ( effs )
        sgm_mu = math.sqrt ( sum (self.nobs ) + sum ( numpy.diag ( self.covb ) ) ) / s_effs / lumi
        #logger.error ( "sgm_mu = %s, sqrt(nobs)=%s, sqrt(cov)=%s" % \
        #               ( sgm_mu, math.sqrt ( sum ( self.nobs ) ), math.sqrt ( sum ( numpy.diag ( self.covb ) ) ) ) )
        # logger.error ( "num=%s" % num )
        return sgm_mu


    def convert ( self, x ):
        """ turn everything into numpy arrays """
        if type ( x ) in ( list, tuple ):
            return array ( x )
        return x

    #Define integrand (gaussian_(bg+signal)*poisson(nobs)):
    # def prob(x0, x1 )
    def probMV( self, *thetaA ):
        """ probability, for nuicance parameters theta """
        theta = numpy.array ( thetaA )
        # ntot = self.nb + self.nsig
        xtot = theta + self.ntot
        for ctr,i in enumerate ( xtot ):
            if i==0. or i<0.:
                # logger.debug ( "encountered a zero in probMV at position %d. Will set to small value." % ctr )
                xtot[ctr]=1e-30
        poisson = numpy.exp(self.nobs*numpy.log(xtot) - theta - self.ntot - self.gammaln)
        # poisson = numpy.exp(self.nobs*numpy.log(xtot) - theta - ntot - special.gammaln(self.nobs + 1))
        gaussian = stats.multivariate_normal.pdf(theta,mean=[0.]*len(theta),cov=(self.covb+numpy.diag(self.deltas**2)))
        ret = gaussian * ( reduce(lambda x, y: x*y, poisson) )
        return ret

    def nll ( self, theta ):
        """ probability, for nuicance parameters theta, 
        as a negative log likelihood. """
        xtot = theta + self.ntot
        xtot[xtot<=0.] = 1e-30 ## turn zeroes to small values
        # logger.info  ( "theta=%s" % theta )
        logpoisson = self.nobs*numpy.log(xtot) - xtot - self.gammaln
        loggaussian = stats.multivariate_normal.logpdf(theta,mean=[0.]*len(theta),cov=self.cov_tot)
        nll_ = - loggaussian - sum(logpoisson)
        return nll_

    def nllprime ( self, theta ):
        """ the derivative of nll as a function of the thetas. 
        Makes it easier to find the maximum likelihood. """
        xtot = theta + self.ntot
        xtot[xtot<=0.] = 1e-30 ## turn zeroes to small values
        nllp_ = self.ones - self.nobs / xtot + numpy.dot( theta , self.weight )
        # nllp_ = 1. - self.nobs / xtot + numpy.dot( theta , weight )
        # logger.debug ( "nllp_=%s" % nllp_ )
        #logger.info ( "nllp: covb^-1=%s, nllp=%s, theta=%s" % ( weight, nllp_, theta ) )
        return nllp_

    def nllHess ( self, theta ):
        """ the Hessian of nll as a function of the thetas. 
        Makes it easier to find the maximum likelihood. """
        xtot = theta + self.ntot
        xtot[xtot<=0.] = 1e-30 ## turn zeroes to small values
        nllh_ = self.weight + numpy.diag ( self.nobs / (xtot**2) )
        return nllh_

    #Define integrand (gaussian_(bg+signal)*poisson(nobs)):
    def prob( self, x ):
    # def prob( self, theta ):
        """ probability for x = ntot + theta """
        """ FIXME this should be rewritten in terms of theta. """
        ntot = self.nb + self.nsig
        poisson = exp(self.nobs*log(x) - x - math.lgamma(self.nobs + 1))
        # poisson = exp(self.nobs*log(theta+ntot) - theta - ntot - math.lgamma(self.nobs + 1))
        gaussian = stats.norm.pdf( x, loc=ntot,\
                                   scale=sqrt(self.covb+self.deltas**2))
        #gaussian = stats.norm.pdf( theta, loc=0.,\
        #                           scale=sqrt(self.covb+self.deltas**2))
        return poisson*gaussian

    def getThetaHat ( self, nobs, nb, nsig, covb, deltas, max_iterations ):
            """ Compute nuisance parameter theta that
            maximizes our likelihood (poisson*gauss).  """
            self.nsig = nsig
            self.deltas = deltas
            sigma2 = covb + numpy.diag ( deltas**2 )
            ## for now deal with variances only
            ntot = nb + nsig
            cov = numpy.matrix ( sigma2 )
            weight = cov**(-1) ## weight matrix
            diag_cov = numpy.diag(cov)
            # first: no covariances:
            q = diag_cov * ( ntot - nobs )
            p = ntot + diag_cov
            thetamaxes = []
            thetamax = -p/2. * ( 1 - sign(p) * sqrt ( 1. - 4*q / p**2 ) )
            thetamaxes.append ( thetamax )
            ndims = len(p)
            def distance ( theta1, theta2 ):
                for ctr,i in enumerate ( theta1 ):
                    if i == 0.:
                        theta1[ctr]=1e-20
                for ctr,i in enumerate ( theta2 ):
                    if i == 0.:
                        theta2[ctr]=1e-20
                return sum ( numpy.abs(theta1 - theta2) / numpy.abs ( theta1+theta2 ) )

            for ctr in range(max_iterations):
                q = diag_cov * ( ntot - nobs )
                p = ntot + diag_cov
                for i in range(ndims):
                    #q[i] = diag_cov[i] * ( ntot[i] - nobs[i] )
                    #p[i] = ntot[i] + diag_cov[i]
                    for j in range(ndims):
                        if i==j: continue
                        dq = thetamax[j]*ntot[i]*diag_cov[i]*weight[i,j]
                        dp = thetamax[j]*weight[i,j]*diag_cov[i] 
                        if abs ( dq / q[i] ) > .3:
                            #logger.warning ( "too big a change in iteration." )
                            dq=numpy.abs( .3 * q[i] ) * numpy.sign ( dq )
                        if abs ( dp / p[i] ) > .3:
                            #logger.warning ( "too big a change in iteration." )
                            dp=numpy.abs( .3 * p[i] ) * numpy.sign ( dp )
                        q[i] += dq
                        p[i] += dp
                    thetamax = -p/2. * ( 1 - sign(p) * sqrt ( 1. - 4*q / p**2 ) )
                thetamaxes.append ( thetamax )
                if len(thetamaxes)>2:
                    d1 = distance ( thetamaxes[-2], thetamax )
                    d2 = distance ( thetamaxes[-3], thetamaxes[-2] )
                    if d1 > d2:
                        logger.error ( "diverging when computing thetamax: %f > %f" % ( d1, d2 ) )
                        sys.exit()
                    if d1 < 1e-5:
                        return thetamax
            return thetamax

    def findThetaHat ( self, nsig, deltas=None ):
            """ Compute nuisance parameter theta that maximizes our likelihood
                (poisson*gauss).
            """
            if type(deltas) == type(None):
                 deltas = self.deltas_default * nsig
                 if type(nsig) in [ list, numpy.ndarray, numpy.array ]:
                    # deltas = numpy.array ( [0.] * len(nsig) )
                    deltas = numpy.array ( self.deltas_default * nsig )
            ## first step is to disregard the covariances and solve the
            ## quadratic equations
            self.timer["theta_hat_ini"]-=time.time()
            ini = self.getThetaHat ( self.nobs, self.nb, nsig, self.covb, deltas, 0 ) 
            self.timer["theta_hat_ini"]+=time.time()
            self.cov_tot = self.covb+numpy.diag(self.deltas**2)
            self.weight = numpy.linalg.inv ( self.cov_tot )
            self.ntot = self.nb + self.nsig
            self.ones = 1.
            if type ( self.nobs) in [ list, numpy.ndarray ]:
                self.ones = numpy.ones ( len (self.nobs) )
            self.gammaln = special.gammaln(self.nobs + 1)
            try:
                self.timer["fmin2"]-=time.time()
                # first use NCG
                ret_c = optimize.fmin_ncg ( self.nll, ini, fprime=self.nllprime, fhess=self.nllHess, full_output=True, disp=0 )
                # then always continue with TNC
                if type ( self.nobs ) in [ int, float ]:
                    bounds = [ ( -10*self.nobs, 10*self.nobs ) ]
                else:
                    bounds = [ ( -10*x, 10*x ) for x in self.nobs ]
                ini = ret_c
                ret_c = optimize.fmin_tnc ( self.nll, ret_c[0], fprime=self.nllprime, disp=0, bounds=bounds )
                if ret_c[-1] not in [ 0, 1, 2 ]:
                    is_expected=False
                    if ( self.nobs == self.nb ).all():
                        is_expected=True
                    # logger.debug ( "tnc failed also: %s [%d]" % ( str(ret_c[-2:]), is_expected ) )
                    #logger.info ( "ini was %s" % str(ini[:]) )
                    return ret_c[0],ret_c[-1]
                    # ret_c = optimize.fmin_l_bfgs_b ( self.nll, ret_c[0], fprime=self.nllprime, disp=5, bounds=bounds )
                    # logger.info ( "ret3 %s" % str(ret3[:]) )
                    # sys.exit()
                else:
                    return ret_c[0],0
                    logger.debug ( "tnc worked." )
                
                ret = ret_c[0]
                self.timer["fmin2"]+=time.time()
                return ret,-2
            except Exception as e:
                logger.error ( "exception: %s. ini[-3:]=%s" % (e,ini[-3:]) )
                logger.error ( "cov-1=%s" % (self.covb+numpy.diag(deltas))**(-1) )
                sys.exit()
            return ini,-1

    def plotLTheta ( self, nsig=None ):
        """ plot the likelihood, but as a function of the nuisance theta! 
        :param nsig: plot it at nsig. If None, plot at mu_hat=0.
        """
        if not self.debug_mode: return
        if ( self.nobs == self.nb ).all():
            # only plot we this isnt the 'expected' case.
            return 
        import ROOT
        # oldsig = self.nsig
        s_mu="#hat{#mu}"
        if type(nsig) == type(None):
            nsig = [0.] * len(self.nobs)
            s_mu="#mu=0"
        self.nsig = nsig
        theta_hat,err = self.findThetaHat ( nsig, self.deltas )
        theta_hat_1 = .9 * theta_hat
        v_tm = self.probMV ( *theta_hat )
        v_tm1 = self.probMV ( *theta_hat_1 )
        #logger.error ( "point  at %s %s" % ( theta_hat[:10], v_tm ) )
        #logger.error ( "point0 at %s %s" % ( theta_hat_0[:10], v_tm0 ) )
        #logger.error ( "point1 at %s %s" % ( theta_hat_1[:10], v_tm1 ) )
        rng = list ( numpy.arange ( 0.3, 1.7, 0.05 ) )
        t=ROOT.TGraph(len(rng))
        t.SetTitle ( "-ln L(%s,#theta/#hat{#theta}), %d SRs" % (s_mu, len(theta_hat) ) )
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
                t.SetPoint( ctr, i, nll )
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
        # self.nsig = oldsig

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
            theta_hat,err = self.findThetaHat ( nsig, deltas )
            #print ( "found xmax at", xmax, "l=",self.probMV(*xmax) )

            #Define initial integration range:
            nrange = 5.
            low = theta_hat-nrange*sqrt(dsigma2)
            low[low<0.] = 0. ## set all negatives to zero!
            up = theta_hat+nrange*sqrt(dsigma2)
            ## first integral can be coarse, we will anyhow do another round
            opts = { "epsabs":1e-2, "epsrel":1e-1 }
            logger.error ( "compute llhd for %s" % nsig )
            like = integrate.nquad( self.probMV, list(zip(low,up)),opts=opts )[0]
            logger.error ( "done compute llhd for %s" % nsig )
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
            Warning: not normalized. 
            Returns profile likelihood and error code (0=no error)
        """
        self.timer["profile"]-=time.time()
        nsig = self.convert ( nsig )
        if type( nsig ) in [ int, float, numpy.float64 ]:
            if type(deltas) == type(None):
                deltas = self.deltas_default * nsig 
            ## FIXME write profiled likelihood for 1d case.
            return self._likelihood1d( nsig, deltas ),0
        if type(deltas) == type(None):
            deltas = numpy.array ( self.deltas_default * nsig )
        # compute the profiled (not normalized) likelihood of observing
        # nsig signal events
        theta_hat,err = self.findThetaHat ( nsig, deltas )
        ret = self.probMV ( *theta_hat )
        # logger.error ( "theta_hat, err=%s, %s, %s" % ( theta_hat[0], err, ret ) )
        self.timer["profile"]+=time.time()
        return (ret,err)

    def likelihood ( self, nsig, deltas = None ):
        """ compute likelihood for nsig, profiling the nuisances """
        nsig = self.convert ( nsig )
        if type(deltas) == type(None):
            deltas = self.deltas_default * nsig ## FIXME for backwards compatibility
        if type( nsig ) in [ int, float, numpy.float64 ]:
            return self._likelihood1d( nsig, deltas )
        return self.profileLikelihood ( nsig, deltas )[0]
        # return self._mvLikelihood( nsig, deltas )

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
                deltas = self.deltas_default * nsig
            # Compute the likelhood for the null hypothesis (signal hypothesis) H0:
            llhd = self.likelihood( nsig, deltas )

            #print ( "deltas=",deltas )
            #print ( "nsig=",nsig )
            #Percentual signal error:
            denom = nsig
            if type(denom)==float and denom == 0.:
                denom = 1e10 ## when nsig is zero, then so is deltas
            if type(denom)==numpy.ndarray:
                ## when nsig is zero, then so is deltas
                denom [ denom == 0. ] = 1e10

            deltas_pct = deltas / denom ## float(nsig)

            # Compute the maximum likelihood H1, which sits at nsig = nobs - nb
            # (keeping the same % error on signal):
            dn = self.nobs-self.nb
            deltas = deltas_pct*( dn )
            maxllhd = self.likelihood( dn, deltas )

            # Return infinite likelihood if it is zero
            # This can happen in case e.g. nb >> nobs
            if llhd == 0.:
                ## set to insanely small number
                llhd = 1e-300

            chi2=-2*log(llhd/maxllhd)

            # Return the test statistic -2log(H0/H1)
            return chi2

if __name__ == "__main__":
    computer = LikelihoodComputer ( nobs_, nb_, deltab_**2 )
    print ( "1d, computer:", computer.likelihood( nsig_, deltas_ )  )
    computer = LikelihoodComputer ( [nobs_], [nb_], numpy.diag([deltab_**2]) )
    print ( "mv 1d, computer:",computer.likelihood( [nsig_], deltas_) )
    print ( "mv 1d, chi2:",computer.chi2 ( [nsig_] ) )
    cov = numpy.diag ([deltab_**2,deltab_**2])
    cov[0,1] = .01
    cov[1,0] = .01
    computer = LikelihoodComputer ( [nobs_,nobs_], [nb_,nb_], cov )
    l = computer.likelihood ( array([nsig_,nsig_]), deltas_  )
    print ( "mv 2d, computer:", l )
    print ( "mv 2d, chi2:", computer.chi2 ( [nsig_,nsig_] ) )
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
