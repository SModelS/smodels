#!/usr/bin/env python

"""
.. module:: SimplifiedLikelihoods
   :synopsis: Code that implements the simplified likelihoods as presented
              in CMS-NOTE-2017-001, see https://cds.cern.ch/record/2242860,
              and FIXME insert arXiv reference here.

.. moduleauthor:: Andy Buckley <andy.buckley@cern.ch>
.. moduleauthor:: Sylvain Fichet <sylvain.fichet@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Nickolas Wardle <nicholas.wardle@cern.ch>

"""

from __future__ import print_function
from scipy import stats, optimize, integrate, special
from numpy  import sqrt, exp, log, sign, array, matrix, ndarray, floor
from functools import reduce
import numpy as NP
import math
import sys, copy

def getLogger():
    """ configure the logging facility. Maybe adapted to fit into
        your framework. """
    import logging
    logger = logging.getLogger("SL")
    formatter = logging.Formatter('%(module)s - %(levelname)s: %(message)s')
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    ch.setLevel(logging.DEBUG)
    logger.addHandler(ch)
    return logger

def importUnits():
    """ Define units (fb,pb). Use unum, if possible. """
    from smodels.tools.physicsUnits import fb, pb
    return fb,pb
    try:
        import unum
        try:
            fb = unum.Unum.unit('fb')
            pb = unum.Unum.unit('pb', 1000 * fb)
            return fb,pb
        except unum.NameConflictError as e:
            # already defined? reset!
            unum.Unum("fb").reset()
            unum.Unum("pb").reset()
            fb = unum.Unum.unit('fb')
            pb = unum.Unum.unit('pb', 1000 * fb)
            return fb,pb
    except ImportError as e:
        logger.error ( "unum not installed. For now I will work without units. You"
                       " are however advised to install the unum module." )
        fb,pb=1.,1000.
        return fb,pb

logger=getLogger()
fb,pb=importUnits()

class Model:
    """ A very simple data container to collect all the data
        needed to fully define a specific statistical model """

    def __init__ ( self, data, backgrounds, covariance, third_moment=None,
                         efficiencies=None, name="model", deltas_rel=None ):
        """
        :param data: number of observed events per dataset
        :param backgrounds: expected bg per dataset
        :param covariance: uncertainty in background, as a covariance matrix
        :param efficiencies: dataset effective signal efficiencies
        :param name: give the model a name, just for convenience
        :param deltas_rel: the relative assumed error on signal hypotheses.
        If None, then no error is assumed.
        """
        self.data = self.convert ( data )
        self.backgrounds = self.convert ( backgrounds )
        self.n = len ( self.data )
        self.covariance = self.convertCov ( covariance )
        self.efficiencies = self.convert ( efficiencies )
        self.third_moment = self.convert ( third_moment )
        if type(self.third_moment) != type(None) and NP.sum ( [ abs(x) for x in self.third_moment ] ) < 1e-10:
            ## all zeroes? then we have no third_moment
            self.third_moment = None
        self.name = name
        self.deltas_rel = deltas_rel
        if deltas_rel == None:
            self.deltas_rel = 1e-10
        self.computeABC()

    def isScalar ( self, obj ):
        """ determine if obj is a scalar (float or int) """
        if type(obj) == ndarray:
            ## need to treat separately since casting array([0.]) to float works
            return False
        try:
            _ = float(obj)
            return True
        except:
            pass
        return False

    def convert ( self, obj ):
        """ convert everything to NP arrays """
        if type(obj) == type(None):
            return obj
        if self.isScalar(obj):
            return array ( [ obj ] )
        return array ( obj )

    def __str__ ( self ):
        return self.name + " (%d dims)" % self.n

    def convertCov ( self, obj ):
        if self.isScalar(obj):
            return array ( [ [ obj ] ] )
        if type(obj[0]) == float:
            ## if the matrix is flattened, unflatten it.
            return array ( [ obj[self.n*i:self.n*(i+1)] for i in range(self.n) ] )
        return obj

    def computeABC ( self ):
        """ compute the terms A, B, C, rho, V. Corresponds with
            Eqs. 1.27-1.30 in the second paper. """
        self.V = self.covariance
        if type(self.third_moment) == type(None):
            self.A = None
            self.B = None
            self.C = None
            return

        covD = self.diagCov()
        C=[]
        for m2,m3 in zip(covD, self.third_moment ):
            if m3 == 0.: m3 = 1e-30
            k=-NP.sign(m3)*sqrt(2.*m2 )
            dm = sqrt ( 8.*m2**3/m3**2 - 1. )
            C.append ( k*NP.cos ( 4.*NP.pi/3. + NP.arctan(dm) / 3. ) )
        self.C=NP.array(C)
        self.B = sqrt ( covD - 2*self.C**2 )
        self.A = self.backgrounds - self.C
        self.rho = NP.array ( [ [0.]*self.n ]*self.n )
        for x in range(self.n):
            for y in range(x,self.n):
                bxby=self.B[x]*self.B[y]
                cxcy=self.C[x]*self.C[y]
                e=(4.*cxcy)**(-1)*(sqrt( bxby**2+8*cxcy*self.covariance[x][y])-bxby)
                self.rho[x][y]=e
                self.rho[y][x]=e

        self.sandwich()
        # self.V = sandwich ( self.B, self.rho )

    def sandwich ( self ):
        """ sandwich product """
        ret = NP.array ( [ [0.]*len(self.B) ]*len(self.B) )
        for x in range(len(self.B)):
            for y in range(x,len(self.B)):
                T=self.B[x]*self.B[y]*self.rho[x][y]
                ret[x][y]=T
                ret[y][x]=T
        self.V = ret

    def isLinear(self):
        """ model is linear, i.e. no quadratic term in poissonians """
        return type(self.C) == type(None)

    def diagCov ( self ):
        """ diagonal elements of covariance matrix. Convenience function. """
        return NP.diag ( self.covariance )

    def correlations ( self ):
        """ correlation matrix, computed from covariance matrix.
            convenience function. """
        if hasattr ( self, "corr" ):
            return self.corr
        self.corr = copy.deepcopy ( self.covariance )
        for x in range(self.n):
            self.corr[x][x]=1.
            for y in range(x+1,self.n):
                rho=self.corr[x][y]/sqrt(self.covariance[x][x]*self.covariance[y][y])
                self.corr[x][y]=rho
                self.corr[y][x]=rho
        return self.corr

    def signals ( self, mu ):
        """ returns the signal cross sections, for all datasets,
        given signal strength mu. """
        return mu * self.efficiencies

class LikelihoodComputer:

    """ the default value for delta_s, the assumed relative error
    on signal yields. """
    debug_mode = False

    def __init__ ( self, model, ntoys = 1000 ):
        """
        :param model: a Model object.
        :param ntoys: number of toys when marginalizing
        """
        self.model = model
        self.ntoys = ntoys

    def dLdMu ( self, mu, effs, theta_hat ):
        """ d (ln L)/d mu, if L is the likelihood. The function
            whose root gives us muhat, i.e. the mu that maximizes
            the likelihood. """
        denominator = mu*effs + self.model.backgrounds + theta_hat
        for ctr,d in enumerate ( denominator ):
            if d==0.:
                if (self.model.data[ctr]*effs[ctr]) == 0.:
                    logger.debug ( "zero denominator, but numerator also zero, so we set denom to 1." )
                    denominator[ctr]=1.
                else:
                    logger.error ( "we have a zero value in the denominator at pos %d, with a non-zero numerator. dont know how to handle." % ctr )
        ret = self.model.data * effs / denominator - effs
        if type (ret ) in [ array, ndarray, list ]:
            ret = sum ( ret )
        return ret

    def findMuHat ( self, effs, lumi=1. ):
        """
        find the most likely signal strength mu
        :param lumi: return yield (lumi=1.) or cross section ("real" lumi)
        :returns: mu_hat, either as signal yield (lumi=1.), or as cross section.
        """
        if ( self.model.backgrounds == self.model.data ).all():
            return 0. / lumi
        if type ( effs ) in [ list, ndarray ]:
            effs = array ( effs )
        effs[effs==0.]=1e-20
        if sum ( effs<0. ):
            logger.error ( "Negative efficiencies!" )
            sys.exit()
        ## we need a very rough initial guess for mu(hat), to come
        ## up with a first theta
        self.nsig = array ( [0.] * len(self.model.data ) )
        self.deltas = array ( [0.] * len(self.model.data) )
        ## we start with theta_hat being all zeroes
        theta_hat = array ( [0.] * len(self.model.data ) )
        mu_hat_old, mu_hat = 0., 1.
        ctr=0
        while abs ( mu_hat - mu_hat_old )/ mu_hat > 1e-2 and ctr < 20:
            ctr+=1
            mu_hat_old = mu_hat
            #logger.info ( "theta hat[%d]=%s" % (ctr,list( theta_hat[:11] ) ) )
            #logger.info ( "   mu hat[%d]=%s" % (ctr, mu_hat ) )
            mu_c = NP.abs ( self.model.data - self.model.backgrounds - theta_hat ) / effs
            ## find mu_hat by finding the root of 1/L dL/dmu. We know
            ## that the zero has to be between min(mu_c) and max(mu_c).
            lower,upper = 0.*max(mu_c),3.*max(mu_c)
            lower_v = self.dLdMu ( lower, effs, theta_hat )
            upper_v = self.dLdMu ( upper, effs, theta_hat )
            total_sign = NP.sign ( lower_v * upper_v )
            if total_sign > -.5:
                if upper_v < lower_v < 0.:
                    ## seems like we really want to go for mu_hat = 0.
                    return 0. / lumi
                logger.debug ( "weird. cant find a zero in the Brent bracket "\
                               "for finding mu(hat). Let me try with a very small"
                               " value." )
                lower = 1e-4*max(mu_c)
                lower_v = self.dLdMu ( lower, effs, theta_hat )
                total_sign = NP.sign ( lower_v * upper_v )
                if total_sign > -.5:
                    logger.error ( "cant find zero in Brentq bracket. l,u=%s,%s" % \
                                   ( lower, upper ) )
                    sys.exit()
            mu_hat = optimize.brentq ( self.dLdMu, lower, upper, args=(effs, theta_hat ) )
            theta_hat,err = self.findThetaHat( mu_hat * effs )
            ctr+=1

        return mu_hat / lumi

    def getSigmaMu ( self, effs, lumi, mu_hat, theta_hat ):
        """
        get a rough estimate for the variance of mu around mu_max
        FIXME need to do this more thorougly.
        """
        if type ( effs ) in [ list, ndarray ]:
            s_effs = sum ( effs )
        sgm_mu = sqrt ( sum (self.model.data ) + sum ( NP.diag ( self.model.covariance ) ) ) / s_effs / lumi
        return sgm_mu

    #Define integrand (gaussian_(bg+signal)*poisson(nobs)):
    # def prob(x0, x1 )
    def probMV( self, nll, *thetaA ):
        """ probability, for nuicance parameters theta
        :params nll: compute negative log likelihood """
        theta = array ( thetaA )
        # ntot = self.model.backgrounds + self.nsig
        # lmbda = theta + self.ntot ## the lambda for the Poissonian
        if self.model.isLinear():
            lmbda = self.model.backgrounds + self.nsig + theta
        else:
            lmbda = self.nsig + self.model.A + theta + self.model.C * theta**2 / self.model.B**2
        lmbda[lmbda<=0.] = 1e-30 ## turn zeroes to small values
        if nll:
            #poisson = self.model.data * log ( lmbda ) - lmbda #- self.gammaln
            poisson = stats.poisson.logpmf ( self.model.data, lmbda )
            #print ( "p",poisson,poisson2 )
        else:
            poisson = stats.poisson.pmf ( self.model.data, lmbda )
            #print ( "nonll",poisson )
        try:
            if nll:
                # gaussian = -.5 * (theta*NP.linalg.inv ( self.model.V ) * theta)[0][0]
                gaussian = stats.multivariate_normal.logpdf(theta,mean=[0.]*len(theta),cov=self.model.V)
                # print ( "gaussian=",gaussian,gaussian2)
                ret = - gaussian - sum(poisson)
            else:
                gaussian = stats.multivariate_normal.pdf(theta,mean=[0.]*len(theta),cov=self.model.V)
                ret = gaussian * ( reduce(lambda x, y: x*y, poisson) )
            return ret
        except ValueError as e:
            logger.error ( "ValueError %s, %s" % ( e, self.model.V ) )
            sys.exit()

    def nll ( self, theta ):
        """ probability, for nuicance parameters theta,
        as a negative log likelihood. """
        return self.probMV(True,*theta)

    def nllprime ( self, theta ):
        """ the derivative of nll as a function of the thetas.
        Makes it easier to find the maximum likelihood. """
        if self.model.isLinear():
            xtot = theta + self.model.backgrounds + self.nsig
            xtot[xtot<=0.] = 1e-30 ## turn zeroes to small values
            nllp_ = self.ones - self.model.data / xtot + NP.dot( theta , self.weight )
            return nllp_
        lmbda = self.nsig + self.model.A + theta + self.model.C * theta**2 / self.model.B**2
        lmbda[lmbda<=0.] = 1e-30 ## turn zeroes to small values
        # nllp_ = ( self.ones - self.model.data / lmbda + NP.dot( theta , self.weight ) ) * ( self.ones + 2*self.model.C * theta / self.model.B**2 )
        T=self.ones + 2*self.model.C/self.model.B**2*theta
        nllp_ = T - self.model.data / lmbda * ( T ) + NP.dot( theta , self.weight )
        return nllp_

    def nllHess ( self, theta ):
        """ the Hessian of nll as a function of the thetas.
        Makes it easier to find the maximum likelihood. """
        # xtot = theta + self.ntot
        if self.model.isLinear():
            xtot = theta + self.model.backgrounds + self.nsig
            xtot[xtot<=0.] = 1e-30 ## turn zeroes to small values
            nllh_ = self.weight + NP.diag ( self.model.data / (xtot**2) )
            return nllh_
        lmbda = self.nsig + self.model.A + theta + self.model.C * theta**2 / self.model.B**2
        lmbda[lmbda<=0.] = 1e-30 ## turn zeroes to small values
        T=self.ones + 2*self.model.C/self.model.B**2*theta
        # T_i = 1_i + 2*c_i/B_i**2 * theta_i
        nllh_ = self.weight + NP.diag ( self.model.data * T**2 / (lmbda**2) ) - NP.diag ( self.model.data / lmbda * 2 * self.model.C / self.model.B**2 ) + NP.diag ( 2*self.model.C/self.model.B**2 )
        return nllh_

    def getThetaHat ( self, nobs, nb, nsig, covb, deltas, max_iterations ):
            """ Compute nuisance parameter theta that
            maximizes our likelihood (poisson*gauss).  """
            self.nsig = nsig
            self.deltas = deltas
            sigma2 = covb + NP.diag ( deltas**2 )
            ## for now deal with variances only
            ntot = nb + nsig
            cov = NP.matrix ( sigma2 )
            weight = cov**(-1) ## weight matrix
            diag_cov = NP.diag(cov)
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
                return sum ( NP.abs(theta1 - theta2) / NP.abs ( theta1+theta2 ) )

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
                            dq=NP.abs( .3 * q[i] ) * NP.sign ( dq )
                        if abs ( dp / p[i] ) > .3:
                            #logger.warning ( "too big a change in iteration." )
                            dp=NP.abs( .3 * p[i] ) * NP.sign ( dp )
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
                 deltas = self.model.deltas_rel * nsig
                 if type(nsig) in [ list, ndarray, array ]:
                    # deltas = array ( [0.] * len(nsig) )
                    deltas = array ( self.model.deltas_rel * nsig )
            ## first step is to disregard the covariances and solve the
            ## quadratic equations
            ini = self.getThetaHat ( self.model.data, self.model.backgrounds, nsig, self.model.covariance, deltas, 0 )
            self.cov_tot = self.model.covariance+NP.diag(self.deltas**2)
            self.ntot = self.model.backgrounds + self.nsig
            if not self.model.isLinear():
                self.cov_tot = self.model.V+NP.diag(self.deltas**2)
                self.ntot = None
            self.weight = NP.linalg.inv ( self.cov_tot )
            self.ones = 1.
            if type ( self.model.data) in [ list, ndarray ]:
                self.ones = NP.ones ( len (self.model.data) )
            self.gammaln = special.gammaln(self.model.data + 1)
            try:
                ret_c = optimize.fmin_ncg ( self.nll, ini, fprime=self.nllprime,
                                       fhess=self.nllHess, full_output=True, disp=0 )
                # then always continue with TNC
                if type ( self.model.data ) in [ int, float ]:
                    bounds = [ ( -10*self.model.data, 10*self.model.data ) ]
                else:
                    bounds = [ ( -10*x, 10*x ) for x in self.model.data ]
                ini = ret_c
                ret_c = optimize.fmin_tnc ( self.nll, ret_c[0], fprime=self.nllprime,
                                            disp=0, bounds=bounds )
                # print ( "[findThetaHat] mu=%s bg=%s data=%s V=%s, nsig=%s theta=%s, nll=%s" % ( self.nsig[0]/self.model.efficiencies[0], self.model.backgrounds, self.model.data,self.model.covariance, self.nsig, ret_c[0], self.nll(ret_c[0]) ) )
                if ret_c[-1] not in [ 0, 1, 2 ]:
                    is_expected=False
                    if ( self.model.data == self.model.backgrounds ).all():
                        is_expected=True
                    return ret_c[0],ret_c[-1]
                else:
                    return ret_c[0],0
                    logger.debug ( "tnc worked." )

                ret = ret_c[0]
                return ret,-2
            except Exception as e:
                logger.error ( "exception: %s. ini[-3:]=%s" % (e,ini[-3:]) )
                logger.error ( "cov-1=%s" % (self.model.covariance+NP.diag(deltas))**(-1) )
                sys.exit()
            return ini,-1

    def marginalizedLLHD1D(self, nsig, deltas, nll):
            """
            Return the likelihood (of 1 signal region) to observe nobs events given the
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

            self.deltas = deltas
            self.sigma2 = self.model.covariance + deltas**2
            self.sigma_tot = sqrt( self.sigma2 )
            self.lngamma = math.lgamma(self.model.data[0] + 1)
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
            def prob( x, nsig ):
                poisson = exp(self.model.data*log(x) - x - self.lngamma )
                gaussian = stats.norm.pdf(x,loc=self.model.backgrounds+nsig,scale=self.sigma_tot)

                return poisson*gaussian

            #Compute maximum value for the integrand:
            xm = self.model.backgrounds + nsig - self.sigma2
            #If nb + nsig = sigma2, shift the values slightly:
            if xm == 0.:
                xm = 0.001
            xmax = xm*(1.+sign(xm)*sqrt(1. + 4.*self.model.data*self.sigma2/xm**2))/2.

            #Define initial integration range:
            nrange = 5.
            a = max(0.,xmax-nrange*sqrt(self.sigma2))
            b = xmax+nrange*self.sigma_tot
            like = integrate.quad(prob,a,b,(nsig), epsabs=0.,epsrel=1e-3)[0]

            #Increase integration range until integral converges
            err = 1.
            while err > 0.01:
                like_old = like
                nrange = nrange*2
                a = max(0.,xmax-nrange*self.sigma_tot)
                b = xmax+nrange*self.sigma_tot
                like = integrate.quad(prob,a,b,(nsig),
                                          epsabs=0.,epsrel=1e-3)[0]
                err = abs(like_old-like)/like

            #Renormalize the likelihood to account for the cut at x = 0.
            #The integral of the gaussian from 0 to infinity gives:
            #(1/2)*(1 + Erf(mu/sqrt(2*sigma2))), so we need to divide by it
            #(for mu - sigma >> 0, the normalization gives 1.)
            norm = (1./2.)*(1. + special.erf((self.model.backgrounds+nsig)/sqrt(2.*self.sigma2)))
            like = like/norm

            if nll:
                like = - log ( like )

            return like[0][0]

    def marginalizedLikelihood( self, nsig, nll, deltas ):
            """ compute the marginalized likelihood of observing nsig signal event"""
            if self.model.isLinear() and self.model.n == 1: ## 1-dimensional non-skewed llhds we can integrate analytically
                return self.marginalizedLLHD1D ( nsig, deltas, nll )

            if type(deltas) == type(None):
                deltas = array ( self.model.deltas_rel * nsig )
            vals=[]
            self.gammaln = special.gammaln(self.model.data + 1)
            thetas = stats.multivariate_normal.rvs(mean=[0.]*self.model.n,
                          cov=(self.model.V+NP.diag(deltas**2)),
                          size=self.ntoys ) ## get ntoys values
            for theta in thetas :
                if self.model.isLinear():
                    lmbda = nsig + self.model.backgrounds + theta
                else:
                    lmbda = nsig + self.model.A + theta + self.model.C*theta**2/self.model.B**2
                if self.model.isScalar ( lmbda ): lmbda = array ( [ lmbda ] )
                for ctr,v in enumerate ( lmbda ):
                    if v<=0.: lmbda[ctr]=1e-30
#                    print ( "lmbda=",lmbda )
                poisson = self.model.data*NP.log(lmbda) - lmbda - self.gammaln
                # poisson = NP.exp(self.model.data*NP.log(lmbda) - lmbda - self.model.backgrounds - self.gammaln)
                vals.append ( NP.exp ( sum(poisson) ) )
                #vals.append ( reduce(lambda x, y: x*y, poisson) )
            mean = NP.mean ( vals )
            if nll:
                mean = - log ( mean )
            return mean


    def profileLikelihood ( self, nsig, nll, deltas ):
        """ compute the profiled likelihood for nsig.
            Warning: not normalized.
            Returns profile likelihood and error code (0=no error)
        """
        if type(deltas) == type(None):
            deltas = array ( self.model.deltas_rel * nsig )
        # compute the profiled (not normalized) likelihood of observing
        # nsig signal events
        theta_hat,err = self.findThetaHat ( nsig, deltas )
        ret = self.probMV ( nll, *theta_hat )
        # logger.error ( "theta_hat, err=%s, %s, %s" % ( theta_hat[0], err, ret ) )
        return ret

    def likelihood ( self, nsig, deltas = None, marginalize=False, nll=False ):
        """ compute likelihood for nsig, profiling the nuisances
        :param deltas: error on signal
        :param marginalize: if true, marginalize, if false, profile
        :param nll: return nll instead of likelihood
        """
        nsig = self.model.convert ( nsig )
        self.ntot = self.model.backgrounds + nsig
        deltas = self.model.convert ( deltas )
        if marginalize:
            # p,err = self.profileLikelihood ( nsig, deltas )
            return self.marginalizedLikelihood( nsig, nll, deltas )
            # print ( "p,l=",p,l,p/l )
        else:
            return self.profileLikelihood ( nsig, nll, deltas )

    def chi2( self, nsig, deltas=None ):
            """
            Computes the chi2 for a given number of observed events nobs given
            the predicted background nb, error on this background deltab,
            expected number of signal events nsig and, if given, the error on
            signal (deltas).
            :return: chi2 (float)

            """
            nsig = self.model.convert ( nsig )
            marg=True
            if deltas == None:
                deltas = self.model.deltas_rel * nsig
            # Compute the likelhood for the null hypothesis (signal hypothesis) H0:
            llhd = self.likelihood( nsig, deltas, marginalize=marg, nll=True )

            #print ( "deltas=",deltas )
            #print ( "nsig=",nsig )
            #Percentual signal error:
            denom = nsig
            if type(denom)==float and denom == 0.:
                denom = 1e10 ## when nsig is zero, then so is deltas
            if type(denom)==ndarray:
                ## when nsig is zero, then so is deltas
                denom [ denom == 0. ] = 1e10

            deltas_pct = deltas / denom ## float(nsig)

            # Compute the maximum likelihood H1, which sits at nsig = nobs - nb
            # (keeping the same % error on signal):
            dn = self.model.data-self.model.backgrounds
            deltas = deltas_pct*( dn )
            maxllhd = self.likelihood( dn, deltas, marginalize=marg, nll=True )

            # Return infinite likelihood if it is zero
            # This can happen in case e.g. nb >> nobs
            #if llhd == 0.:
                ## set to insanely small number
            #    llhd = 1e-300
            # print ( "llhd,maxllhd=",llhd,maxllhd)

            chi2=2*(llhd-maxllhd)

            # Return the test statistic -2log(H0/H1)
            return chi2

class UpperLimitComputer:
    debug_mode = False

    def __init__ ( self, lumi, ntoys=10000, cl=.95):
        """
        :param lumi: integrated luminosity
        :param ntoys: number of toys when marginalizing
        :param cl: desired quantile for limits
        """
        self.lumi = lumi
        self.ntoys = ntoys
        self.cl = cl

    def ulSigma ( self, model, marginalize=True, toys=None, expected=False ):
        """ upper limit obtained from combined efficiencies, by using
            the q_mu test statistic from the CCGV paper (arXiv:1007.1727).

        :params marginalize: if true, marginalize nuisances, else profile them
        :params toys: specify number of toys. Use default is none
        :params expected: compute the expected value, not the observed.
        :returns: upper limit on *production* xsec (efficiencies unfolded)
        """
        if toys==None:
            toys=self.ntoys
        oldmodel = model
        if expected:
            model = copy.deepcopy ( oldmodel )
            model.data = model.backgrounds
        computer = LikelihoodComputer ( model, toys )
        mu_hat = computer.findMuHat ( model.efficiencies )
        theta_hat = computer.findThetaHat ( mu_hat * model.efficiencies )
        theta_hat0,e = computer.findThetaHat ( 0 * model.efficiencies )
        sigma_mu = computer.getSigmaMu ( model.efficiencies, 1., mu_hat, theta_hat )

        aModel = copy.deepcopy ( model )
        aModel.data = array ( [ round(x+y) for x,y in zip(model.backgrounds,theta_hat0 ) ] )
        #print ( "aModeldata=", aModel.data )
        #aModel.data = array ( [ round(x) for x in model.backgrounds ] )
        aModel.name = aModel.name + "A"
        compA = LikelihoodComputer ( aModel, toys )
        ## compute
        mu_hatA = compA.findMuHat ( aModel.efficiencies )
        # print ( "mu_hat=", mu_hat )
        if mu_hat < 0.: mu_hat = 0.
        nll0 = computer.likelihood ( model.efficiencies * mu_hat,
                                     marginalize=marginalize, nll=True )
        nll0A = compA.likelihood ( aModel.efficiencies * mu_hatA,
                                   marginalize=marginalize, nll=True )
        def root_func ( mu ):
            ## the function to minimize.
            nsig = mu*model.efficiencies
            computer.ntot = model.backgrounds + nsig
            nll = computer.likelihood ( nsig, marginalize=marginalize, nll=True )
            nllA = compA.likelihood ( nsig, marginalize=marginalize, nll=True )
            qmu =  2*( nll - nll0 )
            if qmu<0.: qmu=0.
            sqmu = sqrt (qmu)
            qA =  2*( nllA - nll0A )
            # print ( "mu: %s, qMu: %s, qA: %s nll0A: %s nllA: %s" % ( mu, qmu, qA, nll0A, nllA ) )
            if qA<0.: qA=0.
            sqA = sqrt ( qA )
            CLsb = 1. - stats.multivariate_normal.cdf ( sqmu )
            CLb = 0.
            if qA >= qmu:
                CLb =  stats.multivariate_normal.cdf ( sqA - sqmu )
            else:
                if qA == 0.:
                    CLSb = 1.
                    CLb  = 1.
                else:
                    CLsb = 1. - stats.multivariate_normal.cdf ( (qmu + qA)/(2*sqA) )
                    CLb = 1. - stats.multivariate_normal.cdf ( (qmu - qA)/(2*sqA) )
            CLs = 0.
            if CLb>0.:
                CLs = CLsb / CLb
            root = CLs - 1. + self.cl
            # print ( "mu=%s nll=%s nllo=%s qA=%s clb=%s clsb=%s cls=%s" % ( mu, nll, nll0, qA, CLb, CLsb, CLs ) )
            return root
        # print ( "model=%s " % model )

        a,b=1.5*mu_hat,2.5*mu_hat+2*sigma_mu
        ctr=0
        while True:
            while ( NP.sign ( root_func(a)* root_func(b) ) > -.5 ):
                b=1.2*b  ## widen bracket
                a=a-(b-a)*.2 ## widen bracket
                if a < 0.: a=0.
                ctr+=1
                if ctr>20: ## but stop after 20 trials
                    if toys > 2000:
                       logger.error("cannot find brent bracket after 20 trials.")

                       return None
                    else:
                       logger.debug("cannot find brent bracket after 20 trials. but very low number of toys")
                       return self.ulSigma ( model, marginalize, 4*toys )
            try:
                mu_lim = optimize.brentq ( root_func, a, b, rtol=1e-03, xtol=1e-06 )
                return mu_lim / self.lumi
            except ValueError as e: ## it could still be that the signs arent opposite
                # in that case, try again
                pass

if __name__ == "__main__":
    C=[ 18774.2, -2866.97, -5807.3, -4460.52, -2777.25, -1572.97, -846.653, -442.531,
       -2866.97, 496.273, 900.195, 667.591, 403.92, 222.614, 116.779, 59.5958,
       -5807.3, 900.195, 1799.56, 1376.77, 854.448, 482.435, 258.92, 134.975,
       -4460.52, 667.591, 1376.77, 1063.03, 664.527, 377.714, 203.967, 106.926,
       -2777.25, 403.92, 854.448, 664.527, 417.837, 238.76, 129.55, 68.2075,
       -1572.97, 222.614, 482.435, 377.714, 238.76, 137.151, 74.7665, 39.5247,
       -846.653, 116.779, 258.92, 203.967, 129.55, 74.7665, 40.9423, 21.7285,
       -442.531, 59.5958, 134.975, 106.926, 68.2075, 39.5247, 21.7285, 11.5732]
    m=Model ( data=[1964,877,354,182,82,36,15,11],
              backgrounds=[2006.4,836.4,350.,147.1,62.0,26.2,11.1,4.7],
              covariance= C,
#              third_moment = [ 0.1, 0.02, 0.1, 0.1, 0.003, 0.0001, 0.0002, 0.0005 ],
              third_moment = [ 0. ] * 8,
              efficiencies=[ x/100. for x in [47,29.4,21.1,14.3,9.4,7.1,4.7,4.3] ],
              name="CMS-NOTE-2017-001 model" )
    ulComp = UpperLimitComputer ( lumi = 1. / fb, ntoys=500, cl=.95 )
    #uls = ulComp.ulSigma ( Model ( 15,17.5,3.2,0.00454755 ) )
    #print ( "uls=", uls )
    ul_old = 131.828*fb ## ulComp.ulSigmaOld ( m )
    print ( "old ul=", ul_old )
    ul = ulComp.ulSigma ( m )
    print ( "ul (marginalized)", ul )
    ul = ulComp.ulSigma ( m, marginalize=False )
    print ( "ul (profiled)", ul )
