#!/usr/bin/env python3

"""
.. module:: simplifiedLikelihoods
   :synopsis: Code that implements the simplified likelihoods as presented
              in CMS-NOTE-2017-001, see https://cds.cern.ch/record/2242860.
              In collaboration with Andy Buckley, Sylvain Fichet and Nickolas Wardle.
              Simplified likelihoods v2 (JHEP_021P_1018) are partly implemented.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Jack Y. Araz <jack.araz@durham.ac.uk>
"""

from scipy import stats, optimize, integrate, special, linalg
from numpy import sqrt, exp, log, sign, array, ndarray
from functools import reduce
from smodels.tools.statistics import CLsfromNLL, determineBrentBracket
from typing import Text, Optional, Union, Tuple

import numpy as np
import math
import copy

try:
    from smodels.tools.smodelsLogging import logger
except ModuleNotFoundError:

    def getLogger():
        """
        Configure the logging facility. Maybe adapted to fit into
        your framework.
        """

        import logging

        logger = logging.getLogger("SL")
        formatter = logging.Formatter("%(module)s - %(levelname)s: %(message)s")
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        # ch.setLevel(logging.DEBUG)
        logger.addHandler(ch)
        return logger

    logger = getLogger()


class Data:
    """A very simple observed container to collect all the data
    needed to fully define a specific statistical model"""

    def __init__(
        self,
        observed,
        backgrounds,
        covariance,
        third_moment=None,
        nsignal=None,
        name="model",
        deltas_rel=0.2,
        lumi=None,
    ):
        """
        :param observed: number of observed events per dataset
        :param backgrounds: expected bg per dataset
        :param covariance: uncertainty in background, as a covariance matrix
        :param nsignal: number of signal events in each dataset
        :param name: give the model a name, just for convenience
        :param deltas_rel: the assumed relative error on the signal hypotheses.
                           The default is 20%.
        :param lumi: luminosity of dataset in 1/fb, or None
        """
        self.observed = self.convert(observed)  # Make sure observed number of events are integers
        ## self.observed = np.around(self.convert(observed)) #Make sure observed number of events are integers
        self.backgrounds = self.convert(backgrounds)
        self.n = len(self.observed)
        self.covariance = self._convertCov(covariance)
        self.nsignal = self.convert(nsignal)
        self.lumi = lumi
        if self.nsignal is None:
            if len(self.backgrounds) == 1:
                # doesnt matter, does it?
                self.nsignal = np.array([1.0])
            self.signal_rel = self.convert(1.0)
        elif self.nsignal.sum():
            self.signal_rel = self.nsignal / self.nsignal.sum()
        else:
            self.signal_rel = array([0.0] * len(self.nsignal))

        self.third_moment = self.convert(third_moment)
        if (
            type(self.third_moment) != type(None)
            and np.sum([abs(x) for x in self.third_moment]) < 1e-10
        ):
            self.third_moment = None
        self.name = name
        self.deltas_rel = deltas_rel
        self._computeABC()

    def totalCovariance(self, nsig):
        """get the total covariance matrix, taking into account
        also signal uncertainty for the signal hypothesis <nsig>.
        If nsig is None, the predefined signal hypothesis is taken.
        """
        if self.isLinear():
            cov_tot = self.V + self.var_s(nsig)
        else:
            cov_tot = self.covariance + self.var_s(nsig)
        return cov_tot

    def zeroSignal(self):
        """
        Is the total number of signal events zero?
        """
        if self.nsignal is None:
            return True
        return len(self.nsignal[self.nsignal > 0.0]) == 0

    def var_s(self, nsig=None):
        """
        The signal variances. Convenience function.

        :param nsig: If None, it will use the model expected number of signal events,
                    otherwise will return the variances for the input value using the relative
                    signal uncertainty defined for the model.

        """

        if nsig is None:
            nsig = self.nsignal
        else:
            nsig = self.convert(nsig)
        return np.diag((nsig * self.deltas_rel) ** 2)

    def isScalar(self, obj):
        """
        Determine if obj is a scalar (float or int)
        """

        if isinstance(obj, ndarray):
            ## need to treat separately since casting array([0.]) to float works
            return False
        try:
            _ = float(obj)
            return True
        except (ValueError, TypeError):
            pass
        return False

    def convert(self, obj):
        """
        Convert object to numpy arrays.
        If object is a float or int, it is converted to a one element
        array.
        """

        if type(obj) == type(None):
            return obj
        if self.isScalar(obj):
            return array([obj])
        return array(obj)

    def __str__(self):
        return self.name + " (%d dims)" % self.n

    def _convertCov(self, obj):

        if self.isScalar(obj):
            return array([[obj]])
        if isinstance(obj[0], list):
            return array(obj)
        if isinstance(obj[0], float):
            ## if the matrix is flattened, unflatten it.
            return array([obj[self.n * i : self.n * (i + 1)] for i in range(self.n)])

        return obj

    def _computeABC(self):
        """
        Compute the terms A, B, C, rho, V. Corresponds with
        Eqs. 1.27-1.30 in the second paper.
        """
        self.V = self.covariance
        if self.third_moment is None:
            self.A = None
            self.B = None
            self.C = None
            return

        covD = self.diagCov()
        C = []
        for m2, m3 in zip(covD, self.third_moment):
            if m3 == 0.0:
                m3 = 1e-30
            k = -np.sign(m3) * sqrt(2.0 * m2)
            dm = sqrt(8.0 * m2**3 / m3**2 - 1.0)
            C.append(k * np.cos(4.0 * np.pi / 3.0 + np.arctan(dm) / 3.0))

        self.C = np.array(C)  ## C, as define in Eq. 1.27 (?) in the second paper
        self.B = sqrt(covD - 2 * self.C**2)  ## B, as defined in Eq. 1.28(?)
        self.A = self.backgrounds - self.C  ## A, Eq. 1.30(?)
        self.rho = np.array([[0.0] * self.n] * self.n)  ## Eq. 1.29 (?)
        for x in range(self.n):
            for y in range(x, self.n):
                bxby = self.B[x] * self.B[y]
                cxcy = self.C[x] * self.C[y]
                e = (4.0 * cxcy) ** (-1) * (
                    sqrt(bxby**2 + 8 * cxcy * self.covariance[x][y]) - bxby
                )
                self.rho[x][y] = e
                self.rho[y][x] = e

        self.sandwich()
        # self.V = sandwich ( self.B, self.rho )

    def sandwich(self):
        """
        Sandwich product
        """

        ret = np.array([[0.0] * len(self.B)] * len(self.B))
        for x in range(len(self.B)):
            for y in range(x, len(self.B)):
                T = self.B[x] * self.B[y] * self.rho[x][y]
                ret[x][y] = T
                ret[y][x] = T
        self.V = ret

    def isLinear(self):
        """
        Statistical model is linear, i.e. no quadratic term in poissonians
        """

        return type(self.C) == type(None)

    def diagCov(self):
        """
        Diagonal elements of covariance matrix. Convenience function.
        """

        return np.diag(self.covariance)

    def correlations(self):
        """
        Correlation matrix, computed from covariance matrix.
        Convenience function.
        """

        if hasattr(self, "corr"):
            return self.corr

        self.corr = copy.deepcopy(self.covariance)
        for x in range(self.n):
            self.corr[x][x] = 1.0
            for y in range(x + 1, self.n):
                rho = self.corr[x][y] / sqrt(self.covariance[x][x] * self.covariance[y][y])
                self.corr[x][y] = rho
                self.corr[y][x] = rho
        return self.corr

    def rel_signals(self, mu):
        """
        Returns the number of expected relative signal events, for all datasets,
        given total signal strength mu. For mu=1, the sum of the numbers = 1.

        :param mu: Total number of signal events summed over all datasets.
        """

        return mu * self.signal_rel

    def nsignals(self, mu):
        """
        Returns the number of expected signal events, for all datasets,
        given total signal strength mu.

        :param mu: Total number of signal events summed over all datasets.
        """

        return mu * self.nsignal


class LikelihoodComputer:

    debug_mode = False

    def __init__(self, data, toys=30000):
        """
        :param data: a Data object.
        :param toys: number of toys when marginalizing
        """

        self.model = data
        self.toys = toys

    def dNLLdMu(self, mu, theta_hat = None ):
        """
        d (- ln L)/d mu, if L is the likelihood. The function
        whose root gives us muhat, i.e. the mu that maximizes
        the likelihood.

        :param mu: total number of signal events
        :param theta_hat: array with nuisance parameters, if None then
                          compute them

        """
        if type(mu) in [ list, np.ndarray ] and len(mu)==1:
            mu = float(mu[0])
        if theta_hat == None:
            theta_hat, _ = self.findThetaHat ( mu )
        nsig = self.model.nsignal
        if not self.model.isLinear():
            logger.debug("implemented only for linear model")
        # n_pred^i := mu s_i + b_i + theta_i
        # NLL = sum_i [ - n_obs^i * ln ( n_pred^i ) + n_pred^i ]
        # d NLL / d mu = sum_i [ - ( n_obs^i * s_ i ) / n_pred_i + s_i ]

        # Define relative signal strengths:
        n_pred = mu * nsig + self.model.backgrounds + theta_hat

        for ctr, d in enumerate(n_pred):
            if d == 0.0:
                if (self.model.observed[ctr] * nsig[ctr]) == 0.0:
                    #    logger.debug("zero denominator, but numerator also zero, so we set denom to 1.")
                    n_pred[ctr] = 1e-5
                else:
                    # n_pred[ctr]=1e-5
                    raise Exception(
                        "we have a zero value in the denominator at pos %d, with a non-zero numerator. dont know how to handle."
                        % ctr
                    )
        ret = -self.model.observed * nsig / n_pred + nsig

        if type(ret) in [array, ndarray, list]:
            ret = sum(ret)
        return ret

    def extendedOutput(self, extended_output, default=None):
        if extended_output:
            ret = {"muhat": default, "sigma_mu": default, "lmax": default}
            return ret
        return default

    def findAvgr(self, theta_hat ):
        """from the difference observed - background, find got inital
        values for lower and upper"""
        mu_c = self.model.observed - self.model.backgrounds - theta_hat
        mu_r, wmu_r = [], []
        hessian = self.d2NLLdMu2 ( 1., theta_hat )
        wtot = 0.0
        for s in zip(mu_c, self.model.nsignal, hessian):
            if s[1] > 1e-10:
                w = 1. # 1e-5
                if s[2] > 0.0:
                    w = s[2]
                wtot += w
                r = s[0] / s[1]
                mu_r.append( r )
                wmu_r.append(w * r )
        ret = min(mu_r), sum(wmu_r) / wtot, max(mu_r)
        return ret

    def d2NLLdMu2 ( self, mu, theta_hat, allowZeroHessian=True ):
        """ the hessian of the likelihood of mu, at mu,
        which is the Fisher information
        which is approximately the inverse of the covariance
        :param allowZeroHessian: if false and sum(observed)==0, then replace
                                 observed with expected
        """
        # nll=-nobs*ln(mu*s + b + theta) + ( mu*s + b + theta)
        # d nll / d mu = - nobs * s / ( mu*s + b + theta) + s
        # d2nll / dmu2 = nobs * s**2 / ( mu*s + b + theta )**2
        n_pred = mu * self.model.nsignal + self.model.backgrounds + theta_hat
        for i, s in enumerate(n_pred):
            if s == 0.0: # the denominator in the hessian is 0?
                if (self.model.observed[i] * self.model.nsignal[i]) == 0.0:
                    #    logger.debug("zero denominator, but numerator also zero, so we set denom to 1.")
                    n_pred[i] = 1.0
                else:
                    raise Exception( f"we have a zero value in the denominator at pos {i}, with a non-zero numerator. dont know how to handle." )
        obs = self.model.observed
        if sum(obs)==0 and not allowZeroHessian:
            obs = self.model.backgrounds
        hessian = obs * self.model.nsignal**2 / n_pred**2
        if sum(hessian) == 0.0 and not allowZeroHessian:
            # if all observations are zero, we replace them by the expectations
            if sum(self.model.observed) == 0:
                hessian = self.model.nsignal**2 / n_pred
        return hessian

    #def findMuHat(
    def findMuHatViaBracketing( self, allowNegativeSignals=False,
        extended_output=False, nll=False, marginalize=False):
        """
        Find the most likely signal strength mu via a brent bracketing technique
        given the relative signal strengths in each dataset (signal region).

        :param allowNegativeSignals: if true, then also allow for negative values
        :param extended_output: if true, return also sigma_mu, the estimate of the error of mu_hat,
         and lmax, the likelihood at mu_hat
        :param nll: if true, return nll instead of lmax in the extended output

        :returns: mu_hat, i.e. the maximum likelihood estimate of mu, if extended output is
        requested, it returns a dictionary with mu_hat, sigma_mu -- the standard deviation around mu_hat, and lmax, i.e. the likelihood at mu_hat
        """
        if (self.model.backgrounds == self.model.observed).all():
            return self.extendedOutput(extended_output, 0.0)
        nsig = self.model.nsignal

        if type(nsig) in [list, ndarray]:
            nsig = array(nsig)

        nsig[nsig == 0.0] = 1e-20
        if sum(nsig < 0.0):
            raise Exception("Negative relative signal strengths!")

        ## we need a very rough initial guess for mu(hat), to come
        ## up with a first theta
        # self.nsig = array([0.]*len(self.model.observed))
        self.mu = 1.
        ## we start with theta_hat being all zeroes
        # theta_hat = array([0.]*len(self.model.observed))
        mu_hat_old, mu_hat = 0.0, 1.0
        ctr = 0
        widener = 3.0
        while (
            abs(mu_hat - mu_hat_old) > 1e-10
            and abs(mu_hat - mu_hat_old) / (mu_hat + mu_hat_old) > 0.5e-2
            and ctr < 20
        ):
            theta_hat, _ = self.findThetaHat(mu_hat)
            ctr += 1
            mu_hat_old = mu_hat
            minr, avgr, maxr = self.findAvgr( theta_hat )
            # for i,s in enumerate ( signal_rel ):
            #    if abs(s) < 1e-19:
            #        mu_c[i]=0.
            ## find mu_hat by finding the root of 1/L dL/dmu. We know
            ## that the zero has to be between min(mu_c) and max(mu_c).
            lstarters = [ avgr - 0.2 * abs(avgr), minr, 0.0, -1.0, 1.0, 10.0, -0.1, \
                          0.1, -100.0, 100.0, -1000.0, ]
            closestl, closestr = None, float("inf")
            for lower in lstarters:
                lower_v = self.dNLLdMu(lower, theta_hat)
                if lower_v < 0.0:
                    break
                if lower_v < closestr:
                    closestl, closestr = lower, lower_v
            if lower_v > 0.0:
                logger.debug(
                    f"did not find a lower value with rootfinder(lower) < 0. Closest: f({closestl})={closestr}"
                )
                return self.extendedOutput(extended_output, 0.0)
            ustarters = [ avgr + 0.2 * abs(avgr), maxr, 0.0, 1.0, 10.0, -1.0 - 0.1,\
                0.1, 100.0, -100.0, 1000.0, -1000.0, 0.01, -0.01, ]
            closestl, closestr = None, float("inf")
            for upper in ustarters:
                upper_v = self.dNLLdMu(upper, theta_hat)
                if upper_v > 0.0:
                    break
                if upper_v < closestr:
                    closestl, closestr = upper, upper_v
            if upper_v < 0.0:
                logger.debug("did not find an upper value with rootfinder(upper) > 0.")
                return self.extendedOutput(extended_output, 0.0)
            mu_hat = optimize.brentq(self.dNLLdMu, lower, upper, args=(theta_hat, ), rtol=1e-9)
            if not allowNegativeSignals and mu_hat < 0.0:
                mu_hat = 0.0
                theta_hat, _ = self.findThetaHat(mu_hat)
            self.theta_hat = theta_hat

        if extended_output:
            sigma_mu = self.getSigmaMu(mu_hat, theta_hat)
            llhd = self.likelihood( mu_hat, marginalize=marginalize, nll=nll)
            # print ( f"returning {allowNegativeSignals}: mu_hat {mu_hat}+-{sigma_mu} llhd {llhd}" )
            ret = {"muhat": mu_hat, "sigma_mu": sigma_mu, "lmax": llhd}
            return ret
        return mu_hat

    def getSigmaMu(self, mu, theta_hat):
        """
        Get an estimate for the standard deviation of mu at <mu>, from
        the inverse hessian
        """
        if not self.model.isLinear():
            logger.debug("implemented only for linear model")
        # d^2 mu NLL / d mu^2 = sum_i [ n_obs^i * s_i**2 / n_pred^i**2 ]
        hessian = self.d2NLLdMu2 ( mu, theta_hat, allowZeroHessian = False )
        hessian = sum ( hessian )
        """
        if hessian == 0.0:
            # if all observations are zero, we replace them by the expectations
            if sum(self.model.observed) == 0:
                hessian = sum(nsig**2 / n_pred)
        """
        stderr = float(np.sqrt(1.0 / hessian))
        return stderr

    # Define integrand (gaussian_(bg+signal)*poisson(nobs)):
    # def prob(x0, x1 )
    def llhdOfTheta(self, theta, nll = True ):
        """ likelihood for nuicance parameters theta, given signal strength
            self.mu. notice, by default it returns nll
        :param theta: nuisance parameters
        :params nll: if True, compute negative log likelihood
        """
        # theta = array ( thetaA )
        # ntot = self.model.backgrounds + self.nsig
        nsig = self.mu * self.model.nsignal
        if self.model.isLinear():
            lmbda = self.model.backgrounds + nsig + theta
        else:
            lmbda = nsig + self.model.A + theta + self.model.C * theta**2 / self.model.B**2
        lmbda[lmbda <= 0.0] = 1e-30  ## turn zeroes to small values
        obs = self.model.observed

        def is_integer(x):
            if type(x) in [int, np.int64]:
                return True
            if type(x) in [float]:
                return x.is_integer()
            return False

        ## not needed for now
        allintegers = np.all([is_integer(i) for i in obs])
        if nll:
            if allintegers:
                poisson = stats.poisson.logpmf(obs, lmbda)
            else:
                poisson = -lmbda + obs * np.log(lmbda) - special.loggamma(obs + 1)
        else:
            if allintegers:
                poisson = stats.poisson.pmf(obs, lmbda)
            else:
                # poisson = np.exp(-lmbda)*np.power(lmbda,obs)/special.gamma(obs+1)
                logpoiss = -lmbda + obs * np.log(lmbda) - special.loggamma(obs + 1)
                poisson = np.exp(logpoiss)
        try:
            M = [0.0] * len(theta)
            C = self.model.V
            # if self.model.n == 1: I think not a good idea
            #    C = self.model.totalCovariance(self.nsig)
            dTheta = theta - M
            expon = -0.5 * np.dot(np.dot(dTheta, self.weight), dTheta) + self.logcoeff
            # print ( "expon", expon, "coeff", self.coeff )
            if nll:
                gaussian = expon  #  np.log ( self.coeff )
                # gaussian2 = stats.multivariate_normal.logpdf(theta,mean=M,cov=C)
                ret = -gaussian - sum(poisson)
            else:
                gaussian = np.exp(expon)
                # gaussian = self.coeff * np.exp ( expon )
                # gaussian2 = stats.multivariate_normal.pdf(theta,mean=M,cov=C)
                ret = gaussian * (reduce(lambda x, y: x * y, poisson))
            return ret
        except ValueError as e:
            raise Exception("ValueError %s, %s" % (e, self.model.V))
            # raise Exception("ValueError %s, %s" % ( e, self.model.totalCovariance(self.nsig) ))
            # raise Exception("ValueError %s, %s" % ( e, self.model.V ))

    def dNLLdTheta(self, theta):
        """the derivative of nll as a function of the thetas.
        Makes it easier to find the maximum likelihood."""
        # print ( f"nsig {self.nsig} {self.model.nsignal}" )
        nsig = self.mu * self.model.nsignal
        if self.model.isLinear():
            xtot = theta + self.model.backgrounds + nsig
            xtot[xtot <= 0.0] = 1e-30  ## turn zeroes to small values
            nllp_ = self.ones - self.model.observed / xtot + np.dot(theta, self.weight)
            return nllp_
        lmbda = nsig + self.model.A + theta + self.model.C * theta**2 / self.model.B**2
        lmbda[lmbda <= 0.0] = 1e-30  ## turn zeroes to small values
        # nllp_ = ( self.ones - self.model.observed / lmbda + np.dot( theta , self.weight ) ) * ( self.ones + 2*self.model.C * theta / self.model.B**2 )
        T = self.ones + 2 * self.model.C / self.model.B**2 * theta
        nllp_ = T - self.model.observed / lmbda * (T) + np.dot(theta, self.weight)
        return nllp_

    def d2NLLdTheta2(self, theta):
        """the Hessian of nll as a function of the thetas.
        Makes it easier to find the maximum likelihood."""
        # xtot = theta + self.ntot
        nsig = self.mu * self.model.nsignal
        if self.model.isLinear():
            xtot = theta + self.model.backgrounds + nsig
            xtot[xtot <= 0.0] = 1e-30  ## turn zeroes to small values
            nllh_ = self.weight + np.diag(self.model.observed / (xtot**2))
            return nllh_
        lmbda = nsig + self.model.A + theta + self.model.C * theta**2 / self.model.B**2
        lmbda[lmbda <= 0.0] = 1e-30  ## turn zeroes to small values
        T = self.ones + 2 * self.model.C / self.model.B**2 * theta
        # T_i = 1_i + 2*c_i/B_i**2 * theta_i
        nllh_ = (
            self.weight
            + np.diag(self.model.observed * T**2 / (lmbda**2))
            - np.diag(self.model.observed / lmbda * 2 * self.model.C / self.model.B**2)
            + np.diag(2 * self.model.C / self.model.B**2)
        )
        return nllh_

    def getThetaHat(self, nobs, nb, mu, covb, max_iterations):
        """ Compute nuisance parameter theta that
            maximizes our likelihood (poisson*gauss) -- by setting dNLL/dTheta
            to zero
        :param mu: signal strength
        :returns: theta_hat
        """
        # ntot = mu * nsig + nb
        # nll = - nobs ln(ntot + theta) - ntot - theta - theta**2/(2 delta**2)
        # dnll/dtheta = - nobs / ( ntot + theta ) + 1 + theta / delta**2
        # theta**2 + theta * ( delta**2 + ntot ) + delta**2 * ( ntot-nobs) = 0
        nsig = mu * self.model.nsignal
        self.mu = mu
        sigma2 = covb + self.model.var_s(nsig)  ## np.diag ( (self.model.deltas)**2 )
        ## for now deal with variances only
        ntot = nb + nsig
        cov = np.array(sigma2)
        # weight = cov**(-1) ## weight matrix
        weight = linalg.inv(cov)
        diag_cov = np.diag(cov)
        # first: no covariances:
        q = diag_cov * (ntot - nobs)
        p = ntot + diag_cov
        thetamaxes = []
        # thetamax = -p / 2.0 * (1 - sign(p) * sqrt(1.0 - 4 * q / p**2))
        thetamax = -p / 2.0 + sign(p) * sqrt( p**2/4 - q)
        thetamaxes.append(thetamax)
        ndims = len(p)

        def distance(theta1, theta2):
            for ctr, i in enumerate(theta1):
                if i == 0.0:
                    theta1[ctr] = 1e-20
            for ctr, i in enumerate(theta2):
                if i == 0.0:
                    theta2[ctr] = 1e-20
            return sum(np.abs(theta1 - theta2) / np.abs(theta1 + theta2))

        ictr = 0
        while ictr < max_iterations:
            ictr += 1
            q = diag_cov * (ntot - nobs)
            p = ntot + diag_cov
            for i in range(ndims):
                # q[i] = diag_cov[i] * ( ntot[i] - nobs[i] )
                # p[i] = ntot[i] + diag_cov[i]
                for j in range(ndims):
                    if i == j:
                        continue
                    dq = thetamax[j] * ntot[i] * diag_cov[i] * weight[i, j]
                    dp = thetamax[j] * weight[i, j] * diag_cov[i]
                    if abs(dq / q[i]) > 0.3:
                        # logger.warning ( "too big a change in iteration." )
                        dq = np.abs(0.3 * q[i]) * np.sign(dq)
                    if abs(dp / p[i]) > 0.3:
                        # logger.warning ( "too big a change in iteration." )
                        dp = np.abs(0.3 * p[i]) * np.sign(dp)
                    q[i] += dq
                    p[i] += dp
                # thetamax = -p / 2.0 * (1 - sign(p) * sqrt(1.0 - 4 * q / p**2))
                thetamax = -p / 2.0 + sign(p) * sqrt( p**2/4 - q)
            thetamaxes.append(thetamax)
            if len(thetamaxes) > 2:
                d1 = distance(thetamaxes[-2], thetamax)
                d2 = distance(thetamaxes[-3], thetamaxes[-2])
                if d1 > d2:
                    raise Exception("diverging when computing thetamax: %f > %f" % (d1, d2))
                if d1 < 1e-5:
                    return thetamax
        return thetamax

    def findThetaHat(self, mu ):
        """Compute nuisance parameters theta that maximize our likelihood
        (poisson*gauss).
        """
        nsig = mu * self.model.nsignal

        ## first step is to disregard the covariances and solve the
        ## quadratic equations
        ini = self.getThetaHat(
            self.model.observed, self.model.backgrounds, mu, self.model.covariance, 0
        )
        self.cov_tot = self.model.V
        # if self.model.n == 1:
        #    self.cov_tot = self.model.totalCovariance ( nsig )
        # if not self.model.isLinear():
        # self.cov_tot = self.model.V + self.model.var_s(nsig)
        # self.cov_tot = self.model.totalCovariance (nsig)
        self.weight = np.linalg.inv(self.cov_tot)
        # self.coeff = 1.
        logdet = np.linalg.slogdet(self.cov_tot)
        self.logcoeff = -self.model.n / 2 * np.log(2 * np.pi) - 0.5 * logdet[1]
        # self.coeff = (2*np.pi)**(-self.model.n/2) * np.exp(-.5* logdet[1] )
        # print ( "coeff", self.coeff, "n", self.model.n, "det", np.linalg.slogdet ( self.cov_tot ) )
        # print ( "cov_tot", self.cov_tot[:10] )
        self.ones = 1.0
        if type(self.model.observed) in [list, ndarray]:
            self.ones = np.ones(len(self.model.observed))
        self.gammaln = special.gammaln(self.model.observed + 1)
        try:
            ret_c = optimize.fmin_ncg(
                self.llhdOfTheta,
                ini,
                fprime=self.dNLLdTheta,
                fhess=self.d2NLLdTheta2,
                full_output=True,
                disp=0,
            )
            # then always continue with TNC
            if type(self.model.observed) in [int, float]:
                bounds = [(-10 * self.model.observed, 10 * self.model.observed)]
            else:
                bounds = [(-10 * x, 10 * x) for x in self.model.observed]
            ini = ret_c
            ret_c = optimize.fmin_tnc(
                self.llhdOfTheta, ret_c[0], fprime=self.dNLLdTheta, disp=0, bounds=bounds
            )
            # print ( "[findThetaHat] mu=%s bg=%s observed=%s V=%s, nsig=%s theta=%s, nll=%s" % ( self.nsig[0]/self.model.efficiencies[0], self.model.backgrounds, self.model.observed,self.model.covariance, self.nsig, ret_c[0], self.nllOfNuisances(ret_c[0]) ) )
            if ret_c[-1] not in [0, 1, 2]:
                return ret_c[0], ret_c[-1]
            else:
                return ret_c[0], 0
                logger.debug("tnc worked.")

            ret = ret_c[0]
            return ret, -2
        except (IndexError, ValueError) as e:
            logger.error("exception: %s. ini[-3:]=%s" % (e, ini[-3:]))
            raise Exception("cov-1=%s" % (self.model.covariance + self.model.var_s(nsig)) ** (-1))
        return ini, -1

    def marginalizedLLHD1D(self, mu, nll):
        """
        Return the likelihood (of 1 signal region) to observe nobs events given the
        predicted background nb, error on this background (deltab),
        signal strength of mu and the relative error on the signal (deltas_rel).

        :param mu: predicted signal strength (float)
        :param nobs: number of observed events (float)
        :param nb: predicted background (float)
        :param deltab: uncertainty on background (float)

        :return: likelihood to observe nobs events (float)

        """
        nsig = self.model.nsignal * mu
        self.sigma2 = self.model.covariance + self.model.var_s(nsig)  ## (self.model.deltas)**2
        self.sigma_tot = sqrt(self.sigma2)
        self.lngamma = math.lgamma(self.model.observed[0] + 1)
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

        # Define integrand (gaussian_(bg+signal)*poisson(nobs)):
        def prob(x, nsig):
            poisson = exp(self.model.observed * log(x) - x - self.lngamma)
            gaussian = stats.norm.pdf(x, loc=self.model.backgrounds + nsig, scale=self.sigma_tot)

            return poisson * gaussian

        # Compute maximum value for the integrand:
        xm = self.model.backgrounds + nsig - self.sigma2
        # If nb + nsig = sigma2, shift the values slightly:
        if xm == 0.0:
            xm = 0.001
        xmax = (
            xm
            * (1.0 + sign(xm) * sqrt(1.0 + 4.0 * self.model.observed * self.sigma2 / xm**2))
            / 2.0
        )

        # Define initial integration range:
        nrange = 5.0
        a = max(0.0, xmax - nrange * sqrt(self.sigma2))
        b = xmax + nrange * self.sigma_tot
        like = integrate.quad(prob, a, b, (nsig), epsabs=0.0, epsrel=1e-3)[0]
        if like == 0.0:
            return 0.0

        # Increase integration range until integral converges
        err = 1.0
        ctr = 0
        while err > 0.01:
            ctr += 1
            if ctr > 10.0:
                raise Exception("Could not compute likelihood within required precision")

            like_old = like
            nrange = nrange * 2
            a = max(0.0, (xmax - nrange * self.sigma_tot)[0][0])
            b = (xmax + nrange * self.sigma_tot)[0][0]
            like = integrate.quad(prob, a, b, (nsig), epsabs=0.0, epsrel=1e-3)[0]
            if like == 0.0:
                continue
            err = abs(like_old - like) / like

        # Renormalize the likelihood to account for the cut at x = 0.
        # The integral of the gaussian from 0 to infinity gives:
        # (1/2)*(1 + Erf(mu/sqrt(2*sigma2))), so we need to divide by it
        # (for mu - sigma >> 0, the normalization gives 1.)
        norm = (1.0 / 2.0) * (
            1.0 + special.erf((self.model.backgrounds + nsig) / sqrt(2.0 * self.sigma2))
        )
        like = like / norm

        if nll:
            like = -log(like)

        return like[0][0]

    def marginalizedLikelihood(self, mu, nll):
        """compute the marginalized likelihood of observing nsig signal event"""
        if (
            self.model.isLinear() and self.model.n == 1
        ):  ## 1-dimensional non-skewed llhds we can integrate analytically
            return self.marginalizedLLHD1D(mu, nll)
        nsig = mu * self.model.nsignal

        vals = []
        self.gammaln = special.gammaln(self.model.observed + 1)
        thetas = stats.multivariate_normal.rvs(
            mean=[0.0] * self.model.n,
            # cov=(self.model.totalCovariance(nsig)),
            cov=self.model.V,
            size=self.toys,
        )  ## get ntoys values
        for theta in thetas:
            if self.model.isLinear():
                lmbda = nsig + self.model.backgrounds + theta
            else:
                lmbda = nsig + self.model.A + theta + self.model.C * theta**2 / self.model.B**2
            if self.model.isScalar(lmbda):
                lmbda = array([lmbda])
            for ctr, v in enumerate(lmbda):
                if v <= 0.0:
                    lmbda[ctr] = 1e-30
                # print ( "lmbda=",lmbda )
            poisson = self.model.observed * np.log(lmbda) - lmbda - self.gammaln
            # poisson = np.exp(self.model.observed*np.log(lmbda) - lmbda - self.model.backgrounds - self.gammaln)
            vals.append(np.exp(sum(poisson)))
            # vals.append ( reduce(lambda x, y: x*y, poisson) )
        mean = np.mean(vals)
        if nll:
            if mean == 0.0:
                mean = 1e-100
            mean = -log(mean)
        return mean

    def profileLikelihood(self, mu, nll ):
        """compute the profiled likelihood for mu.
        Warning: not normalized.
        Returns profile likelihood and error code (0=no error)
        """
        # compute the profiled (not normalized) likelihood of observing
        # nsig signal events
        theta_hat, _ = self.findThetaHat(mu)
        if self.debug_mode:
            self.theta_hat = theta_hat
        ret = self.llhdOfTheta( theta_hat, nll )

        return ret

    def likelihood(self, mu: float, marginalize: bool = False, nll: bool = False) -> float:
        """compute likelihood for mu, profiling or marginalizing the nuisances
        :param mu: float Parameter of interest, signal strength
        :param marginalize: if true, marginalize, if false, profile
        :param nll: return nll instead of likelihood
        """
        # print ( f"likelihood {nsig[:2]} {self.model.nsignal[:2]} mu={mu}" )
        if marginalize:
            # p,err = self.profileLikelihood ( nsig, deltas )
            return self.marginalizedLikelihood(mu, nll)
            # print ( "p,l=",p,l,p/l )
        else:
            return self.profileLikelihood(mu, nll)

    def lmax(self, marginalize=False, nll=False, allowNegativeSignals=False):
        """convenience function, computes likelihood for nsig = nobs-nbg,
        :param marginalize: if true, marginalize, if false, profile nuisances.
        :param nll: return nll instead of likelihood
        :param allowNegativeSignals: if False, then negative nsigs are replaced with 0.
        """
        if len(self.model.observed) == 1:
            dn = self.model.observed - self.model.backgrounds
            if not allowNegativeSignals and dn[0] < 0.0:
                dn = [0.0]
            self.muhat = float(dn[0])
            if abs(self.model.nsignal) > 1e-100:
                self.muhat = float(dn[0] / self.model.nsignal[0])
            self.sigma_mu = np.sqrt(self.model.observed[0] + self.model.covariance[0][0])
            return self.likelihood(marginalize=marginalize, nll=nll, mu = self.muhat )
        fmh = self.findMuHat( allowNegativeSignals=allowNegativeSignals,
                              extended_output=True, nll=nll
        )
        muhat_, sigma_mu, lmax = fmh["muhat"], fmh["sigma_mu"], fmh["lmax"]
        self.muhat = muhat_
        self.sigma_mu = sigma_mu
        return self.likelihood ( marginalize=marginalize, nll=nll, mu=muhat_ )

    def findMuHat(
    #def findMuHatViaGradientDescent(
        self,
        allowNegativeSignals=False,
        extended_output=False,
        nll=False,
        marginalize=False,
    ):
        """
        Find the most likely signal strength mu via gradient descent
        given the relative signal strengths in each dataset (signal region).

        :param allowNegativeSignals: if true, then also allow for negative values
        :param extended_output: if true, return also sigma_mu, the estimate of the error of mu_hat,
         and lmax, the likelihood at mu_hat
        :param nll: if true, return nll instead of lmax in the extended output

        :returns: mu_hat, i.e. the maximum likelihood estimate of mu, if extended output is
        requested, it returns mu_hat, sigma_mu -- the standard deviation around mu_hat, and llhd,
        the likelihood at mu_hat
        """
        theta_hat, _ = self.findThetaHat( 0. )
        minr, avgr, maxr = self.findAvgr( theta_hat )
        theta_hat, _ = self.findThetaHat( avgr )
        minr, avgr, maxr = self.findAvgr( theta_hat )

        def myllhd(mu):
            theta = self.findThetaHat ( mu=mu )
            ret = self.likelihood(nll=True, marginalize=marginalize, mu = mu )
            return ret

        import scipy.optimize
        ominr = minr
        if minr > 0.:
            minr = .5 * minr
        if minr < 0.:
            minr = 2.*minr
        if maxr > 0.:
            maxr = 3.*maxr + 1e-5
        if maxr <= 0.:
            maxr = .3 * maxr + 1e-5

        bounds = [(minr,maxr)]
        if not allowNegativeSignals:
            bounds = [(0, max(maxr,1e-5))]
        assert bounds[0][1] > bounds[0][0], f"bounds are in wrong order: {bounds}"
        o = scipy.optimize.minimize( myllhd, x0=avgr, bounds=bounds, jac = self.dNLLdMu )
        llhd = o.fun
        if not nll:
            llhd = np.exp(-o.fun)
        """
        hess = o.hess_inv
        try:
            hess = hess.todense()
        except Exception as e:
            pass
        """
        mu_hat = float(o.x[0])
        if extended_output:
            sigma_mu = self.getSigmaMu ( mu_hat, theta_hat )
            llhd = self.likelihood( mu_hat, marginalize=marginalize, nll=nll)
            # sigma_mu = float(np.sqrt(hess[0][0]))
            ret = {"muhat": mu_hat, "sigma_mu": sigma_mu, "lmax": llhd }
            return ret
        return mu_hat

    def chi2(self, marginalize=False):
        """
        Computes the chi2 for a given number of observed events nobs given
        the predicted background nb, error on this background deltab,
        expected number of signal events nsig and the relative error on
        signal (deltas_rel).
        :param marginalize: if true, marginalize, if false, profile
        :param nsig: number of signal events
        :return: chi2 (float)

        """

        # Compute the likelhood for the null hypothesis (signal hypothesis) H0:
        llhd = self.likelihood(1., marginalize=marginalize, nll=True)

        # Compute the maximum likelihood H1, which sits at nsig = nobs - nb
        # (keeping the same % error on signal):
        if len(self.model.observed) == 1:
            # TODO this nsig initiation seems wrong and changing maxllhd to likelihood
            # fails ./testStatistics.py : zero division error in L115
            mu_hat = ( self.model.observed - self.model.backgrounds ) / self.model.nsignal
            maxllhd = self.likelihood (mu_hat, marginalize=marginalize, nll=True )
        else:
            maxllhd = self.lmax( marginalize=marginalize, nll=True, allowNegativeSignals=False)
        chi2 = 2 * (llhd - maxllhd)

        if not np.isfinite(chi2):
            logger.error("chi2 is not a finite number! %s,%s,%s" % (chi2, llhd, maxllhd))
        # Return the test statistic -2log(H0/H1)
        return chi2


class UpperLimitComputer:
    debug_mode = False

    def __init__(self, ntoys: float = 30000, cl: float = 0.95):

        """
        :param ntoys: number of toys when marginalizing
        :param cl: desired quantile for limits
        """
        self.toys = ntoys
        self.cl = cl

    def getUpperLimitOnSigmaTimesEff(
        self, model, marginalize=False, toys=None, expected=False, trylasttime=False
    ):
        """upper limit on the fiducial cross section sigma times efficiency,
            summed over all signal regions, i.e. sum_i xsec^prod_i eff_i
            obtained from the defined Data (using the signal prediction
            for each signal regio/dataset), by using
            the q_mu test statistic from the CCGV paper (arXiv:1007.1727).

        :params marginalize: if true, marginalize nuisances, else profile them
        :params toys: specify number of toys. Use default is none
        :params expected: if false, compute observed,
                          true: compute a priori expected, "posteriori":
                          compute a posteriori expected
        :params trylasttime: if True, then dont try extra
        :returns: upper limit on fiducial cross section
        """
        ul = self.getUpperLimitOnMu(
            model, marginalize=marginalize, toys=toys, expected=expected,
            trylasttime=trylasttime)

        if ul == None:
            return ul
        if model.lumi is None:
            logger.error(f"asked for upper limit on fiducial xsec, but no lumi given with the data")
            return ul
        xsec = sum(model.nsignal) / model.lumi
        return ul * xsec

    def getCLsRootFunc(
        self,
        model: Data,
        marginalize: Optional[bool] = False,
        toys: Optional[float] = None,
        expected: Optional[Union[bool, Text]] = False,
        trylasttime: Optional[bool] = False,
    ) -> Tuple:
        """
        Obtain the function "CLs-alpha[0.05]" whose root defines the upper limit,
        plus mu_hat and sigma_mu
        :param model: statistical model
        :param marginalize: if true, marginalize nuisances, else profile them
        :param toys: specify number of toys. Use default is none
        :param expected: if false, compute observed,
                          true: compute a priori expected, "posteriori":
                          compute a posteriori expected
        :param trylasttime: if True, then dont try extra
        :return: mu_hat, sigma_mu, CLs-alpha
        """
        # if expected:
        #    marginalize = True
        if model.zeroSignal():
            """only zeroes in efficiencies? cannot give a limit!"""
            return None, None, None
        if toys == None:
            toys = self.toys
        oldmodel = model
        if expected:
            model = copy.deepcopy(oldmodel)
            if expected == "posteriori":
                tempc = LikelihoodComputer(oldmodel, toys)
                theta_hat_, _ = tempc.findThetaHat(0 * oldmodel.nsignal )
            for i, d in enumerate(model.backgrounds):
                if expected == "posteriori":
                    d += theta_hat_[i]
                model.observed[i] = float(d)
        computer = LikelihoodComputer(model, toys)
        mu_hat = computer.findMuHat( allowNegativeSignals=False, extended_output=False)
        theta_hat0, _ = computer.findThetaHat(0 * model.nsignal )
        sigma_mu = computer.getSigmaMu(mu_hat, theta_hat0)

        nll0 = computer.likelihood( mu_hat, marginalize=marginalize, nll=True)
        # print ( f"SL nll0 {nll0:.3f} muhat {mu_hat:.3f} sigma_mu {sigma_mu:.3f} {signal_type} {sum(model.nsignal):.3f}" )
        if np.isinf(nll0) and not marginalize and not trylasttime:
            logger.warning(
                "nll is infinite in profiling! we switch to marginalization, but only for this one!"
            )
            marginalize = True
            # TODO convert rel_signals to signals
            nll0 = computer.likelihood( mu = mu_hat, marginalize=True, nll=True)
            if np.isinf(nll0):
                logger.warning("marginalization didnt help either. switch back.")
                marginalize = False
            else:
                logger.warning("marginalization worked.")
        aModel = copy.deepcopy(model)
        aModel.observed = array([x + y for x, y in zip(model.backgrounds, theta_hat0)])
        aModel.name = aModel.name + "A"
        # print ( f"SL finding mu hat with {aModel.signal_rel}: mu_hatA, obs: {aModel.observed}" )
        compA = LikelihoodComputer(aModel, toys)
        ## compute
        mu_hatA = compA.findMuHat()
        # TODO convert rel_signals to signals
        nll0A = compA.likelihood( mu=mu_hatA, marginalize=marginalize, nll=True)
        # return 1.

        def clsRoot(mu: float, return_type: Text = "CLs-alpha") -> float:
            """
            Calculate the root
            :param mu: float POI
            :param return_type: (Text) can be "CLs-alpha", "1-CLs", "CLs"
                        CLs-alpha: returns CLs - 0.05
                        1-CLs: returns 1-CLs value
                        CLs: returns CLs value
            """
            nll = computer.likelihood(mu, marginalize=marginalize, nll=True)
            nllA = compA.likelihood(mu, marginalize=marginalize, nll=True)
            return CLsfromNLL(nllA, nll0A, nll, nll0, return_type=return_type)

        return mu_hat, sigma_mu, clsRoot

    def getUpperLimitOnMu(
        self, model, marginalize=False, toys=None, expected=False, trylasttime=False
    ):
        """upper limit on the signal strength multiplier mu
            obtained from the defined Data (using the signal prediction

            for each signal regio/dataset), by using
            the q_mu test statistic from the CCGV paper (arXiv:1007.1727).

        :params marginalize: if true, marginalize nuisances, else profile them
        :params toys: specify number of toys. Use default is none
        :params expected: if false, compute observed,
                          true: compute a priori expected, "posteriori":
                          compute a posteriori expected
        :params trylasttime: if True, then dont try extra
        :returns: upper limit on the signal strength multiplier mu
        """
        mu_hat, sigma_mu, clsRoot = self.getCLsRootFunc(
            model, marginalize, toys, expected, trylasttime
        )
        if mu_hat == None:
            return None
        a, b = determineBrentBracket(mu_hat, sigma_mu, clsRoot)
        mu_lim = optimize.brentq(clsRoot, a, b, rtol=1e-03, xtol=1e-06)
        return mu_lim

    def computeCLs(
        self,
        model: Data,
        marginalize: bool = False,
        toys: float = None,
        expected: Union[bool, Text] = False,
        trylasttime: bool = False,
        return_type: Text = "1-CLs",
    ) -> float:
        """
        Compute the exclusion confidence level of the model (1-CLs)
        :param model: statistical model
        :param marginalize: if true, marginalize nuisances, else profile them
        :param toys: specify number of toys. Use default is none
        :param expected: if false, compute observed,
                          true: compute a priori expected, "posteriori":
                          compute a posteriori expected
        :param trylasttime: if True, then dont try extra
        :param return_type: (Text) can be "CLs-alpha", "1-CLs", "CLs"
                        CLs-alpha: returns CLs - 0.05 (alpha)
                        1-CLs: returns 1-CLs value
                        CLs: returns CLs value
        """
        _, _, clsRoot = self.getCLsRootFunc(model, marginalize, toys, expected, trylasttime )
        ret = clsRoot(1.0, return_type=return_type)
        # its not an uppser limit on mu, its on nsig
        return ret


if __name__ == "__main__":
    C = [
        18774.2,
        -2866.97,
        -5807.3,
        -4460.52,
        -2777.25,
        -1572.97,
        -846.653,
        -442.531,
        -2866.97,
        496.273,
        900.195,
        667.591,
        403.92,
        222.614,
        116.779,
        59.5958,
        -5807.3,
        900.195,
        1799.56,
        1376.77,
        854.448,
        482.435,
        258.92,
        134.975,
        -4460.52,
        667.591,
        1376.77,
        1063.03,
        664.527,
        377.714,
        203.967,
        106.926,
        -2777.25,
        403.92,
        854.448,
        664.527,
        417.837,
        238.76,
        129.55,
        68.2075,
        -1572.97,
        222.614,
        482.435,
        377.714,
        238.76,
        137.151,
        74.7665,
        39.5247,
        -846.653,
        116.779,
        258.92,
        203.967,
        129.55,
        74.7665,
        40.9423,
        21.7285,
        -442.531,
        59.5958,
        134.975,
        106.926,
        68.2075,
        39.5247,
        21.7285,
        11.5732,
    ]
    nsignal = [x / 100.0 for x in [47, 29.4, 21.1, 14.3, 9.4, 7.1, 4.7, 4.3]]
    m = Data(
        observed=[1964, 877, 354, 182, 82, 36, 15, 11],
        backgrounds=[2006.4, 836.4, 350.0, 147.1, 62.0, 26.2, 11.1, 4.7],
        covariance=C,
        # third_moment = [ 0.1, 0.02, 0.1, 0.1, 0.003, 0.0001, 0.0002, 0.0005 ],
        third_moment=[0.0] * 8,
        nsignal=nsignal,
        name="CMS-NOTE-2017-001 model",
    )
    ulComp = UpperLimitComputer(ntoys=500, cl=0.95)
    # uls = ulComp.ulSigma ( Data ( 15,17.5,3.2,0.00454755 ) )
    # print ( "uls=", uls )
    ul_old = 131.828 * sum(
        nsignal
    )  # With respect to the older refernece value one must normalize the xsec
    print("old ul=", ul_old)
    ul = ulComp.getUpperLimitOnMu(m, marginalize=True)

    cls = ulComp.computeCLs(m, marginalize=True)
    print("ul (marginalized)", ul)
    print("CLs (marginalized)", cls)

    ul = ulComp.getUpperLimitOnMu(m, marginalize=False)
    cls = ulComp.computeCLs(m, marginalize=False)
    print("ul (profiled)", ul)
    print("CLs (profiled)", cls)

    """
    results:
    old ul= 180.999844
    ul (marginalized) 184.8081186162269
    CLs (marginalized) 1.0
    ul (profiled) 180.68039063387553
    CLs (profiled) 0.75
    """
