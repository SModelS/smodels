#!/usr/bin/env python3

"""
.. module:: analysesCombinations
   :synopsis: a module with the methods required for computing the likelihood and
              upper limits for a combination of theory predictions (analyses combination).
              Mostly used by the TheoryPredictionsCombiner class.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Jamie Yellen <j.yellen.1@research.gla.ac.uk>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

import numpy as np
from smodels.base.physicsUnits import fb
from smodels.base.smodelsLogging import logger
from smodels.statistics.basicStats import CLsfromNLL, determineBrentBracket, \
     findRoot
import scipy.optimize as optimize
from smodels.statistics.exceptions import SModelSStatisticsError as SModelSError
from typing import Text, Tuple, Callable, Union, Dict
from smodels.tools.caching import roundCache, lru_cache
from smodels.statistics.basicStats import observed, apriori, aposteriori, NllEvalType
from smodels.matching.theoryPrediction import mu_digits

class AnaCombLikelihoodComputer(object):

    def __init__(self, theoryPredictions: list, deltas_rel=None):
        """constructor.

        :param theoryPredictions: the List of theory predictions
        :param deltas_rel: relative uncertainty in signal (float). \
                           Default value is 20%.
        """
        if len(theoryPredictions) == 0:
            raise SModelSError("asking to combine zero theory predictions")
        self.theoryPredictions = theoryPredictions
        if deltas_rel is None:
            from smodels.tools.runtime import _deltas_rel_default

            deltas_rel = _deltas_rel_default
        self.deltas_rel = deltas_rel

    def nll(
        self,
        mu: float = 1.0,
        evaluationType : NllEvalType = observed,
        asimov: Union[None,float] = None,
    ) -> float:
        """
        Compute the negative log likelihood at a given mu. For now this is
        just a proxy method to ease transition of codes that depend on SModelS.

        :param mu: signal strength
        :param expected: one of: observed, apriori, aposteriori
        :param asimov: if not None, compute llhd for asimov data with mu=asimov
        """
        return self.likelihood ( mu, evaluationType, return_nll=True, asimov=asimov )

    @roundCache(argname='mu',argpos=1,digits=mu_digits)
    def likelihood(
        self,
        mu: float = 1.0,
        evaluationType : NllEvalType = observed,
        return_nll: bool = False,
        asimov: Union[None,float] = None,
    ) -> float:
        """
        Compute the likelihood at a given mu

        :param mu: signal strength
        :param expected: one of: observed, apriori, aposteriori
        :param return_nll: if True, return negative log likelihood, else likelihood
        :param asimov: if not None, compute llhd for asimov data with mu=asimov
        """
        try:
            mu = mu[0]  # some of these methods use arrays with a single element
        except:
            pass

        nll = 0.0
        changed = False
        for tp in self.theoryPredictions:
            tmp = tp.likelihood(mu, evaluationType=evaluationType, return_nll=True, asimov=asimov)
            if tmp != None:
                nll = nll + tmp     #Add neg log llhds
                changed = True
            else:
                return None
        if changed == False:
            return None
        if not return_nll:
            llhd = np.exp(-nll)
            return llhd
# print ( f"@@AC mu={mu} evaluationType {expected} asimov {asimov} nll {nll}" )

        return nll

    @lru_cache
    def lmax(
        self,
        allowNegativeSignals: bool = False,
        evaluationType : NllEvalType = observed,
        return_nll: bool = False,
        asimov: Union[None,float] = None,
    ) -> Union[Dict, None]:
        """find muhat and lmax.

        :param allowNegativeSignals: if true, then also allow for negative values
        :param expected: one of: observed, apriori, aposteriori
        :param return_nll: if true, return negative log max likelihood instead of lmax
        :returns: mu_hat, i.e. the maximum likelihood estimate of mu, if extended \
                  output is requested, it returns a dictionary with mu_hat, \
                  sigma_mu -- the standard deviation around mu_hat, and lmax, \
                  i.e. the likelihood at mu_hat
        """


        muhats, weighted = [], []
        totweight = 0.0
        for tp in self.theoryPredictions:
            muhat = asimov
            if asimov is None:
                muhat = tp.muhat(evaluationType=evaluationType)
            sigma_mu = tp.sigma_mu(evaluationType=evaluationType)
            if sigma_mu in [None, 0.0]:
                sigma_mu = 1.0  # unity weights if no weights
            if muhat != None:
                muhats.append(muhat)
                w = 1.0 / sigma_mu**2
                weighted.append(w * muhat)
                totweight += w
        # for a single theory prediction, we return just that
        if len(muhats)==1:
            if muhat < 0. and not allowNegativeSignals:
                muhat = 0.
            retllh = self.theoryPredictions[0].likelihood ( muhat, return_nll = return_nll, evaluationType = evaluationType, asimov=asimov )
            ret = {"muhat": muhat, "sigma_mu": sigma_mu, "lmax": retllh}
            return ret

        if len(muhats) == 0:
            logger.error(f"asked to compute muhat for combination, but no individual values")
            ret = {"muhat": None, "sigma_mu": None, "lmax": None}
            return  ret


        toTry = [sum(weighted) / totweight]

        #Add additional initialization points near the first one to make sure that minima is not lost due to small differences in the initialization values
        fluc = max(1e-6, 1e-5*toTry[0])  #incase toTry[0] = 0.0, take 1e-06
        toTry += [toTry[0] - fluc, toTry[0] + fluc]

        def fun(mu):
            # Make sure to always compute the correct llhd value (from theoryPrediction)
            # and not used the cached value (which is constant for mu~=1 an mu~=0)
            if isinstance(mu,np.ndarray):
                x = float(mu[0])
            else:
                x = float(mu)
            return self.likelihood(x, evaluationType=evaluationType, return_nll=True, asimov=asimov)

        if allowNegativeSignals:
            toTry += [1.0, 0.0, 3.0, -1.0, 10.0, -3.0, 0.1, -0.1]
        else:
            toTry += [1.0, 0.0, 3.0, 10.0, 0.1]
        for mu0 in toTry:
            # Minimize with a delta_mu = 1e-3*mu0 when computing the first derivative
            # (if delta_mu is too small, the derivative might give zero and terminate the
            # minimization) o  = scipy.optimize.minimize(fun, mu0 ) ## the unbounded method
            mxm = float(max(muhats))
            upper = 3.0 * mxm
            if upper <= 1.5:  #
                upper = 1.5
            bounds = [0, upper]
            if allowNegativeSignals:
                m = float(min(muhats))
                if m < 0.0:
                    bounds[0] = 3.0 * m
            bounds = [
                tuple(bounds),
            ]
            # if bounds[1] < bounds[0]:
            #    logger.error ( f"bounds are reversed, this should not happen" )
            o = optimize.minimize(fun, mu0, bounds=bounds, tol=1e-9)
            if not o.success:
                logger.debug(
                    f"combiner.lmax did not terminate successfully: {o.message} "
                    f"mu_hat={o.x} hess={o.hess_inv}"
                )
            # the inverted hessian is a good approximation for the variance at the
            # minimum
            invh = o.hess_inv
            try:
                invh = invh.todense()
            except AttributeError as e:
                pass
            hessian = invh[0][0]
            nll_ = o.fun
            #removed if condition: "and nll_ < 998.0", because sometimes nll_ > 998.0. Not necessarily a bad thing?
            #However lmax = np.exp(-nll_) for large nll -> 0
            if hessian > 0.0 and o.success:  # found a maximum. all good. #and nll_ < 998.0
                break
            # the hessian is negative meaning we found a maximum, not a minimum
            if hessian <= 0.0:
                logger.debug(
                    f"combiner.lmax the hessian {hessian} is negative at mu_hat={o.x}. "+
                    "try again with different initialisation."
                )
        mu_hat = o.x[0]
        lmax = np.exp(-o.fun)  # fun is *always* nll
        if hessian < 0.0:    #or nll_ > 998.0: again remove nll_ cut
            logger.error(
                "tried with several starting points to find maximum, always ended up in minimum. "
                "bailing out."
            )
            return None

        if not allowNegativeSignals and mu_hat < 0.0:
            mu_hat = 0.0  # fixme for this case we should reevaluate the hessian!
        sigma_mu = np.sqrt(hessian)
        retllh = lmax
        if return_nll:
            retllh = nll_
        ret = {"muhat": mu_hat, "sigma_mu": sigma_mu, "lmax": retllh}       #return nll_ too in the dictionary?
        return ret

    @lru_cache
    def getUpperLimitOnMu(self, evaluationType : NllEvalType=observed, allowNegativeSignals = False ):
        """get upper limit on signal strength multiplier, i.e. value for mu for \
           which CLs = 0.95

        :param expected: one of: observed, apriori, aposteriori
        :returns: upper limit on signal strength multiplier mu
        """
        mu_hat, sigma_mu, clsRoot = self.getCLsRootFunc(evaluationType=evaluationType,
                                                        allowNegativeSignals=allowNegativeSignals)

        a, b = determineBrentBracket(mu_hat, sigma_mu, clsRoot,
                                     allowNegative = allowNegativeSignals )
        mu_lim = findRoot(clsRoot, a, b, rtol=1e-03, xtol=1e-06 )
        return mu_lim

    def getUpperLimitOnSigmaTimesEff(self, evaluationType : NllEvalType=observed, allowNegativeSignals= False):
        """upper limit on the fiducial cross section sigma times efficiency,
            summed over all signal regions, i.e. sum_i xsec^prod_i eff_i
            obtained from the defined Data (using the signal prediction
            for each signal region/dataset), by using
            the q_mu test statistic from the CCGV paper (arXiv:1007.1727).

        :param expected: one of: observed, apriori, aposteriori
        :returns: upper limit on fiducial cross section
        """
        ul = self.getUpperLimitOnMu(evaluationType=evaluationType,
                                    allowNegativeSignals=allowNegativeSignals)

        if ul == None:
            return ul
        xsec = 0.0*fb
        for tp in self.theoryPredictions:
            xsec += tp.xsection
        return ul * xsec

    def getCLsRootFunc(self, evaluationType : NllEvalType=observed,
			    allowNegativeSignals : bool = False ) -> Tuple[float, float, Callable]:
        """
        Obtain the function "CLs-alpha[0.05]" whose root defines the upper limit,
        plus mu_hat and sigma_mu

        :param expected: one of: observed, apriori, aposteriori
        """
        fmh = self.lmax(evaluationType=evaluationType, allowNegativeSignals=allowNegativeSignals)
        mu_hat, sigma_mu, _ = fmh["muhat"], fmh["sigma_mu"], fmh["lmax"]
        mu_hat = mu_hat if mu_hat is not None else 0.0
        nll0 = self.likelihood(mu_hat, evaluationType=evaluationType, return_nll=True,asimov=None )
        # a posteriori evaluationType is needed here
        # mu_hat is mu_hat for signal_rel
        fmh = self.lmax(evaluationType=evaluationType, allowNegativeSignals=allowNegativeSignals,
                             return_nll=True, asimov = 0. )
        _, _, nll0A = fmh["muhat"], fmh["sigma_mu"], fmh["lmax"]

        # logger.error ( f"COMB nll0A {nll0A:.3f} mu_hatA {mu_hatA:.3f}" )
        # return 1.

        def clsRoot(mu: float, return_type: Text = "CLs-alpha") -> float:
            # at - infinity this should be .95,
            # at + infinity it should -.05
            # Make sure to always compute the correct llhd value (from theoryPrediction)
            # and not used the cached value (which is constant for mu~=1 an mu~=0)
            nll = self.likelihood(mu, return_nll=True, evaluationType=evaluationType, asimov = None)
            nllA = self.likelihood(mu, evaluationType=evaluationType, return_nll=True, asimov = 0. )

            return CLsfromNLL(nllA, nll0A, nll, nll0, (mu_hat>mu),
                     return_type=return_type) if nll and nllA is not None else None

        return mu_hat, sigma_mu, clsRoot

    @roundCache(argname='mu',argpos=1,digits=mu_digits)
    def CLs( self, mu : float = 1., evaluationType : NllEvalType=observed,
             return_type: Text = "CLs" ):
        """
        Compute the exclusion confidence level of the model

        :param mu: compute for the parameter of interest mu
        :param expected: one of: observed, apriori, aposteriori
        :param return_type: (Text) can be "CLs-alpha", "1-CLs", "CLs" \
                        CLs-alpha: returns CLs - 0.05 \
                        alpha-CLs: returns 0.05 - CLs \
                        1-CLs: returns 1-CLs value \
                        CLs: returns CLs value
        """
        assert return_type in ["CLs-alpha", "alpha-CLs", "1-CLs", "CLs"], f"Unknown return type: {return_type}."
        _, _, clsRoot = self.getCLsRootFunc(evaluationType=evaluationType)

        return float(clsRoot(mu, return_type=return_type))

    def getLlhds(self,muvals,evaluationType : NllEvalType=observed,
			     normalize : bool =True):
        """
        Compute the likelihoods for the individual analyses and the combined
        likelihood.
        Returns a dictionary with the analysis IDs as keys and the likelihood values as values.

        :param muvals: List with values for the signal strenth for which the likelihoods must
                       be evaluated.
        :param expected: one of: observed, apriori, aposteriori
        :param normalize: If True normalizes the likelihood by its integral over muvals.
        """

        llhds = {}
        llhds['combined'] = np.array([self.likelihood(mu,evaluationType=evaluationType) for mu in muvals])
        tpreds = self.theoryPredictions
        for t in tpreds:
            Id = t.analysisId()
            t.computeStatistics( evaluationType = evaluationType )
            l = np.array([t.likelihood(mu,evaluationType=evaluationType) for mu in muvals])
            llhds[Id]=l

        if normalize:
            # Compute delta mu
            dmuvals = np.diff(muvals)
            dmuvals = np.append(dmuvals,dmuvals[-1])
            for Id,l in llhds.items():
                # Compute norm (integral over mu)
                norm = np.sum(l*dmuvals)
                llhds[Id] = l/norm

        return llhds
