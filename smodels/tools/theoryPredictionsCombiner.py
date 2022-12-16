#!/usr/bin/env python3

"""
.. module:: theoryPredictionsCombiner
   :synopsis: a module that deals with combining signal regions from
              different analyses, offering the same API as simplified likelihoods.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Jamie Yellen <j.yellen.1@research.gla.ac.uk>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

import numpy as np
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import fb
from smodels.tools.statistics import CLsfromNLL, determineBrentBracket
import scipy
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from typing import Text, Tuple, Callable, Union, Dict


class TheoryPredictionsCombiner(object):
    """
    Facility used to combine theory predictions from different analyes.
    """

    def __init__(self, theoryPredictions: list, slhafile=None, marginalize=False, deltas_rel=None):
        """constructor.
        :param theoryPredictions: the List of theory predictions
        :param slhafile: optionally, the slhafile can be given, for debugging
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        """
        if len(theoryPredictions) == 0:
            raise SModelSError("asking to combine zero theory predictions")
        self.theoryPredictions = theoryPredictions
        self.slhafile = slhafile
        self.marginalize = marginalize
        if deltas_rel is None:
            from smodels.tools.runtime import _deltas_rel_default

            deltas_rel = _deltas_rel_default
        self.deltas_rel = deltas_rel
        self.cachedObjs = {False: {}, True: {}, "posteriori": {}}
        self.cachedLlhds = {False: {}, True: {}, "posteriori": {}}

    def __iter__(self):
        # Iterate over theory predictions
        for model in self.theoryPredictions:
            yield model

    def __repr__(self):
        return f"{self.analysisId()} : {self.totalXsection():.3f} [fb]"

    def __str__(self):
        return self.__repr__()

    @classmethod
    def selectResultsFrom(cls, theoryPredictions, anaIDs):
        """
        Select the results from theoryPrediction list which match one
        of the IDs in anaIDs. If there are multiple predictions for the
        same ID for which a likelihood is available, it gives priority
        to the ones with largest expected r-values.

        :param theoryPredictions: list of TheoryPrediction objects
        :param anaIDs: list with the analyses IDs (in string format) to be combined
        :return: a TheoryPredictionsCombiner object for the selected predictions.
                 If no theory prediction was selected, return None.
        """

        # First select the theory predictions which correspond to the analyses to be combined
        filteredTPs = [tp for tp in theoryPredictions if tp.analysisId() in anaIDs]
        filteredIDs = set([tp.analysisId() for tp in filteredTPs])
        # Now remove results with no likelihood available
        selectedTPs = [tp for tp in filteredTPs if tp.likelihood() is not None]
        selectedIDs = set([tp.analysisId() for tp in selectedTPs])
        # Warn the user concerning results with no likelihoods:
        for anaID in filteredIDs.difference(selectedIDs):
            logger.info(
                "No likelihood available for %s. This analysis will not be used in analysis combination."
                % anaID
            )
        # If no results are available, return None
        if len(selectedTPs) == 0:
            return None

        # Define a hierarchy for the results:
        priority = {"combined": 2, "efficiencyMap": 1, "upperLimit": 0}
        # Now sort by highest priority and then by highest expected r-value:
        selectedTPs = sorted(
            selectedTPs, key=lambda tp: (priority[tp.dataType()], tp.getRValue(expected=True))
        )
        # Now get a single TP for each result
        # (the highest ranking analyses come last and are kept in the dict)
        uniqueTPs = {tp.analysisId(): tp for tp in selectedTPs}
        uniqueTPs = list(uniqueTPs.values())

        combiner = cls(uniqueTPs)
        return combiner

    def singleDecorator(function):
        """
        If the combiner holds a single theoryPrediction, calls
        the theoryPrediction function.
        """

        def wrapper(self, *args, **kwargs):
            # Check if obj contains a single prediction
            if len(self.theoryPredictions) == 1:
                tp = self.theoryPredictions[0]
                # Check if the theoryPrediction object has the function
                if hasattr(tp, function.__name__):
                    newF = getattr(tp, function.__name__)
                    # Check if the attribute is callable:
                    if callable(newF):
                        # Return value from theoryPrediction.function call
                        return newF(*args, **kwargs)

            # If anything failed, call method from combiner:
            return function(self, *args, **kwargs)

        return wrapper

    @singleDecorator
    def dataId(self):
        """
        Return a string with the IDs of all the datasets used in the combination.
        """
        ids = [str(tp.dataset.getID()) for tp in self]
        ret = ",".join(ids)

        return ret

    @singleDecorator
    def analysisId(self):
        """
        Return a string with the IDs of all the experimental results used in the combination.
        """

        ret = ",".join(sorted([tp.analysisId() for tp in self]))

        return ret

    @singleDecorator
    def dataType(self, short=False):
        """
        Return the type of dataset (=combined)
        :param: short, if True, return abbreviation (comb)
        """
        if short:
            return "comb"
        else:
            return "combined"

    @singleDecorator
    def getUpperLimit(self, expected=False, trylasttime=False):
        """get upper limit on *fiducial* cross sections
            which CLs = 0.95
        :param expected: if True, compute expected likelihood, else observed
        :returns: upper limit on *fiducial* cross sections (list)
        """
        clmu = self.getUpperLimitOnMu(expected=expected)
        ret = []
        for tp in self:
            ret.append(tp.xsection.value * clmu)
        return ret

    @singleDecorator
    def getRValue(self, expected=False):
        """obtain r-value, i.e. predicted_xsec / ul
        :param expected: if True, compute expected likelihood, else observed
        :returns: r-value
        """
        if "r" in self.cachedObjs[expected]:
            return self.cachedObjs[expected]["r"]
        clmu = self.getUpperLimitOnMu(expected=expected)
        if clmu == 0.0:
            # has no meaning, has it?
            self.cachedObjs[expected]["r"] = None
            return None
        r = 1.0 / clmu
        self.cachedObjs[expected]["r"] = r
        return r

    @singleDecorator
    def lsm(self, expected=False):
        """compute the likelihood at mu = 0.
        :param expected: if True, compute expected likelihood, else observed.
                         If "posteriori", compute posterior expected.
        """
        llhd = 1.0
        for tp in self:
            llhd = llhd * tp.likelihood(0.0, expected=expected)
        return llhd

    @singleDecorator
    def lmax(self, expected=False, allowNegativeSignals=False):
        if "lmax" not in self.cachedObjs[expected]:
            self.computeStatistics(expected, allowNegativeSignals)
        if "lmax" not in self.cachedObjs[expected]:
            self.cachedObjs[expected]["lmax"][allowNegativeSignals] = None
        return self.cachedObjs[expected]["lmax"][allowNegativeSignals]

    @singleDecorator
    def muhat(self, expected=False, allowNegativeSignals=False):
        if "muhat" not in self.cachedObjs[expected]:
            self.computeStatistics(expected, allowNegativeSignals)
        if "muhat" not in self.cachedObjs[expected]:
            self.cachedObjs[expected]["muhat"][allowNegativeSignals] = None
        return self.cachedObjs[expected]["muhat"][allowNegativeSignals]

    @singleDecorator
    def chi2(self, expected=False, allowNegativeSignals=False):
        if "llhd" not in self.cachedObjs[expected] or "lmax" not in self.cachedObjs[expected]:
            logger.error("call computeStatistics before calling chi2")
            return
        llhd = self.cachedObjs[expected]["llhd"]
        lmax = self.cachedObjs[expected]["lmax"][allowNegativeSignals]

        if llhd == 0.0:
            return 2000.0  # we cut off at > 1e-300 or so ;)
        return -2 * np.log(llhd / lmax)

    def describe(self):
        """returns a string containing a list of all analysisId and dataIds"""
        ids = [f"{tp.analysisId()}:{tp.dataId()}" for tp in self]
        return f"SRs: {', '.join(ids)}"

    @singleDecorator
    def likelihood(
        self,
        mu: float = 1.0,
        expected: Union[bool, Text] = False,
        nll: bool = False,
        useCached: bool = False,
    ) -> Union[float, None]:
        """
        Compute the likelihood at a given mu
        :param mu: signal strength
        :param expected: if True, compute expected likelihood, else observed
        :param nll: if True, return negative log likelihood, else likelihood
        :param useCached: if True, will return the cached value from theoryPrediction (if available)
        """
        if isinstance(mu, (list, np.ndarray)):
            mu = mu[0]
        # For backwards compatibility
        return_nll = nll

        if useCached and mu in self.cachedLlhds[expected]:
            llhd = self.cachedLlhds[expected][mu]
            if return_nll:
                return 999.0 if llhd == 0.0 else -np.log(llhd)
            return llhd

        # @JACK: Instead of multiplying likelihoods, sum nlls.
        # This allows for more numeric stability
        nll, changed = 0.0, False
        for tp in self:
            # Set the theory marginalize attribute to the combiner value:
            tp_marginalize = tp.marginalize
            tp.marginalize = self.marginalize
            tmp = tp.likelihood(mu, expected=expected, useCached=useCached, nll=True)
            if tmp is not None:
                nll += tmp
                changed = True
            # Restore marginalize setting:
            tp.marginalize = tp_marginalize

        self.cachedLlhds[expected][mu] = np.exp(-nll) if changed else None
        if not changed:
            return None

        return nll if return_nll else np.exp(-nll)

    @singleDecorator
    def computeStatistics(self, expected=False, allowNegativeSignals=False):
        """
        Compute the likelihoods, chi2, and upper limit for this combination
        :param expected: if True, compute expected likelihood, else observed
        """
        if self.theoryPredictions is None:
            return
        if expected == "both":
            for e in [False, True]:
                self.computeStatistics(e, allowNegativeSignals)
            return

        # Negative log likelihood at mu=1 and nll at mu=0
        nll1, nll0 = 0.0, 0.0
        for tp in self:
            tp_marginalize = tp.marginalize
            tp.marginalize = self.marginalize
            tmp1 = tp.likelihood(1.0, expected=expected, useCached=True, nll=True)
            nll1 += tmp1 if tmp1 is not None else 0.0

            tmp0 = tp.likelihood(0.0, expected=expected, useCached=True, nll=True)
            nll0 += tmp0 if tmp0 is not None else 0.0
            # Restore marginalize setting:
            tp.marginalize = tp_marginalize

        self.cachedObjs[expected]["llhd"] = nll1
        self.cachedObjs[expected]["lsm"] = nll0
        if expected:
            if not "lmax" in self.cachedObjs[expected]:
                self.cachedObjs[expected]["lmax"] = {}
            self.cachedObjs[expected]["lmax"][allowNegativeSignals] = nll0
        else:
            fmh = self.findMuHat(expected=False, extended_output=True)
            self.mu_hat, self.sigma_mu, lmax = fmh["muhat"], fmh["sigma_mu"], fmh["lmax"]
            if not "lmax" in self.cachedObjs[expected]:
                self.cachedObjs[expected]["lmax"] = {}
            self.cachedObjs[expected]["lmax"][allowNegativeSignals] = lmax
            if not "muhat" in self.cachedObjs[expected]:
                self.cachedObjs[expected]["muhat"] = {}
            self.cachedObjs[expected]["muhat"][allowNegativeSignals] = self.mu_hat
            if not "sigma_mu" in self.cachedObjs[expected]:
                self.cachedObjs[expected]["sigma_mu"] = {}
            self.cachedObjs[expected]["sigma_mu"][allowNegativeSignals] = self.mu_hat

    @singleDecorator
    def totalXsection(self):
        ret = 0.0 * fb
        if self.theoryPredictions is not None:
            for tp in self:
                ret += tp.xsection.value
        return ret

    @singleDecorator
    def getmaxCondition(self):
        """
        Returns the maximum xsection from the list conditions

        :returns: maximum condition xsection (float)
        """
        conditions = [tp.getmaxCondition() for tp in self]
        return max(conditions)

    def findMuHat(
        self,
        allowNegativeSignals: bool = False,
        expected: Union[bool, Text] = False,
        extended_output: bool = False,
        nll: bool = False,
    ) -> Union[Dict, float, None]:
        """
        Finds the POI (signal strength; mu) that maximizes the likelihood.

        :param allowNegativeSignals: (bool) if true, then also allow for negative values
        :param expected: (Union[bool, Text]) if true, compute expected prior (=lsm),
                         if "posteriori" compute posteriori expected.
        :param extended_output: if true, return also sigma_mu, the estimate of the
                         error of mu_hat, and lmax, the likelihood at mu_hat
        :param nll: (bool) if true, return negative log max likelihood instead of lmax
        :return:  mu_hat, i.e. the maximum likelihood estimate of mu, if extended
                  output is requested, it returns a dictionary with mu_hat,
                  sigma_mu -- the standard deviation around mu_hat, and lmax,
                  i.e. the likelihood at mu_hat
        """

        muhats, weighted_mu, weights, sigma_mus = [], [], [], []
        for tp in self:
            muhat = tp.muhat(expected=expected, allowNegativeSignals=allowNegativeSignals)
            sigma_mu = tp.sigma_mu(expected=expected, allowNegativeSignals=allowNegativeSignals)
            # Set sigma_mu = 1 and mu = 0. in case computation fails
            sigma_mu = 1.0 if sigma_mu in [None, 0.0] else sigma_mu
            muhat = 0.0 if muhat is None else muhat

            # Do not allow mu=0 at this stage breaks the computation
            if muhat > 0.0:
                muhats.append(muhat)
                sigma_mus.append(sigma_mu)
                weighted_mu.append(muhat / sigma_mu**2)
                weights.append(1.0 / sigma_mu**2)

        nll_ = None
        if len(muhats) > 1:
            combined_muhat = (
                sum(weighted_mu) / sum(weights)
                if allowNegativeSignals
                else max(sum(weighted_mu) / sum(weights), 0.0)
            )

            # @JACK: for multi-POI analyses this function will not work,
            # likelihood function needs to be adapted for such cases
            negloglikelihood = lambda mu: self.likelihood(
                mu[0], expected=expected, nll=True, useCached=False
            )

            # @JACK: It is possible to allow user to modify the optimiser properties in the future
            opt = scipy.optimize.minimize(
                negloglikelihood,
                combined_muhat,
                method="COBYLA",
                tol=1e-6,
                options={"maxiter": 10000},
            )
            if not opt.success:
                logger.debug(
                    f"combiner.findMuHat did not terminate successfully: {opt.message} \n"
                    f"mu_hat={opt.x[0]:.5e} slhafile={self.slhafile}"
                )

            nll_ = opt.fun if nll else np.exp(-opt.fun)
            combined_muhat = opt.x[0] if allowNegativeSignals else max(opt.x[0], 0.0)
            sigma_mu = combined_muhat * np.sqrt(
                sum((s / m) ** 2 for s, m in zip(sigma_mus, muhats))
            )  # add in quadrature
        elif len(muhats) == 1:
            combined_muhat = muhats[-1] if allowNegativeSignals else max(muhats[-1], 0.0)
            sigma_mu = sigma_mus[-1]
        else:
            combined_muhat, sigma_mu = 0.0, 0.0

        if extended_output:
            if nll_ is None:
                nll_ = self.likelihood(combined_muhat, expected=expected, nll=nll)
            return {"muhat": combined_muhat, "sigma_mu": sigma_mu, "lmax": nll_}

        return combined_muhat

    def getCLsRootFunc(self, expected: bool = False) -> Tuple[float, float, Callable]:
        """
        Obtain the function "CLs-alpha[0.05]" whose root defines the upper limit,
        plus mu_hat and sigma_mu
        :param expected: if True, compute expected likelihood, else observed
        """
        fmh = self.findMuHat(expected=expected, allowNegativeSignals=False, extended_output=True)
        mu_hat, sigma_mu, lmax = fmh["muhat"], fmh["sigma_mu"], fmh["lmax"]
        mu_hat = mu_hat if mu_hat is not None else 0.0
        nll0 = self.likelihood(mu_hat, expected=expected, nll=True)
        # a posteriori expected is needed here
        # mu_hat is mu_hat for signal_rel
        fmh = self.findMuHat(
            expected="posteriori", allowNegativeSignals=False, nll=True, extended_output=True
        )
        mu_hatA, _, nll0A = fmh["muhat"], fmh["sigma_mu"], fmh["lmax"]

        def clsRoot(mu: float, return_type: Text = "CLs-alpha") -> float:
            # at - infinity this should be .95,
            # at + infinity it should -.05
            # Make sure to always compute the correct llhd value (from theoryPrediction)
            # and not used the cached value (which is constant for mu~=1 an mu~=0)
            nll = self.likelihood(mu, expected=expected, useCached=False, nll=True)
            nllA = self.likelihood(mu, expected="posteriori", useCached=False, nll=True)
            return CLsfromNLL(nllA, nll0A, nll, nll0, return_type=return_type)

        return mu_hat, sigma_mu, clsRoot

    @singleDecorator
    def getUpperLimitOnMu(self, expected=False):
        """get upper limit on signal strength multiplier, i.e. value for mu for
            which CLs = 0.95
        :param expected: if True, compute expected likelihood, else observed
        :returns: upper limit on signal strength multiplier mu
        """
        mu_hat, sigma_mu, clsRoot = self.getCLsRootFunc(expected=expected)

        a, b = determineBrentBracket(mu_hat, sigma_mu, clsRoot, allowNegative=False)
        mu_lim = scipy.optimize.brentq(clsRoot, a, b, rtol=1e-03, xtol=1e-06)
        self.cachedObjs[expected]["UL"] = mu_lim
        return mu_lim

    def computeCLs(self, expected: bool = False, return_type: Text = "1-CLs"):
        """
        Compute the exclusion confidence level of the model (1-CLs)
        :param expected: if false, compute observed, true: compute a priori expected
        :param return_type: (Text) can be "CLs-alpha", "1-CLs", "CLs"
                        CLs-alpha: returns CLs - 0.05
                        1-CLs: returns 1-CLs value
                        CLs: returns CLs value
        """
        assert return_type in ["CLs-alpha", "1-CLs", "CLs"], f"Unknown return type: {return_type}."
        _, _, clsRoot = self.getCLsRootFunc(expected=expected)

        return clsRoot(1.0, return_type=return_type)

    def getLlhds(self, muvals, expected=False, normalize=True):
        """
        Compute the likelihoods for the individual analyses and the combined
        likelihood.
        Returns a dictionary with the analysis IDs as keys and the likelihood values as values.

        :param muvals: List with values for the signal strenth for which the likelihoods must
                       be evaluated.
        :param expected: If True returns the expected likelihood values.
        :param normalize: If True normalizes the likelihood by its integral over muvals.
        """

        llhds = {}
        llhds["combined"] = np.array([self.likelihood(mu, expected=expected) for mu in muvals])
        tpreds = self.theoryPredictions
        for t in tpreds:
            Id = t.analysisId()
            t.computeStatistics(expected=expected)
            l = np.array([t.likelihood(mu, expected=expected) for mu in muvals])
            llhds[Id] = l

        if normalize:
            # Compute delta mu
            dmuvals = np.diff(muvals)
            dmuvals = np.append(dmuvals, dmuvals[-1])
            for Id, l in llhds.items():
                # Compute norm (integral over mu)
                norm = np.sum(l * dmuvals)
                llhds[Id] = l / norm

        return llhds
