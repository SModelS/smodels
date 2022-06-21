"""
.. module:: theoryPrediction
   :synopsis: Provides a class to encapsulate the results of the computation of
              reference cross sections and related functions.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. autofunction:: _getElementsFrom
"""

from smodels.theory import clusterTools, crossSection, element
from smodels.theory.auxiliaryFunctions import cSim, cGtr, elementsInStr, average
from smodels.tools.physicsUnits import TeV, fb
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.experiment.datasetObj import CombinedDataSet
from smodels.tools.smodelsLogging import logger
from smodels.tools.statistics import TruncatedGaussians
from smodels.tools.srCombinations import getCombinedStatistics, \
            getCombinedUpperLimitFor, getCombinedLikelihood
import itertools
import numpy as np


class TheoryPrediction(object):
    """
    An instance of this class represents the results of the theory prediction
    for an analysis.
    """

    def __init__(self, marginalize=False, deltas_rel=None):
        """a theory prediction. marginalize and deltas_rel are meant to be
            constants
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        """
        self.analysis = None
        self.xsection = None
        self.conditions = None
        self.mass = None
        self.totalwidth = None
        self.marginalize = marginalize
        if deltas_rel is None:
            from smodels.tools.runtime import _deltas_rel_default

            deltas_rel = _deltas_rel_default
        self.deltas_rel = deltas_rel
        self.cachedObjs = {False: {}, True: {}, "posteriori": {}}
        self.cachedLlhds = {False: {}, True: {}, "posteriori": {}}

    def __str__(self):
        ret = "%s:%s" % (self.analysisId(), self.totalXsection())
        return ret

    def dataId(self):
        """
        Return ID of dataset
        """

        return self.dataset.getID()

    def analysisId(self):
        """
        Return experimental analysis ID
        """

        return self.dataset.globalInfo.id

    def dataType(self, short=False):
        """
        Return the type of dataset
        :param: short, if True, return abbreviation (ul,em,comb)
        """
        if short:
            t = self.dataset.getType()
            D = {"upperLimit": "ul", "efficiencyMap": "em", "combined": "comb"}
            if t in D.keys():
                return D[t]
            return "??"

        return self.dataset.getType()

    def getUpperLimit(self, expected=False):
        """
        Get the upper limit on sigma*eff.
        For UL-type results, use the UL map. For EM-Type returns
        the corresponding dataset (signal region) upper limit.
        For combined results, returns the upper limit on the
        total sigma*eff (for all signal regions/datasets).

        :param expected: return expected Upper Limit, instead of observed.
        :return: upper limit (Unum object)
        """

        # First check if the upper-limit and expected upper-limit have already been computed.
        # If not, compute it and store them.
        if "UL" not in self.cachedObjs[expected]:
            if self.dataType() == "efficiencyMap":
                self.cachedObjs[expected]["UL"] = self.dataset.getSRUpperLimit(
                    expected=expected, deltas_rel=self.deltas_rel
                )
            if self.dataType() == "upperLimit":
                self.cachedObjs[expected]["UL"] = self.dataset.getUpperLimitFor(
                    element=self.avgElement, txnames=self.txnames, expected=expected
                )
            if self.dataType() == "combined":
                # Create a list of signal events in each dataset/SR sorted according to datasetOrder
                # lumi = self.dataset.getLumi()
                if hasattr(self.dataset.globalInfo, "covariance"):
                    srNsigDict = dict(
                        [
                            [
                                pred.dataset.getID(),
                                (pred.xsection.value * pred.dataset.getLumi()).asNumber(),
                            ]
                            for pred in self.datasetPredictions
                        ]
                    )
                    srNsigs = [
                        srNsigDict[dataID] if dataID in srNsigDict else 0.0
                        for dataID in self.dataset.globalInfo.datasetOrder
                    ]
                elif hasattr(self.dataset.globalInfo, "jsonFiles"):
                    srNsigDict = dict(
                        [
                            [
                                pred.dataset.getID(),
                                (pred.xsection.value * pred.dataset.getLumi()).asNumber(),
                            ]
                            for pred in self.datasetPredictions
                        ]
                    )
                    srNsigs = [
                        srNsigDict[ds.getID()] if ds.getID() in srNsigDict else 0.0
                        for ds in self.dataset._datasets
                    ]
                self.cachedObjs[expected]["UL"] = getCombinedUpperLimitFor(
                    self.dataset, srNsigs, expected=expected, deltas_rel=self.deltas_rel
                )

        # Return the expected or observed UL:
        # if not self.cachedObjs[expected]["UL"]:
        #    self.cachedObjs[expected]["UL"]=None
        return self.cachedObjs[expected]["UL"]

    def getRValue(self, expected=False):
        """
        Get the r value = theory prediction / experimental upper limit
        """
        if "r" not in self.cachedObjs[expected]:
            upperLimit = self.getUpperLimit(expected)
            if upperLimit is None or upperLimit.asNumber(fb) == 0.0:
                r = None
                self.cachedObjs[expected]["r"] = r
                return r
            else:
                r = (self.xsection.value / upperLimit).asNumber()
                self.cachedObjs[expected]["r"] = r
                return r
        return self.cachedObjs[expected]["r"]

    def lsm(self, expected=False):
        """likelihood at SM point, same as .def likelihood( ( mu = 0. )"""
        if "lsm" not in self.cachedObjs[expected]:
            self.computeStatistics(expected)
        if "lsm" not in self.cachedObjs[expected]:
            self.cachedObjs[expected]["lsm"] = None
        return self.cachedObjs[expected]["lsm"]

    def lmax(self, expected=False, allowNegativeSignals=False):
        """likelihood at mu_hat"""
        if not "lmax" in self.cachedObjs[expected]:
            self.computeStatistics(expected, allowNegativeSignals)
        if (
            "lmax" in self.cachedObjs[expected]
            and not allowNegativeSignals in self.cachedObjs[expected]["lmax"]
        ):
            self.computeStatistics(expected, allowNegativeSignals)
        if not "lmax" in self.cachedObjs[expected]:
            self.cachedObjs[expected]["lmax"] = {allowNegativeSignals: None}
        return self.cachedObjs[expected]["lmax"][allowNegativeSignals]

    def sigma_mu(self, expected=False, allowNegativeSignals=False):
        """sigma_mu of mu_hat"""
        if not "sigma_mu" in self.cachedObjs[expected]:
            self.computeStatistics(expected, allowNegativeSignals)
        if not "sigma_mu" in self.cachedObjs[expected]:
            return None
        if not allowNegativeSignals in self.cachedObjs[expected]["sigma_mu"]:
            self.computeStatistics(expected, allowNegativeSignals)
        return self.cachedObjs[expected]["sigma_mu"][allowNegativeSignals]

    def muhat(self, expected=False, allowNegativeSignals=False):
        """position of maximum likelihood"""
        if not "muhat" in self.cachedObjs[expected]:
            self.computeStatistics(expected, allowNegativeSignals)
        if not allowNegativeSignals in self.cachedObjs[expected]["muhat"]:
            self.computeStatistics(expected, allowNegativeSignals)
        if (
            "muhat" in self.cachedObjs[expected]
            and not allowNegativeSignals in self.cachedObjs[expected]["muhat"]
        ):
            self.computeStatistics(expected, allowNegativeSignals)
        if not "muhat" in self.cachedObjs[expected]:
            self.cachedObjs[expected]["muhat"] = {allowNegativeSignals: None}
        ret = self.cachedObjs[expected]["muhat"][allowNegativeSignals]
        return ret

    def chi2(self, expected=False):
        if not "chi2" in self.cachedObjs[expected]:
            self.computeStatistics(expected)
        if not "chi2" in self.cachedObjs[expected]:
            self.cachedObjs[expected]["chi2"] = None
        return self.cachedObjs[expected]["chi2"]

    def likelihood(self, mu=1.0, expected=False, nll=False, useCached=True):
        """
        get the likelihood for a signal strength modifier mu
        :param expected: compute expected, not observed likelihood. if "posteriori",
                         compute expected posteriori.
        :param nll: if True, return negative log likelihood, else likelihood
        :param useCached: if True, will return the cached value, if available
        """
        if useCached and mu in self.cachedLlhds[expected]:
            llhd = self.cachedLlhds[expected][mu]
            if nll:
                if llhd == 0.0:
                    return 700.0
                return -np.log(llhd)
            return llhd
        # self.computeStatistics(expected=expected)
        if useCached:
            if "llhd" in self.cachedObjs[expected] and abs(mu - 1.0) < 1e-5:
                llhd = self.cachedObjs[expected]["llhd"]
                if nll:
                    return -np.log(llhd)
                return llhd
            if "lsm" in self.cachedObjs[expected] and abs(mu) < 1e-5:
                lsm = self.cachedObjs[expected]["lsm"]
                if nll:
                    return -np.log(lsm)
                return lsm
        lumi = self.dataset.getLumi()
        if self.dataType() == "combined":
            srNsigDict = dict(
                [
                    [pred.dataset.getID(), (pred.xsection.value * lumi).asNumber()]
                    for pred in self.datasetPredictions
                ]
            )
            srNsigs = [
                srNsigDict[ds.getID()] if ds.getID() in srNsigDict else 0.0
                for ds in self.dataset._datasets
            ]
            llhd = getCombinedLikelihood(
                self.dataset, srNsigs, self.marginalize, self.deltas_rel, expected=expected, mu=mu
            )
        if self.dataType() == "efficiencyMap":
            nsig = (mu * self.xsection.value * lumi).asNumber()
            llhd = self.dataset.likelihood(
                nsig, marginalize=self.marginalize, deltas_rel=self.deltas_rel, expected=expected
            )

        if self.dataType() == "upperLimit":
            # these fits only work with negative signals!
            llhd, chi2 = self.likelihoodFromLimits(
                mu, chi2also=True, expected=expected, allowNegativeSignals=True
            )
        self.cachedLlhds[expected][mu] = llhd
        if nll:
            if llhd == 0.0:
                return 700.0
            return -np.log(llhd)
        return llhd

    def likelihoodFromLimits(
        self, mu=1.0, expected=False, chi2also=False, corr=0.6, allowNegativeSignals=False
    ):
        """compute the likelihood from expected and observed upper limits.
        :param expected: compute expected, not observed likelihood
        :param mu: signal strength multiplier, applied to theory prediction. If None,
                   then find muhat
        :param chi2also: if true, return also chi2
        :param corr: correction factor:
                 ULexp_mod = ULexp / (1. - corr*((ULobs-ULexp)/(ULobs+ULexp)))
                 a factor of corr = 0.6 is proposed.
        :param allowNegativeSignals: if False, then negative nsigs are replaced with 0.
        :returns: likelihood; none if no expected upper limit is defined.
        """
        # marked as experimental feature
        from smodels.tools.runtime import experimentalFeatures

        if not experimentalFeatures():
            if chi2also:
                return (None, None)
            return None
        if not hasattr(self, "avgElement"):
            logger.error("theory prediction %s has no average element! why??" % self.analysisId())
            if chi2also:
                return (None, None)
            return None

        eul = self.dataset.getUpperLimitFor(
            element=self.avgElement, txnames=self.txnames, expected=True
        )
        if eul is None:
            if chi2also:
                return (None, None)
            return None
        ul = self.dataset.getUpperLimitFor(
            element=self.avgElement, txnames=self.txnames, expected=False
        )
        lumi = self.dataset.getLumi()
        computer = TruncatedGaussians ( ul, eul, self.xsection.value, lumi=lumi, corr = corr )
        ret = computer.likelihood ( mu )
        llhd, muhat, sigma_mu = ret["llhd"], ret["muhat"], ret["sigma_mu"]
        """
        if False: # muhat < 0.0 and allowNegativeSignals == False:
            oldmuhat = muhat
            oldl = llhd
            muhat = 0.0
            llhd, muhat, sigma_mu = likelihoodFromLimits(
                ulN, eulN, nsig, 0.0, allowNegativeMuhat=True, corr=corr
            )
            nllhd, nmuhat, nsigma_mu = computer.likelihood ( 0. )
        """
        self.muhat_ = muhat
        self.sigma_mu_ = sigma_mu 
        if chi2also:
            return (llhd, computer.chi2 ( ) )
        return llhd

    def computeStatistics(self, expected=False, allowNegativeSignals=False):
        """
        Compute the likelihoods, chi2 and upper limit for this theory prediction.
        The resulting values are stored as the likelihood, lmax, lsm and chi2
        attributes (chi2 being phased out).
        :param expected: computed expected quantities, not observed
        """
        if not "lmax" in self.cachedObjs[expected]:
            self.cachedObjs[expected]["lmax"] = {}
            self.cachedObjs[expected]["muhat"] = {}

        if self.dataType() == "upperLimit":
            llhd, chi2 = self.likelihoodFromLimits(1.0, expected=expected, chi2also=True)
            lsm = self.likelihoodFromLimits(0.0, expected=expected, chi2also=False)
            lmax = self.likelihoodFromLimits(
                None, expected=expected, chi2also=False, allowNegativeSignals=True
            )
            if allowNegativeSignals == False and hasattr(self, "muhat_") and self.muhat_ < 0.0:
                self.muhat_ = 0.0
                lmax = lsm
            self.cachedObjs[expected]["llhd"] = llhd
            self.cachedObjs[expected]["lsm"] = lsm
            self.cachedObjs[expected]["lmax"][allowNegativeSignals] = lmax
            self.cachedObjs[expected]["chi2"] = chi2
            if hasattr(self, "muhat_"):
                self.cachedObjs[expected]["muhat"][allowNegativeSignals] = self.muhat_
            if hasattr(self, "sigma_mu_"):
                if not "sigma_mu" in self.cachedObjs[expected]:
                    self.cachedObjs[expected]["sigma_mu"] = {}
                self.cachedObjs[expected]["sigma_mu"][allowNegativeSignals] = self.sigma_mu_

        elif self.dataType() == "efficiencyMap":
            lumi = self.dataset.getLumi()
            nsig = (self.xsection.value * lumi).asNumber()
            llhd = self.dataset.likelihood(
                nsig, marginalize=self.marginalize, deltas_rel=self.deltas_rel, expected=expected
            )
            llhd_sm = self.dataset.likelihood(
                nsig=0.0,
                marginalize=self.marginalize,
                deltas_rel=self.deltas_rel,
                expected=expected,
            )
            llhd_max = self.dataset.lmax(
                marginalize=self.marginalize,
                deltas_rel=self.deltas_rel,
                allowNegativeSignals=allowNegativeSignals,
                expected=expected,
            )
            muhat = None
            if hasattr(self.dataset, "muhat"):
                muhat = self.dataset.muhat / nsig
            if hasattr(self.dataset, "sigma_mu"):
                sigma_mu = float(self.dataset.sigma_mu / nsig)
                if not "sigma_mu" in self.cachedObjs[expected]:
                    self.cachedObjs[expected]["sigma_mu"] = {}
                self.cachedObjs[expected]["sigma_mu"][allowNegativeSignals] = sigma_mu
            self.cachedObjs[expected]["llhd"] = llhd
            self.cachedObjs[expected]["lsm"] = llhd_sm
            self.cachedObjs[expected]["lmax"][allowNegativeSignals] = llhd_max
            self.cachedObjs[expected]["muhat"][allowNegativeSignals] = muhat
            from smodels.tools.statistics import chi2FromLmax

            self.cachedObjs[expected]["chi2"] = chi2FromLmax(llhd, llhd_max)

        elif self.dataType() == "combined":
            lumi = self.dataset.getLumi()
            # Create a list of signal events in each dataset/SR sorted according to datasetOrder
            srNsigDict = dict(
                [
                    [pred.dataset.getID(), (pred.xsection.value * lumi).asNumber()]
                    for pred in self.datasetPredictions
                ]
            )
            srNsigs = [
                srNsigDict[ds.getID()] if ds.getID() in srNsigDict else 0.0
                for ds in self.dataset._datasets
            ]
            # srNsigs = [srNsigDict[dataID] if dataID in srNsigDict else 0. for dataID in self.dataset.globalInfo.datasetOrder]
            s = getCombinedStatistics(
                self.dataset,
                srNsigs,
                self.marginalize,
                self.deltas_rel,
                expected=expected,
                allowNegativeSignals=allowNegativeSignals,
            )
            llhd, lmax, lsm, muhat, sigma_mu = (
                s["lbsm"],
                s["lmax"],
                s["lsm"],
                s["muhat"],
                s["sigma_mu"],
            )
            self.cachedObjs[expected]["llhd"] = llhd
            self.cachedObjs[expected]["lsm"] = lsm
            self.cachedObjs[expected]["lmax"][allowNegativeSignals] = lmax
            self.cachedObjs[expected]["muhat"][allowNegativeSignals] = muhat
            if not "sigma_mu" in self.cachedObjs[expected]:
                self.cachedObjs[expected]["sigma_mu"] = {}
            self.cachedObjs[expected]["sigma_mu"][allowNegativeSignals] = sigma_mu
            from smodels.tools.statistics import chi2FromLmax

            self.cachedObjs[expected]["chi2"] = chi2FromLmax(llhd, lmax)

    def totalXsection(self):
        return self.xsection.value

    def getmaxCondition(self):
        """
        Returns the maximum xsection from the list conditions

        :returns: maximum condition xsection (float)
        """

        if not self.conditions:
            return 0.0
        # maxcond = 0.
        values = [0.0]
        for value in self.conditions.values():
            if value == "N/A":
                return value
            if value == None:
                continue
            # print ( "value=",value,type(value),float(value) )
            # maxcond = max(maxcond,float(value))
            values.append(float(value))
        return max(values)
        # return maxcond


class TheoryPredictionList(object):
    """
    An instance of this class represents a collection of theory prediction
    objects.
    """

    def __init__(self, theoryPredictions=None, maxCond=None):
        """
        Initializes the list.

        :parameter theoryPredictions: list of TheoryPrediction objects
        :parameter maxCond: maximum relative violation of conditions for valid
        results. If defined, it will keep only the theory predictions with
        condition violation < maxCond.
        """
        self._theoryPredictions = []
        if theoryPredictions and isinstance(theoryPredictions, list):

            if not maxCond:
                self._theoryPredictions = theoryPredictions
            else:
                newPredictions = []
                for theoPred in theoryPredictions:
                    mCond = theoPred.getmaxCondition()
                    if mCond == "N/A" or round(mCond / maxCond, 2) > 1.0:
                        continue
                    else:
                        newPredictions.append(theoPred)
                self._theoryPredictions = newPredictions

    def append(self, theoryPred):
        self._theoryPredictions.append(theoryPred)

    def __str__(self):
        if len(self._theoryPredictions) == 0:
            return "no predictions."
        ret = "%d predictions: " % len(self._theoryPredictions)
        ret += ", ".join([str(s) for s in self._theoryPredictions])
        return ret

    def __iter__(self):
        for theoryPrediction in self._theoryPredictions:
            yield theoryPrediction

    def __getitem__(self, index):
        return self._theoryPredictions[index]

    def __len__(self):
        return len(self._theoryPredictions)

    def __add__(self, theoPredList):
        if isinstance(theoPredList, TheoryPredictionList):
            res = TheoryPredictionList()
            res._theoryPredictions = self._theoryPredictions + theoPredList._theoryPredictions
            return res
        else:
            return None

    def __radd__(self, theoPredList):
        if theoPredList == 0:
            return self
        else:
            return self.__add__(theoPredList)

    def sortTheoryPredictions(self):
        """
        Reverse sort theoryPredictions by R value.
        Used for printer.
        """
        self._theoryPredictions = sorted(
            self._theoryPredictions,
            key=lambda theoPred: (theoPred.getRValue() is not None, theoPred.getRValue()),
            reverse=True,
        )


def theoryPredictionsFor(
    expResult,
    smsTopList,
    maxMassDist=0.2,
    useBestDataset=True,
    combinedResults=True,
    marginalize=False,
    deltas_rel=None,
):
    """
    Compute theory predictions for the given experimental result, using the list of
    elements in smsTopList.
    For each Txname appearing in expResult, it collects the elements and
    efficiencies, combine the masses (if needed) and compute the conditions
    (if exist).

    :parameter expResult: expResult to be considered (ExpResult object), if list of ExpResults is given, produce theory predictions for all
    :parameter smsTopList: list of topologies containing elements
                           (TopologyList object)
    :parameter maxMassDist: maximum mass distance for clustering elements (float)
    :parameter useBestDataset: If True, uses only the best dataset (signal region).
               If False, returns predictions for all datasets (if combinedResults is False),
               or only the combinedResults (if combinedResults is True).
    :parameter combinedResults: add theory predictions that result from
               combining datasets.
    :parameter marginalize: If true, marginalize nuisances. If false, profile them.
    :parameter deltas_rel: relative uncertainty in signal (float). Default value is 20%.

    :returns:  a TheoryPredictionList object containing a list of TheoryPrediction
               objects
    """
    if deltas_rel is None:
        from smodels.tools.runtime import _deltas_rel_default

        deltas_rel = _deltas_rel_default

    if type(expResult) in [list, tuple]:
        ret = []
        for er in expResult:
            tpreds = theoryPredictionsFor(
                er,
                smsTopList,
                maxMassDist,
                useBestDataset,
                combinedResults,
                marginalize,
                deltas_rel,
            )
            if tpreds:
                for tp in tpreds:
                    ret.append(tp)
        return TheoryPredictionList(ret)

    dataSetResults = []
    # Compute predictions for each data set (for UL analyses there is one single set)
    for dataset in expResult.datasets:
        predList = _getDataSetPredictions(dataset, smsTopList, maxMassDist)
        if predList:
            dataSetResults.append(predList)
    if not dataSetResults:  # no results at all?
        return None
    elif len(dataSetResults) == 1:  # only a single dataset? Easy case.
        result = dataSetResults[0]
        for theoPred in result:
            theoPred.expResult = expResult
            theoPred.deltas_rel = deltas_rel
            theoPred.upperLimit = theoPred.getUpperLimit()
        return result

    # For results with more than one dataset, return all dataset predictions
    # if useBestDataSet=False and combinedResults=False:
    if not useBestDataset and not combinedResults:
        allResults = sum(dataSetResults)
        for theoPred in allResults:
            theoPred.expResult = expResult
            theoPred.deltas_rel = deltas_rel
            theoPred.upperLimit = theoPred.getUpperLimit()
        return allResults
    elif combinedResults:  # Try to combine results
        bestResults = TheoryPredictionList()
        combinedRes = _getCombinedResultFor(dataSetResults, expResult, marginalize)
        if combinedRes is None:  # Can not combine, use best dataset:
            combinedRes = _getBestResult(dataSetResults)
        bestResults.append(combinedRes)
    else:  # Use best dataset:
        bestResults = TheoryPredictionList()
        bestResults.append(_getBestResult(dataSetResults))

    for theoPred in bestResults:
        theoPred.expResult = expResult
        theoPred.deltas_rel = deltas_rel
        theoPred.upperLimit = theoPred.getUpperLimit()

    return bestResults


def _getCombinedResultFor(dataSetResults, expResult, marginalize=False):
    """
    Compute the combined result for all datasets, if covariance
    matrices are available. Return a TheoryPrediction object
    with the signal cross-section summed over all the signal regions
    and the respective upper limit.

    :param datasetPredictions: List of TheoryPrediction objects for each signal region
    :param expResult: ExpResult object corresponding to the experimental result
    :parameter marginalize: If true, marginalize nuisances. If false, profile them.

    :return: TheoryPrediction object
    """

    if len(dataSetResults) == 1:
        return dataSetResults[0]
    elif not expResult.hasCovarianceMatrix() and not expResult.hasJsonFile():
        return None

    txnameList = []
    elementList = []
    totalXsec = None
    massList = []
    widthList = []
    PIDList = []
    datasetPredictions = []
    weights = []
    for predList in dataSetResults:
        if len(predList) != 1:
            raise SModelSError(
                "Results with multiple datasets should have a single theory prediction (EM-type)!"
            )
        pred = predList[0]
        datasetPredictions.append(pred)
        txnameList += pred.txnames
        elementList += pred.elements
        if not totalXsec:
            totalXsec = pred.xsection
        else:
            totalXsec += pred.xsection
        massList.append(pred.mass)
        widthList.append(pred.totalwidth)
        weights.append(pred.xsection.value.asNumber(fb))
        PIDList += pred.PIDs

    txnameList = list(set(txnameList))
    if None in massList:
        mass = None
        totalwidth = None
    else:
        mass = average(massList, weights=weights)
        totalwidth = average(widthList, weights=weights)

    # Create a combinedDataSet object:
    combinedDataset = CombinedDataSet(expResult)
    combinedDataset._marginalize = marginalize
    # Create a theory preidction object for the combined datasets:
    theoryPrediction = TheoryPrediction(marginalize)
    theoryPrediction.dataset = combinedDataset
    theoryPrediction.txnames = txnameList
    theoryPrediction.xsection = totalXsec
    theoryPrediction.datasetPredictions = datasetPredictions
    theoryPrediction.conditions = None
    theoryPrediction.elements = elementList
    theoryPrediction.mass = mass
    theoryPrediction.totalwidth = totalwidth
    theoryPrediction.PIDs = [pdg for pdg, _ in itertools.groupby(PIDList)]  # Remove duplicates

    return theoryPrediction


def _getBestResult(dataSetResults):
    """
    Returns the best result according to the expected upper limit

    :param datasetPredictions: list of TheoryPredictionList objects
    :return: best result (TheoryPrediction object)
    """

    # In the case of UL analyses or efficiency-maps with a single signal region
    # return the single result:
    if len(dataSetResults) == 1:
        return dataSetResults[0]

    # For efficiency-map analyses with multipler signal regions,
    # select the best one according to the expected upper limit:
    bestExpectedR = 0.0
    bestXsec = 0.0 * fb
    for predList in dataSetResults:
        if len(predList) != 1:
            logger.error("Multiple clusters should only exist for upper limit results!")
            raise SModelSError()
        dataset = predList[0].dataset
        if dataset.getType() != "efficiencyMap":
            txt = (
                "Multiple data sets should only exist for efficiency map results, but we have them for %s?"
                % (predList[0].analysisId())
            )
            logger.error(txt)
            raise SModelSError(txt)
        pred = predList[0]
        xsec = pred.xsection
        expectedR = (xsec.value / dataset.getSRUpperLimit(0.05, True, False)).asNumber()
        if expectedR > bestExpectedR or (expectedR == bestExpectedR and xsec.value > bestXsec):
            bestExpectedR = expectedR
            bestPred = pred
            bestXsec = xsec.value

    return bestPred


def _getDataSetPredictions(dataset, smsTopList, maxMassDist, marginalize=False, deltas_rel=None):
    """
    Compute theory predictions for a given data set.
    For upper-limit results returns the list of theory predictions for the
    experimental result.  For efficiency-map results returns the list of theory
    predictions for the signal region.  Uses the list of elements in
    smsTopList.
    For each Txname appearing in dataset, it collects the elements and efficiencies,
    combine the masses (if needed) and compute the conditions (if existing).

    :parameter dataset: Data Set to be considered (DataSet object)
    :parameter smsTopList: list of topologies containing elements
                           (TopologyList object)
    :parameter maxMassDist: maximum mass distance for clustering elements (float)
    :returns:  a TheoryPredictionList object containing a list of TheoryPrediction
               objects
    """
    if deltas_rel is None:
        from smodels.tools.runtime import _deltas_rel_default

        deltas_rel = _deltas_rel_default

    predictionList = TheoryPredictionList()
    # Select elements belonging to expResult and apply efficiencies
    elements = _getElementsFrom(smsTopList, dataset)

    # Check dataset sqrts format:
    if (dataset.globalInfo.sqrts / TeV).normalize()._unit:
        ID = dataset.globalInfo.id
        logger.error("Sqrt(s) defined with wrong units for %s" % (ID))
        return False

    # Remove unwanted cross sections
    newelements = []
    for el in elements:
        el.weight = el.weight.getXsecsFor(dataset.globalInfo.sqrts)
        if not el.weight:
            continue
        newelements.append(el)
    elements = newelements
    if len(elements) == 0:
        return None

    # Combine elements according to their respective constraints and masses
    # (For efficiencyMap analysis group all elements)
    clusters = _combineElements(elements, dataset, maxDist=maxMassDist)

    # Collect results and evaluate conditions
    for cluster in clusters:
        theoryPrediction = TheoryPrediction(marginalize, deltas_rel)
        theoryPrediction.dataset = dataset
        theoryPrediction.txnames = cluster.txnames
        theoryPrediction.xsection = _evalConstraint(cluster)
        theoryPrediction.conditions = _evalConditions(cluster)
        theoryPrediction.elements = cluster.elements
        theoryPrediction.avgElement = cluster.averageElement()
        theoryPrediction.mass = theoryPrediction.avgElement.mass
        theoryPrediction.totalwidth = theoryPrediction.avgElement.totalwidth
        PIDs = [el.pdg for el in cluster.elements]
        theoryPrediction.PIDs = [pdg for pdg, _ in itertools.groupby(PIDs)]  # Remove duplicates
        predictionList._theoryPredictions.append(theoryPrediction)

    if len(predictionList) == 0:
        return None
    else:
        return predictionList


def _getElementsFrom(smsTopList, dataset):
    """
    Get elements that belong to any of the TxNames in dataset
    (appear in any of constraints in the result).
    Loop over all elements in smsTopList and returns a copy of the elements belonging
    to any of the constraints (i.e. have efficiency != 0). The copied elements
    have their weights multiplied by their respective efficiencies.

    :parameter dataset:  Data Set to be considered (DataSet object)
    :parameter smsTopList: list of topologies containing elements
                           (TopologyList object)
    :returns: list of elements (Element objects)
    """

    elements = []
    for txname in dataset.txnameList:
        for top in smsTopList:
            itop = txname._topologyList.index(top)  # Check if the topology appear in txname
            if itop is None:
                continue
            for el in top.getElements():
                newEl = txname.hasElementAs(el)  # Check if element appears in txname
                if not newEl:
                    continue
                el.setCoveredBy(dataset.globalInfo.type)
                eff = txname.getEfficiencyFor(newEl)
                if eff == None or abs(eff) < 1e-14:
                    continue
                el.setTestedBy(dataset.globalInfo.type)
                newEl.eff = eff
                newEl.weight *= eff
                newEl.txname = txname
                elements.append(newEl)  # Save element with correct branch ordering

    return elements


def _combineElements(elements, dataset, maxDist):
    """
    Combine elements according to the data set type.
    If expResult == upper limit type, first group elements with different TxNames
    and then into mass clusters.
    If expResult == efficiency map type, group all elements into a single cluster.

    :parameter elements: list of elements (Element objects)
    :parameter expResult: Data Set to be considered (DataSet object)
    :returns: list of element clusters (ElementCluster objects)
    """

    clusters = []

    if dataset.getType() == "efficiencyMap":  # cluster all elements
        clusters += clusterTools.clusterElements(elements, maxDist, dataset)
    elif dataset.getType() == "upperLimit":  # Cluster each txname individually
        txnames = list(set([el.txname for el in elements]))
        for txname in txnames:
            txnameEls = [el for el in elements if el.txname == txname]
            clusters += clusterTools.clusterElements(txnameEls, maxDist, dataset)
    else:
        logger.warning("Unkown data type: %s. Data will be ignored." % dataset.getType())

    return clusters


def _evalConstraint(cluster):
    """
    Evaluate the constraint inside an element cluster.
    If the cluster refers to a specific TxName, sum all the elements' weights
    according to the analysis constraint.
    For efficiency map results, sum all the elements' weights.

    :parameter cluster: cluster of elements (ElementCluster object)
    :returns: cluster cross section
    """

    if cluster.getDataType() == "efficiencyMap":
        return cluster.getTotalXSec()
    elif cluster.getDataType() == "upperLimit":
        if len(cluster.txnames) != 1:
            logger.error("An upper limit cluster should never contain more than one TxName")
            raise SModelSError()
        txname = cluster.txnames[0]
        if not txname.constraint or txname.constraint == "not yet assigned":
            return txname.constraint
        exprvalue = _evalExpression(txname.constraint, cluster)
        return exprvalue
    else:
        logger.error("Unknown data type %s" % (str(cluster.getDataType())))
        raise SModelSError()


def _evalConditions(cluster):
    """
    Evaluate the conditions (if any) inside an element cluster.

    :parameter cluster: cluster of elements (ElementCluster object)
    :returns: list of condition values (floats) if analysis type == upper limit. None, otherwise.
    """

    conditionVals = {}
    for txname in cluster.txnames:
        if not txname.condition or txname.condition == "not yet assigned":
            continue
        # Make sure conditions is always a list
        if isinstance(txname.condition, str):
            conditions = [txname.condition]
        elif isinstance(txname.condition, list):
            conditions = txname.condition
        else:
            logger.error("Conditions should be a list or a string")
            raise SModelSError()

        # Loop over conditions
        for cond in conditions:
            exprvalue = _evalExpression(cond, cluster)
            if isinstance(exprvalue, crossSection.XSection):
                conditionVals[cond] = exprvalue.value
            else:
                conditionVals[cond] = exprvalue

    if not conditionVals:
        return None
    else:
        return conditionVals


def _evalExpression(stringExpr, cluster):
    """
    Auxiliary method to evaluate a string expression using the weights of the elements in the cluster.
    Replaces the elements in stringExpr (in bracket notation) by their weights and evaluate the
    expression.

    e.g. computes the total weight of string expressions such as "[[[e^+]],[[e^-]]]+[[[mu^+]],[[mu^-]]]"
    or ratios of weights of string expressions such as "[[[e^+]],[[e^-]]]/[[[mu^+]],[[mu^-]]]"
    and so on...

    :parameter stringExpr: string containing the expression to be evaluated
    :parameter cluster: cluster of elements (ElementCluster object)
    :returns: xsection for the expression. Can be a XSection object, a float or not numerical (None,string,...)

    """

    # Get model for final state particles (database particles):
    model = cluster.dataset.globalInfo._databaseParticles
    # Get txname final state:
    if not hasattr(cluster.txnames[0], "finalState"):
        finalState = ["MET", "MET"]
    else:
        finalState = cluster.txnames[0].finalState
    if not hasattr(cluster.txnames[0], "intermediateState"):
        intermediateState = None
    else:
        intermediateState = cluster.txnames[0].intermediateState

    # Get cross section info from cluster (to generate zero cross section values):
    infoList = cluster.elements[0].weight.getInfo()
    # Get weights for elements appearing in stringExpr
    weightsDict = {}
    evalExpr = stringExpr.replace("'", "").replace(" ", "")
    for i, elStr in enumerate(elementsInStr(evalExpr)):
        el = element.Element(
            elStr, intermediateState=intermediateState, finalState=finalState, model=model
        )
        weightsDict["w%i" % i] = crossSection.XSectionList(infoList)
        for el1 in cluster.elements:
            if el1 == el:
                weightsDict["w%i" % i] += el1.weight
        evalExpr = evalExpr.replace(elStr, "w%i" % i)

    weightsDict.update({"Cgtr": cGtr, "cGtr": cGtr, "cSim": cSim, "Csim": cSim})
    exprvalue = eval(evalExpr, weightsDict)
    if type(exprvalue) == type(crossSection.XSectionList()):
        if len(exprvalue) != 1:
            logger.error("Evaluation of expression " + evalExpr + " returned multiple values.")
        return exprvalue[0]  # Return XSection object
    return exprvalue
