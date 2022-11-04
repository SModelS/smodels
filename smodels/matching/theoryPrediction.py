"""
.. module:: theoryPrediction
   :synopsis: Provides a class to encapsulate the results of the computation of
              reference cross sections and related functions.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

from smodels.base.physicsUnits import TeV, fb
from smodels.experiment.datasetObj import CombinedDataSet
from smodels.matching.exceptions import SModelSMatcherError as SModelSError
from smodels.matching import clusterTools
from smodels.matching.matcherAuxiliaryFuncs import average
from smodels.base.smodelsLogging import logger
from smodels.statistics.statisticalTools import TruncatedGaussians
from smodels.statistics.srCombinations import getCombinedStatistics
from smodels.statistics.srCombinations import getCombinedUpperLimitFor, getCombinedLikelihood
import itertools
import numpy as np


class TheoryPrediction(object):
    """
    An instance of this class represents the results of the theory prediction
    for an analysis.
    """

    def __init__(self, marginalize=False,
                 deltas_rel=None):
        """ a theory prediction. marginalize and deltas_rel are meant to be
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
            from smodels.base.runtime import _deltas_rel_default
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
            D = {"upperLimit": "ul", "efficiencyMap": "em",
                 "combined": "comb"}
            if t in D.keys():
                return D[t]
            return "??"

        return self.dataset.getType()

    def computeXSection(self):

        xsection = 0*fb
        # Adds the contributions of all txnames
        # (for an UL result there will be a single txname)
        for tx in self.txnames:
            # Filter the SMS
            smsList = [sms for sms in self.smsList
                        if sms.txname is tx]
            # Build dictionary needed for evaluation:
            xsection += tx.evalConstraintFor(smsList)

        self.xsection = xsection

    def computeConditions(self):

        # Adds the contributions of all txnames
        # (for an UL result there will be a single txname
        # and for EM results there should be no conditions)
        allConditions = []
        for tx in self.txnames:
            # Filter the SMS
            smsList = [sms for sms in self.smsList
                        if sms.txname is tx]
            cond = tx.evalConditionsFor(smsList)
            if cond is None:
                continue
            allConditions.append(cond)

        if not allConditions:
            self.conditions = None
        elif len(allConditions) == 1:
            self.conditions = allConditions[0]
        else:
            msgError = "Multiple conditions found for a %s result." % self.dataType()
            logger.error(msgError)
            raise SModelSError(msgError)

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

        #  First check if the upper-limit and expected upper-limit have already been computed.
        #  If not, compute it and store them.
        if "UL" not in self.cachedObjs[expected]:
            if self.dataType() == 'efficiencyMap':
                ul = self.dataset.getSRUpperLimit(expected=expected)
                self.cachedObjs[expected]["UL"] = ul
            if self.dataType() == 'upperLimit':
                ul = self.dataset.getUpperLimitFor(sms=self.avgSMS,
                                                   txnames=self.txnames,
                                                   expected=expected)
                self.cachedObjs[expected]["UL"] = ul
            if self.dataType() == 'combined':
                #  Create a list of signal events in each dataset/SR sorted according to datasetOrder
                #  lumi = self.dataset.getLumi()
                if hasattr(self.dataset.globalInfo, "covariance"):
                    srNsigDict = dict([[pred.dataset.getID(), (pred.xsection*pred.dataset.getLumi()).asNumber()] for pred in self.datasetPredictions])
                    srNsigs = [srNsigDict[dataID] if dataID in srNsigDict else 0. for dataID in self.dataset.globalInfo.datasetOrder]
                elif hasattr(self.dataset.globalInfo, "jsonFiles"):
                    srNsigDict = dict([[pred.dataset.getID(), (pred.xsection*pred.dataset.getLumi()).asNumber()] for pred in self.datasetPredictions])
                    srNsigs = [srNsigDict[ds.getID()] if ds.getID() in srNsigDict else 0. for ds in self.dataset._datasets]
                self.cachedObjs[expected]["UL"] = getCombinedUpperLimitFor(self.dataset, srNsigs, expected=expected, deltas_rel=self.deltas_rel)

        return self.cachedObjs[expected]["UL"]

    def getUpperLimitOnMu(self, expected=False):
        """
        Get upper limit on signal strength multiplier, using the
        theory prediction value and the corresponding upper limit
        (i.e. mu_UL = upper limit/theory xsec)

        :param expected: if True, compute expected upper limit, else observed
        :returns: upper limit on signal strength multiplier mu
        """

        upperLimit = self.getUpperLimit(expected=expected)
        xsec = self.xsection
        if xsec is None or upperLimit is None:
            return None

        muUL = (upperLimit/xsec).asNumber()

        return muUL

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
                r = (self.xsection/upperLimit).asNumber()
                self.cachedObjs[expected]["r"] = r
                return r
        return self.cachedObjs[expected]["r"]

    def lsm(self, expected=False):
        """ likelihood at SM point, same as .def likelihood( ( mu = 0. ) """
        if "lsm" not in self.cachedObjs[expected]:
            self.computeStatistics(expected)
        if "lsm" not in self.cachedObjs[expected]:
            self.cachedObjs[expected]["lsm"] = None
        return self.cachedObjs[expected]["lsm"]

    def lmax(self, expected=False, allowNegativeSignals=False):
        """ likelihood at mu_hat """
        if "lmax" not in self.cachedObjs[expected]:
            self.computeStatistics(expected, allowNegativeSignals)
        if ("lmax" in self.cachedObjs[expected]
            and allowNegativeSignals not in self.cachedObjs[expected]["lmax"]):
            self.computeStatistics(expected, allowNegativeSignals)
        if "lmax" not in self.cachedObjs[expected]:
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
        """ position of maximum likelihood  """
        if "muhat" not in self.cachedObjs[expected]:
            self.computeStatistics(expected, allowNegativeSignals)
        if allowNegativeSignals not in self.cachedObjs[expected]["muhat"]:
            self.computeStatistics(expected, allowNegativeSignals)
        if ("muhat" in self.cachedObjs[expected]
	        and allowNegativeSignals not in self.cachedObjs[expected]["muhat"]):
            self.computeStatistics(expected, allowNegativeSignals)
        if "muhat" not in self.cachedObjs[expected]:
            self.cachedObjs[expected]["muhat"] = {allowNegativeSignals: None}
        ret = self.cachedObjs[expected]["muhat"][allowNegativeSignals]
        return ret

    def chi2(self, expected=False):
        if "chi2" not in self.cachedObjs[expected]:
            self.computeStatistics(expected)
        if "chi2" not in self.cachedObjs[expected]:
            self.cachedObjs[expected]["chi2"] = None
        return self.cachedObjs[expected]["chi2"]

    def likelihood(self, mu=1., expected=False,
                   nll=False, useCached=True):
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
            srNsigDict = dict([[pred.dataset.getID(), (pred.xsection*lumi).asNumber()] for
                              pred in self.datasetPredictions])
            srNsigs = [srNsigDict[ds.getID()] if ds.getID() in srNsigDict else 0.
                       for ds in self.dataset._datasets]
            llhd = getCombinedLikelihood(self.dataset, srNsigs,
                                             self.marginalize,
                                             self.deltas_rel,
                                             expected=expected, mu=mu)
        if self.dataType() == "efficiencyMap":
            nsig = (mu*self.xsection*lumi).asNumber()
            llhd = self.dataset.likelihood(nsig, marginalize=self.marginalize,
                                           deltas_rel=self.deltas_rel,
                                           expected=expected)

        if self.dataType() == "upperLimit":
            # these fits only work with negative signals!
            llhd, chi2 = self.likelihoodFromLimits(mu, chi2also=True,
                                                   expected=expected,
                                                   allowNegativeSignals=True)
        self.cachedLlhds[expected][mu] = llhd
        if nll:
            if llhd == 0.0:
                return 700.0
            return -np.log(llhd)
        return llhd

    def likelihoodFromLimits(self, mu=1.0, expected=False, chi2also=False,
                             corr=0.6, allowNegativeSignals=False):
        """ compute the likelihood from expected and observed upper limits.
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
        from smodels.base.runtime import experimentalFeatures

        if not experimentalFeatures():
            if chi2also:
                return (None, None)
            return None
        if not hasattr(self, "avgSMS"):
            logger.error("theory prediction %s has no average SMS! why??" % self.analysisId())
            if chi2also:
                return (None, None)
            return None

        eul = self.dataset.getUpperLimitFor(sms=self.avgSMS,
                                            txnames=self.txnames,
                                            expected=True)
        if eul is None:
            if chi2also:
                return (None, None)
            return None
        ul = self.dataset.getUpperLimitFor(sms=self.avgSMS,
                                           txnames=self.txnames,
                                           expected=False)
        lumi = self.dataset.getLumi()
        computer = TruncatedGaussians(ul, eul, self.xsection, lumi=lumi, corr = corr)
        ret = computer.likelihood(mu)
        llhd, muhat, sigma_mu = ret["llhd"], ret["muhat"], ret["sigma_mu"]

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
        if "lmax" not in self.cachedObjs[expected]:
            self.cachedObjs[expected]["lmax"] = {}
            self.cachedObjs[expected]["muhat"] = {}

        if self.dataType() == "upperLimit":
            llhd, chi2 = self.likelihoodFromLimits(1.0, expected=expected, chi2also=True)
            lsm = self.likelihoodFromLimits(0.0, expected=expected, chi2also=False)
            lmax = self.likelihoodFromLimits(None, expected=expected, chi2also=False,
                                             allowNegativeSignals=True)
            if allowNegativeSignals is False and hasattr(self, "muhat_") and self.muhat_ < 0.0:
                self.muhat_ = 0.0
                lmax = lsm
            self.cachedObjs[expected]["llhd"] = llhd
            self.cachedObjs[expected]["lsm"] = lsm
            self.cachedObjs[expected]["lmax"][allowNegativeSignals] = lmax
            self.cachedObjs[expected]["chi2"] = chi2
            if hasattr(self, "muhat_"):
                self.cachedObjs[expected]["muhat"][allowNegativeSignals] = self.muhat_
            if hasattr(self, "sigma_mu_"):
                if "sigma_mu" not in self.cachedObjs[expected]:
                    self.cachedObjs[expected]["sigma_mu"] = {}
                self.cachedObjs[expected]["sigma_mu"][allowNegativeSignals] = self.sigma_mu_

        elif self.dataType() == "efficiencyMap":
            lumi = self.dataset.getLumi()
            nsig = (self.xsection*lumi).asNumber()
            llhd = self.dataset.likelihood(nsig, marginalize=self.marginalize,
                                           deltas_rel=self.deltas_rel,
                                           expected=expected)
            llhd_sm = self.dataset.likelihood(nsig=0.0, marginalize=self.marginalize,
                                              deltas_rel=self.deltas_rel,
                                              expected=expected)
            llhd_max = self.dataset.lmax(marginalize=self.marginalize,
                                         deltas_rel=self.deltas_rel,
                                         allowNegativeSignals=allowNegativeSignals,
                                         expected=expected)
            muhat = None
            if hasattr(self.dataset, "muhat"):
                muhat = self.dataset.muhat / nsig
            if hasattr(self.dataset, "sigma_mu"):
                sigma_mu = float(self.dataset.sigma_mu/nsig)
                if not "sigma_mu" in self.cachedObjs[expected]:
                    self.cachedObjs[expected]["sigma_mu"] = {}
                self.cachedObjs[expected]["sigma_mu"][allowNegativeSignals] = sigma_mu
            self.cachedObjs[expected]["llhd"] = llhd
            self.cachedObjs[expected]["lsm"] = llhd_sm
            self.cachedObjs[expected]["lmax"][allowNegativeSignals] = llhd_max
            self.cachedObjs[expected]["muhat"][allowNegativeSignals] = muhat
            from smodels.statistics.statisticalTools import chi2FromLmax
            self.cachedObjs[expected]["chi2"] = chi2FromLmax(llhd, llhd_max)

        elif self.dataType() == "combined":
            lumi = self.dataset.getLumi()
            #  Create a list of signal events in each dataset/SR sorted according to datasetOrder
            srNsigDict = dict([[pred.dataset.getID(), (pred.xsection*lumi).asNumber()] for pred in self.datasetPredictions])
            srNsigs = [srNsigDict[ds.getID()] if ds.getID() in srNsigDict else 0.0 for ds in self.dataset._datasets]

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
            from smodels.statistics.statisticalTools import chi2FromLmax

            self.cachedObjs[expected]["chi2"] = chi2FromLmax(llhd, lmax)

    def totalXsection(self):
        return self.xsection

    def getmaxCondition(self):
        """
        Returns the maximum xsection from the list conditions

        :returns: maximum condition xsection (float)
        """

        if not self.conditions:
            return 0.0

        values = [0.0]
        for value in self.conditions:
            if value == 'N/A':
                return value
            if value is None:
                continue
            values.append(float(value))
        return max(values)


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
                    if mCond == "N/A" or round(mCond/maxCond, 2) > 1.0:
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
        self._theoryPredictions = sorted(self._theoryPredictions,
                                         key=lambda theoPred: (theoPred.getRValue() is not None, theoPred.getRValue()),
                                         reverse=True)


def theoryPredictionsFor(database, smsTopDict, maxMassDist=0.2,
                         useBestDataset=True, combinedResults=True,
                         marginalize=False, deltas_rel=None):
    """
    Compute theory predictions for the given experimental result, using the list of
    SMS in smsTopDict.
    For each Txname appearing in expResult, it collects the SMS and
    efficiencies, combine the SMS and compute the conditions  (if exist).

    :parameter expResult: expResult to be considered (ExpResult object), if list of ExpResults is given, produce theory predictions for all
    :parameter smsTopDict: dictionary of SMS, where the canonical names are keys and the TheorySMS objects are values.
                           (TopologyDict object)
    :parameter maxMassDist: maximum mass distance for clustering SMS (float)
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
        from smodels.base.runtime import _deltas_rel_default
        deltas_rel = _deltas_rel_default


    # Compute matches between TheorySMS from decomposition and the
    # (unique) ExpSMS in the database
    smsDict = database.expSMSDict
    smsMatch = {sms : [] for sms in smsDict}
    cNamesDict = {}
    for sms in smsDict:
        if sms.canonName not in cNamesDict:
            cNamesDict[sms.canonName] = []
        cNamesDict[sms.canonName].append(sms)

    for sms in smsTopDict.getSMSList():
        canonName  = sms.canonName
        # Select txSMS with matching canon name:
        selectedSMSlist = []
        for cName,smsList in cNamesDict.items():
            if cName == canonName:
                selectedSMSlist += smsList

        # Loop over selected SMS and check for matches
        # Store the (correctly orderd) match in smsMatch
        for txsms in selectedSMSlist:
            matchedSMS = txsms.matchesTo(sms)
            if matchedSMS is None:
                continue
            smsMatch[txsms].append((matchedSMS,sms))

    ret = []
    for expResult in database.expResultList:
        dataSetResults = []
        #  Compute predictions for each data set (for UL analyses there is one single set)
        for dataset in expResult.datasets:
            predList = _getDataSetPredictions(dataset, smsMatch, smsDict, maxMassDist)
            if predList:
                dataSetResults.append(predList)
        if not dataSetResults:  # no results at all?
            continue
        
        #  For results with more than one dataset keep all dataset predictions
        if len(dataSetResults) == 1:  # only a single dataset? Easy case.
            expResults = dataSetResults[0]
        
        #  if useBestDataSet=False and combinedResults=False:
        elif not useBestDataset and not combinedResults:
            expResults = sum(dataSetResults)

        elif combinedResults:  # Try to combine results
            expResults = TheoryPredictionList()
            combinedRes = _getCombinedResultFor(dataSetResults,
                                                expResult, marginalize)
            if combinedRes is None:  # Can not combine, use best dataset:
                combinedRes = _getBestResult(dataSetResults)
            expResults.append(combinedRes)
        else:  # Use best dataset:
            expResults = TheoryPredictionList()
            expResults.append(_getBestResult(dataSetResults))

        for theoPred in expResults:
            theoPred.expResult = expResult
            theoPred.deltas_rel = deltas_rel
            theoPred.upperLimit = theoPred.getUpperLimit()
        expResults.sortTheoryPredictions()

        for theoPred in expResults:
            ret.append(theoPred)

    tpList = TheoryPredictionList(ret)
    tpList.sortTheoryPredictions()

    return tpList


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
    smsList = []
    totalXsec = None
    datasetPredictions = []
    for predList in dataSetResults:
        if len(predList) != 1:
            raise SModelSError("Results with multiple datasets should have a single theory prediction (EM-type)!")
        pred = predList[0]
        datasetPredictions.append(pred)
        txnameList += pred.txnames
        smsList += pred.smsList
        if totalXsec is None:
            totalXsec = pred.xsection
        else:
            totalXsec += pred.xsection

    txnameList = list(set(txnameList))

    # Create a combinedDataSet object:
    combinedDataset = CombinedDataSet(expResult)
    combinedDataset._marginalize = marginalize
    # Create a theory preidction object for the combined datasets:
    theoryPrediction = TheoryPrediction(marginalize)
    theoryPrediction.dataset = combinedDataset
    theoryPrediction.txnames = txnameList
    theoryPrediction.avgSMS = None
    theoryPrediction.xsection = totalXsec
    theoryPrediction.datasetPredictions = datasetPredictions
    theoryPrediction.conditions = None
    theoryPrediction.smsList = smsList

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
    bestXsec = 0.0*fb
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
        expectedR = (xsec/dataset.getSRUpperLimit(expected=True)).asNumber()
        if expectedR > bestExpectedR or (expectedR == bestExpectedR and xsec > bestXsec):
            bestExpectedR = expectedR
            bestPred = pred
            bestXsec = xsec

    return bestPred


def _getDataSetPredictions(dataset, smsMatch,smsDict, maxMassDist,
                           marginalize=False, deltas_rel=None,new=False):
    """
    Compute theory predictions for a given data set.
    For upper-limit results returns the list of theory predictions for the
    experimental result.  For efficiency-map results returns the list of theory
    predictions for the signal region.  Uses the list of SMS in
    smsTopDict.
    For each Txname appearing in dataset, it collects the SMS and efficiencies,
    combine the masses (if needed) and compute the conditions (if existing).

    :parameter dataset: Data Set to be considered (DataSet object)
    :parameter smsTopDict: dictionary of SMS, where the canonical names are keys and the TheorySMS objects are values.
                           (TopologyDict object)
    :parameter maxMassDist: maximum mass distance for clustering SMS (float)
    :returns:  a TheoryPredictionList object containing a list of TheoryPrediction
               objects
    """
    if deltas_rel is None:
        from smodels.base.runtime import _deltas_rel_default
        deltas_rel = _deltas_rel_default

    predictionList = TheoryPredictionList()
    #  Select SMS belonging to expResult and apply efficiencies
    smsList = _getSMSFor(dataset,smsMatch,smsDict)

    if len(smsList) == 0:
        return None

    #  Check dataset sqrts format:
    if (dataset.globalInfo.sqrts/TeV).normalize()._unit:
        ID = dataset.globalInfo.id
        logger.error("Sqrt(s) defined with wrong units for %s" % (ID))
        return False

    # Compute relevant SMS weights
    newList = []
    for sms in smsList:
        # Get cross-sections for correct CM energy:
        sms.weight = sms.weightList.getXsecsFor(dataset.globalInfo.sqrts)
        if not sms.weight:
            continue
        # Get largest weight (in case there are LO, NLO,... values)
        sms.weight = sms.weight.getMaxXsec()
        # Multiply weight by SMSs efficiency:
        sms.weight = sms.weight*sms.eff
        newList.append(sms)
    smsList = newList[:]


    #  Combine SMS according to their respective constraints and masses
    #  (For efficiencyMap analysis group all SMS)
    clusters = _combineSMS(smsList, dataset, maxDist=maxMassDist)

    #  Collect results and evaluate conditions
    for cluster in clusters:
        theoryPrediction = TheoryPrediction(marginalize, deltas_rel)
        theoryPrediction.dataset = dataset
        theoryPrediction.txnames = cluster.txnames
        theoryPrediction.smsList = cluster.smsList
        theoryPrediction.avgSMS = cluster.averageSMS
        # Compute relevant cross-section and conditions:
        theoryPrediction.computeXSection()
        # Skip results with too small (invisible) cross-sections
        if theoryPrediction.xsection < 1e-6*fb:
            continue
        theoryPrediction.computeConditions()

        predictionList._theoryPredictions.append(theoryPrediction)

    if len(predictionList) == 0:
        return None
    else:
        return predictionList


def _getSMSFor(dataset,smsMatch,smsDict):
    """
    Get SMS that belong to any of the TxNames in dataset
    (appear in any of constraints in the result).

    :parameter dataset:  Data Set to be considered (DataSet object)
    :parameter smsMatch: dictionary with unique ExpSMS as keys and the corresponding list of
                         (matched TheorySMS, orignal TheorySMS) as values
    :returns: list of SMS (TheorySMS objects)
    """

    smsList = []
    for txname in dataset.txnameList:
        # Loop over unique SMS for the given txname:
        for smsLabel,txsms in smsDict[txname].items():
            for sms,sms_orig in smsMatch[txsms]:
                # Tag the original SMS as covered:
                sms_orig.setCoveredBy(dataset.globalInfo.type)
                newSMS = sms.copy()
                newSMS = smsDict.setTxNodeOrdering(newSMS,txname,smsLabel)
                # Compute efficiency
                eff = txname.getEfficiencyFor(newSMS)
                if eff is None or abs(eff) < 1e-14:
                    continue
                # Tag the original SMS as tested:
                sms_orig.setTestedBy(dataset.globalInfo.type)
                newSMS.eff = eff
                newSMS.txname = txname
                newSMS.txlabel = smsLabel             
                smsList.append(newSMS)

    return smsList


def _combineSMS(smsList, dataset, maxDist):
    """
    Combine SMS according to the data set type.
    If expResult == upper limit type, first group SMS with different TxNames
    and then into mass clusters.
    If expResult == efficiency map type, group all SMS into a single cluster.

    :parameter smsList: list of SMS (TheorySMS objects)
    :parameter expResult: Data Set to be considered (DataSet object)
    :returns: list of SMS clusters (SMSCluster objects)
    """

    clusters = []

    if dataset.getType() == "efficiencyMap":  # cluster all SMS
        clusters += clusterTools.clusterSMS(smsList, maxDist, dataset)
    elif dataset.getType() == "upperLimit":  # Cluster each txname individually
        txnames = list(set([sms.txname for sms in smsList]))
        for txname in txnames:
            txnameSMS = [sms for sms in smsList if sms.txname is txname]
            clusters += clusterTools.clusterSMS(txnameSMS, maxDist, dataset)
    else:
        logger.warning("Unkown data type: %s. Data will be ignored."
                       % dataset.getType())

    return clusters
