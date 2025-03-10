"""
.. module:: theoryPrediction
   :synopsis: Provides a class to encapsulate the results of the computation of
              reference cross sections and related functions.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

from smodels.base.physicsUnits import TeV, fb
from smodels.experiment.datasetObj import CombinedDataSet
from smodels.experiment.databaseObj import Database
from smodels.matching.exceptions import SModelSMatcherError as SModelSError
from smodels.matching import clusterTools
from smodels.base.smodelsLogging import logger
from typing import Union, Text, Dict
import numpy as np

__all__ = [ "TheoryPrediction", "theoryPredictionsFor", "TheoryPredictionsCombiner" ]

class TheoryPrediction(object):
    """
    An instance of this class represents the results of the theory prediction
    for an analysis.
    """

    def __init__(self, deltas_rel=None):
        """
        Initialize the theory prediction object. deltas_rel is meant to be a constant.

        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        """
        self.analysis = None
        self.xsection = None
        self.conditions = None
        self.mass = None
        self.totalwidth = None
        if deltas_rel is None:
            from smodels.base.runtime import _deltas_rel_default
            deltas_rel = _deltas_rel_default
        self.deltas_rel = deltas_rel
        self.cachedObjs = {False: {}, True: {}, "posteriori": {}}
        self.cachedNlls = {False: {}, True: {}, "posteriori": {}}
        self._statsComputer = None

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
            allConditions += cond

        if not allConditions:
            self.conditions = None
        else:
            self.conditions = allConditions[:]

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
            if value == "N/A" or value is None:
                continue
            values.append(float(value))
        return max(values)

    @property
    def statsComputer(self):
        if self._statsComputer is None:
            self.setStatsComputer()
        return self._statsComputer

    def setStatsComputer(self):
        """
        Creates and instance of StatsComputer depending on the
        type of TheoryPrediction/dataset. In case it is not possible
        to define a statistical computer (upper limit result or no expected
        upper limits), set the computer to 'N/A'.
        """
        from smodels.statistics.statsTools import getStatsComputerModule
        StatsComputer = getStatsComputerModule()

        if self.dataType() == "upperLimit":
            from smodels.base.runtime import experimentalFeature
            if not experimentalFeature( "truncatedGaussians" ):
                computer = 'N/A'
            else:
                computer = StatsComputer.forTruncatedGaussian(self)
                if computer is None: # No expected UL available
                    computer = 'N/A'

        elif self.dataType() == "efficiencyMap":
            nsig = (self.xsection * self.dataset.getLumi()).asNumber()
            computer = StatsComputer.forSingleBin(dataset=self.dataset,
                                                  nsig=nsig,
                                                  deltas_rel=self.deltas_rel )

        elif self.dataType() == "combined":
            # Get dictionary with dataset IDs and signal yields
            srNsigDict = {ds.getID() : 0.0 for ds in self.dataset.origdatasets}
            # Update with theory predictions
            srNsigDict.update({pred.dataset.getID() :
                          (pred.xsection*pred.dataset.getLumi()).asNumber()
                          for pred in self.datasetPredictions})

            # Get ordered list of datasets:
            if hasattr(self.dataset.globalInfo, "covariance"):
                datasetList = self.dataset.globalInfo.datasetOrder[:]
                # Get list of signal yields corresponding to the dataset order:
                srNsigs = [srNsigDict[dataID] for dataID in datasetList]
                # Get computer
                computer = StatsComputer.forMultiBinSL(dataset=self.dataset,
                                                       nsig=srNsigs,
                                                       deltas_rel = self.deltas_rel)

            elif hasattr(self.dataset.globalInfo, "jsonFiles"):
                computer = StatsComputer.forPyhf(dataset=self.dataset,
                                                       nsig=srNsigDict,
                                                       deltas_rel = self.deltas_rel)

        self._statsComputer = computer

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
            ul = None
            if self.dataType() == "efficiencyMap":
                ul = self.dataset.getSRUpperLimit(expected=expected)
            if self.dataType() == "upperLimit":
                ul = self.dataset.getUpperLimitFor(
                    sms=self.avgSMS, txnames=self.txnames, expected=expected
                )
            if self.dataType() == "combined":
                ul = self.statsComputer.poi_upper_limit(expected = expected,
                                                        limit_on_xsec = True)
            self.cachedObjs[expected]["UL"] = ul

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
        xsec = self.totalXsection()
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
                r = (self.totalXsection()/upperLimit).asNumber()
                self.cachedObjs[expected]["r"] = r
                return r
        return self.cachedObjs[expected]["r"]

    def whenDefined(function):
        """
        Returns the function whenever the statistical
        calculation is possible (i.e. when it is possible to define
        self.StatsComputer)
        """

        def wrapper(self, *args, **kwargs):
            if self.statsComputer == 'N/A':
                return None
            else:
                return function(self, *args, **kwargs)

        return wrapper

    @whenDefined
    def lsm(self, expected=False, return_nll : bool = False ):
        """likelihood at SM point, same as .def likelihood( ( mu = 0. )"""
        if "nll_sm" not in self.cachedObjs[expected]:
            self.computeStatistics(expected)
        if "nll_sm" not in self.cachedObjs[expected]:
            self.cachedObjs[expected]["lsm"] = None
        return self.nllToLikelihood ( self.cachedObjs[expected]["nll_sm"],
               return_nll )

    @whenDefined
    def lmax(self, expected=False, return_nll : bool = False ):
        """likelihood at mu_hat"""

        if not "nllmax" in self.cachedObjs[expected]:
            self.computeStatistics(expected)
        return self.nllToLikelihood ( self.cachedObjs[expected]["nllmax"],
                return_nll )

    @whenDefined
    def CLs(self, mu : float = 1., expected : Union[Text,bool] = False ) -> \
                    Union[float,None]:
        """ obtain the CLs value of the combination for a given poi value "mu" """
        if not "CLs" in self.cachedObjs[expected]:
            self.cachedObjs[expected]["CLs"] = {}
        if mu in self.cachedObjs[expected]["CLs"]:
            return self.cachedObjs[expected]["CLs"][mu]
        cls = self.statsComputer.CLs ( poi_test = mu, expected = expected )
        self.cachedObjs[expected]["CLs"][mu] = cls
        return cls

    @whenDefined
    def sigma_mu(self, expected=False):
        """sigma_mu of mu_hat"""

        if not "sigma_mu" in self.cachedObjs[expected]:
            self.computeStatistics(expected)

        return self.cachedObjs[expected]["sigma_mu"]

    @whenDefined
    def muhat(self, expected=False):
        """position of maximum likelihood"""

        if not "muhat" in self.cachedObjs[expected]:
            self.computeStatistics(expected)

        return self.cachedObjs[expected]["muhat"]

    @whenDefined
    def likelihood(self, mu=1.0, expected=False, return_nll=False, useCached=True):
        """
        get the likelihood for a signal strength modifier mu
        :param expected: compute expected, not observed likelihood. if "posteriori",
                         compute expected posteriori.
        :param return_nll: if True, return negative log likelihood, else likelihood
        :param useCached: if True, will return the cached value, if available
        """
        if useCached and mu in self.cachedNlls[expected]:
            nll = self.cachedNlls[expected][mu]
            return self.nllToLikelihood ( nll, return_nll )

        if useCached:
            if "nll" in self.cachedObjs[expected] and abs(mu - 1.0) < 1e-5:
                nll = self.cachedObjs[expected]["nll"]
                return self.nllToLikelihood ( nll, return_nll )
            if "nll_sm" in self.cachedObjs[expected] and abs(mu) < 1e-5:
                nllsm = self.cachedObjs[expected]["nll_sm"]
                return self.nllToLikelihood ( nllsm, return_nll )

        # for truncated gaussians the fits only work with negative signals!
        nll = self.statsComputer.likelihood(poi_test = mu,
                                             expected = expected,
                                             return_nll = True )
        self.cachedNlls[expected][mu] = nll

        if abs(mu) < 1e-5:
            self.cachedObjs[expected]["nll_sm"] = nll

        return self.nllToLikelihood ( nll, return_nll )

    def nllToLikelihood ( self, nll : Union[None,float], return_nll : bool ):
        """ if not return_nll, then compute likelihood from nll """
        if return_nll:
            return nll
        return np.exp ( - nll ) if nll is not None else None

    @whenDefined
    def computeStatistics(self, expected=False):
        """
        Compute the likelihoods, and upper limit for this theory prediction.
        The resulting values are stored as the likelihood, lmax, and lsm
        attributes.
        :param expected: computed expected quantities, not observed
        """

        if not "lmax" in self.cachedObjs[expected]:
            self.cachedObjs[expected]["lmax"] = {}
            self.cachedObjs[expected]["muhat"] = {}
            self.cachedObjs[expected]["sigma_mu"] = {}

        # Compute likelihoods and related parameters:
        llhdDict = self.statsComputer.get_five_values(expected = expected,
                     return_nll = True )
        if llhdDict not in [ None, {} ]:
            self.cachedObjs[expected]["nll"] = llhdDict["lbsm"]
            self.cachedObjs[expected]["nll_sm"] = llhdDict["lsm"]
            self.cachedObjs[expected]["nllmax"] = llhdDict["lmax"]
            self.cachedObjs[expected]["muhat"] = llhdDict["muhat"]
            if "sigma_mu" in llhdDict:
                self.cachedObjs[expected]["sigma_mu"] = llhdDict["sigma_mu"]


class TheoryPredictionsCombiner(TheoryPrediction):
    """
    Facility used to combine theory predictions from different analyes.
    If a list with a single TheoryPrediction is given, return the TheoryPrediction
    object.
    """

    def __new__(cls,theoryPredictions: list, slhafile=None, deltas_rel=None):
        """
        If called with a list containing a single TheoryPrediction, return the TheoryPrediction object.
        Otherwise, create a TheoryPredictionsCombiner object.
        """

        if len(theoryPredictions) == 1:
            return theoryPredictions[0]
        else:
            tpCombiner = super(TheoryPredictionsCombiner, cls).__new__(cls)
            return tpCombiner

    def __init__(self, theoryPredictions: list, slhafile=None, deltas_rel=None):
        """
        Constructor.

        :param theoryPredictions: the List of theory predictions
        :param slhafile: optionally, the slhafile can be given, for debugging
        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        """

        if len(theoryPredictions) == 0:
            raise SModelSError("asking to combine zero theory predictions")

        self.theoryPredictions = theoryPredictions
        self.slhafile = slhafile
        if deltas_rel is None:
            from smodels.base.runtime import _deltas_rel_default

            deltas_rel = _deltas_rel_default
        self.deltas_rel = deltas_rel
        self.cachedObjs = {False: {}, True: {}, "posteriori": {}}
        self.cachedNlls = {False: {}, True: {}, "posteriori": {}}
        self._statsComputer = None

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
            selectedTPs, key=lambda tp: (priority[tp.dataType()],
                                         tp.getRValue(expected=True) is not None,
                                         tp.getRValue(expected=True))
        )
        # Now get a single TP for each result
        # (the highest ranking analyses with r != None come last and are kept in the dict)
        uniqueTPs = {tp.analysisId(): tp for tp in selectedTPs}
        uniqueTPs = list(uniqueTPs.values())

        combiner = cls(uniqueTPs)
        return combiner

    def dataId(self):
        """
        Return a string with the IDs of all the datasets used in the combination.
        """
        ids = [str(tp.dataset.getID()) for tp in self.theoryPredictions]
        ret = ",".join(ids)

        return ret

    def analysisId(self):
        """
        Return a string with the IDs of all the experimental results used in the combination.
        """

        ret = ",".join(sorted([tp.analysisId() for tp in self.theoryPredictions]))

        return ret

    def dataType(self, short=False):
        """
        Return its type (combined)
        :param: short, if True, return abbreviation (anacomb)
        """
        if short:
            return "comb"
        else:
            return "combined"

    def totalXsection(self):
        ret = 0.0 * fb
        if self.theoryPredictions is not None:
            for tp in self.theoryPredictions:
                ret += tp.xsection
        return ret

    def getmaxCondition(self):
        """
        Returns the maximum xsection from the list conditions

        :returns: maximum condition xsection (float)
        """
        conditions = [tp.getmaxCondition() for tp in self.theoryPredictions]
        return max(conditions)

    def setStatsComputer(self):
        """
        Creates and instance of StatsComputer depending on the
        type of TheoryPrediction/dataset. In case it is not possible
        to define a statistical computer (upper limit result or no expected
        upper limits), set the computer to 'N/A'.
        """

        # First make sure all theory predictions in the combiner
        # have well-defined stats models
        if any(tp.statsComputer == 'N/A' for tp in self.theoryPredictions):
            computer = 'N/A'
        else:
            from smodels.statistics.statsTools import getStatsComputerModule
            StatsComputer = getStatsComputerModule()
            computer = StatsComputer.forAnalysesComb(self.theoryPredictions, self.deltas_rel)

        self._statsComputer = computer

    def getLlhds(self,muvals,expected=False,normalize=True):
        """
        Facility to access the likelihoods for the individual analyses and the combined
        likelihood.
        Returns a dictionary with the analysis IDs as keys and the likelihood values as values.
        Mostly used for plotting the likelihoods.

        :param muvals: List with values for the signal strenth for which the likelihoods must
                       be evaluated.
        :param expected: If True returns the expected likelihood values.
        :param normalize: If True normalizes the likelihood by its integral over muvals.
        """

        return self.statsComputer.likelihoodComputer.getLlhds(muvals,expected,normalize)

    def describe(self):
        """returns a string containing a list of all analysisId and dataIds"""
        ids = []
        for tp in self.theoryPredictions:
            ids.append(f"{tp.analysisId()}:{tp.dataId()}")
        return f"SRs: {', '.join(ids)}"


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
        if theoryPredictions and isinstance(theoryPredictions,
                                            (TheoryPredictionList,list)):

            if not maxCond:
                self._theoryPredictions = theoryPredictions
            else:
                newPredictions = []
                for theoPred in theoryPredictions:
                    mCond = theoPred.getmaxCondition()
                    if mCond != "N/A" and round(mCond/maxCond, 2) > 1.0:
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

    def removeRNone(self):
        """
        Remove predictions for which r-value = None
        (such as when the UL computer fails due to convergence issues).
        """

        tpList = [tp for tp in self._theoryPredictions
                  if tp.getRValue() is not None]
        self._theoryPredictions = tpList[:]

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


def theoryPredictionsFor(database : Database, smsTopDict : Dict,
        maxMassDist : float = 0.2, useBestDataset : bool = True,
        combinedResults : bool = True, deltas_rel : Union[None,float] = None):
    """
    Compute theory predictions for the given experimental result, using the list of
    SMS in smsTopDict.
    For each Txname appearing in expResult, it collects the SMS and
    efficiencies, combine the SMS and compute the conditions  (if exist).

    :parameter database: the database with the selected experimental results
    :parameter smsTopDict: dictionary of SMS, where the canonical names are keys and the TheorySMS objects are values.
                           (TopologyDict object)
    :parameter maxMassDist: maximum mass distance for clustering SMS (float)
    :parameter useBestDataset: If True, uses only the best dataset (signal region).
               If False, returns predictions for all datasets (if combinedResults is False),
               or only the combinedResults (if combinedResults is True).
    :parameter combinedResults: add theory predictions that result from
               combining datasets.
    :parameter deltas_rel: relative uncertainty in signal (float). Default value is 20%.

    :returns:  a TheoryPredictionList object containing a list of TheoryPrediction
               objects
    """
    if deltas_rel is None:
        from smodels.base.runtime import _deltas_rel_default
        deltas_rel = _deltas_rel_default

    if not isinstance(database,Database):
        errorMsg = "The argument for theoryPredictionsFor should be a"
        errorMsg += " database with the selected experimental results and not %s" %type(database)
        logger.error(errorMsg)
        raise SModelSError(errorMsg)

    # Compute matches between TheorySMS from decomposition and the
    # (unique) ExpSMS in the database
    expSMSDict = database.expSMSDict
    # Compute dictionary with matches:
    smsMatch = expSMSDict.getMatchesFrom(smsTopDict)

    ret = []
    for expResult in database.expResultList:
        dataSetResults = []
        #  Compute predictions for each data set (for UL analyses there is one single set)
        for dataset in expResult.datasets:
            predList = _getDataSetPredictions(dataset, smsMatch, expSMSDict, maxMassDist)
            if predList:
                dataSetResults.append(predList)
        if not dataSetResults:  # no results at all?
            continue

        #  For results with more than one dataset keep all dataset predictions
        if len(dataSetResults) == 1:  # only a single dataset? Easy case.
            expResults = dataSetResults[0]
        else:
            if combinedResults: # Include combination
                combinedRes = _getCombinedResultFor(dataSetResults,
                                                expResult)
                if combinedRes is not None:
                    dataSetResults.append(TheoryPredictionList([combinedRes]))
            if not useBestDataset: #  Report all datasets
                expResults = sum(dataSetResults)
            else:
                expResults = TheoryPredictionList()
                bestRes = _getBestResult(dataSetResults)
                if not bestRes is None:
                    expResults.append(bestRes) # Best result = combination if available

        for theoPred in expResults:
            theoPred.expResult = expResult
            theoPred.deltas_rel = deltas_rel
            tpe = None
            if isinstance(theoPred.dataset,CombinedDataSet): # Individual CRs shouldn't give results
                theoPred.upperLimit = theoPred.getUpperLimit()
                continue
            else:
                if hasattr(theoPred.dataset.globalInfo, "jsonFiles"): # Only signal in CRs for jsonFiles so far
                    for regionSet in theoPred.dataset.globalInfo.jsonFiles.values():
                        for region in regionSet:
                            if region["smodels"] == theoPred.dataset.dataInfo.dataId:
                                tpe = region["type"]
                                break
                else:
                    tpe = "SR"

                if tpe is None:
                    logger.error(f"Could not find type of region {theoPred.dataType()} from {theoPred.analysisId()}")
                    raise SModelSError()

                if tpe == "SR":
                    theoPred.upperLimit = theoPred.getUpperLimit()
                else:
                    theoPred.upperLimit = None

        expResults.sortTheoryPredictions()

        for theoPred in expResults:
            ret.append(theoPred)

    tpList = TheoryPredictionList(ret)
    # tpList.removeRNone() # Remove results with r = None
    tpList.sortTheoryPredictions()

    return tpList


def _getCombinedResultFor(dataSetResults, expResult):
    """
    Compute the combined result for all datasets, if covariance
    matrices or jsonFiles are available. Return a TheoryPrediction object
    with the signal cross-section summed over all the signal regions
    and the respective upper limit.

    :param datasetPredictions: List of TheoryPrediction objects for each signal region
    :param expResult: ExpResult object corresponding to the experimental result

    :return: TheoryPrediction object
    """
    # Don't give combined result if all regions are CRs
    isNotSR = []
    for predList in dataSetResults:
        if hasattr ( expResult.globalInfo, "jsonFiles" ):
            for regionSet in expResult.globalInfo.jsonFiles.values():
                for region in regionSet:
                    if region['smodels'] == predList[0].dataset.dataInfo.dataId:
                        if region['type'] == 'SR':
                            isNotSR.append(False)
                        else:
                            isNotSR.append(True)
        else:
            isNotSR = [ False ]

    if all(isNotSR):
        return None

    if len(dataSetResults) == 1:
        return dataSetResults[0]
    elif not expResult.hasCovarianceMatrix() and not expResult.hasJsonFile():
        return None

    txnameList = []
    smsList = []
    totalXsec = None
    datasetPredictions = []
    avgSMSlist = []
    for predList in dataSetResults:
        if len(predList) != 1:
            raise SModelSError("Results with multiple datasets should have a single theory prediction (EM-type)!")
        pred = predList[0]
        datasetPredictions.append(pred)
        txnameList += pred.txnames
        smsList += pred.smsList
        avgSMSlist.append(pred.avgSMS)
        if totalXsec is None:
            totalXsec = pred.xsection
        else:
            totalXsec += pred.xsection

    txnameList = list(set(txnameList))

    # Check if all avgSMS are the same:
    uniqueSMS = all(avgSMSlist[0] == avgSMS for avgSMS in avgSMSlist[1:])
    # If so, keep it in the theoryPrediction, otherwise set to None
    if uniqueSMS:
        avgSMS = avgSMSlist[0]
    else:
        avgSMS = None

    # Create a combinedDataSet object:
    combinedDataset = CombinedDataSet(expResult)
    # Create a theory preidction object for the combined datasets:
    theoryPrediction = TheoryPrediction()
    theoryPrediction.dataset = combinedDataset
    theoryPrediction.txnames = txnameList
    theoryPrediction.avgSMS = avgSMS
    theoryPrediction.xsection = totalXsec
    theoryPrediction.datasetPredictions = datasetPredictions
    theoryPrediction.conditions = None
    theoryPrediction.smsList = smsList

    return theoryPrediction


def _getBestResult(dataSetResults):
    """
    Returns the best result according to the expected upper limit.
    If a combined result is included in the list, always return it.

    :param datasetPredictions: list of TheoryPredictionList objects
    :return: best result (TheoryPrediction object)
    """

    # In the case of UL analyses or efficiency-maps with a single signal region
    # return the single result:
    if len(dataSetResults) == 1:
        return dataSetResults[0]

    # If combination is included, always return it
    for predList in dataSetResults:
        for tp in predList:
            if isinstance(tp.dataset,CombinedDataSet):
                return tp


    # For efficiency-map analyses with multipler signal regions,
    # select the best one according to the expected upper limit:
    bestExpectedR = 0.0
    bestXsec = 0.0*fb
    bestPred = None
    for predList in dataSetResults:
        if len(predList) != 1:
            logger.error("Multiple clusters should only exist for upper limit results!")
            raise SModelSError()
        dataset = predList[0].dataset
        globalInfo = dataset.globalInfo

        # Only a SR can be the best SR
        stop = False
        if hasattr(globalInfo,"jsonFiles"):
            for regionSet in globalInfo.jsonFiles.values():
                for region in regionSet:
                    if type(region) == dict and \
                            region['smodels'] == dataset.dataInfo.dataId:
                        if "type" in region and region['type'] != 'SR':
                            stop = True
                    if stop: break
                if stop: break
        if stop:
            continue

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
                           deltas_rel=None):
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
        theoryPrediction = TheoryPrediction(deltas_rel)
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
