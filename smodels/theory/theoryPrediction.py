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
from smodels.tools.statsTools import StatsComputer
from typing import Union
import itertools
import numpy as np

__all__ = [ "TheoryPrediction", "theoryPredictionsFor", "TheoryPredictionsCombiner" ]


class TheoryPrediction(object):
    """
    An instance of this class represents the results of the theory prediction
    for an analysis.
    """

    def __init__(self, deltas_rel=None):
        """a theory prediction. deltas_rel is meant to be a constant
        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        """
        self.analysis = None
        self.xsection = None
        self.conditions = None
        self.mass = None
        self.totalwidth = None
        if deltas_rel is None:
            from smodels.tools.runtime import _deltas_rel_default

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
            D = {"upperLimit": "ul", "efficiencyMap": "em", "combined": "comb"}
            if t in D.keys():
                return D[t]
            return "??"

        return self.dataset.getType()

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
            if value == "N/A" or value is None:
                continue
            # print ( "value=",value,type(value),float(value) )
            # maxcond = max(maxcond,float(value))
            values.append(float(value))
        return max(values)
        # return maxcond

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

        if self.dataType() == "upperLimit":
            from smodels.tools.runtime import experimentalFeatures
            if not experimentalFeatures():
                computer = 'N/A'
            else:
                computer = StatsComputer.forTruncatedGaussian(self)
                if computer is None: # No expected UL available
                    computer = 'N/A'

        elif self.dataType() == "efficiencyMap":
            nsig = (self.xsection.value * self.dataset.getLumi()).asNumber()
            computer = StatsComputer.forSingleBin(dataset=self.dataset,
                                                  nsig=nsig,
                                                  deltas_rel=self.deltas_rel )

        elif self.dataType() == "combined":
            # Get dictionary with dataset IDs and signal yields
            srNsigDict = {pred.dataset.getID() :
                          (pred.xsection.value*pred.dataset.getLumi()).asNumber()
                          for pred in self.datasetPredictions}

            # Get ordered list of datasets:
            if hasattr(self.dataset.globalInfo, "covariance"):
                datasetList = self.dataset.globalInfo.datasetOrder[:]
                # Get list of signal yields corresponding to the dataset order:
                srNsigs = [srNsigDict[dataID] if dataID in srNsigDict else 0.0
                       for dataID in datasetList]
                # Get computer
                computer = StatsComputer.forMultiBinSL(dataset=self.dataset,
                                                       nsig=srNsigs,
                                                       deltas_rel = self.deltas_rel)

            elif hasattr(self.dataset.globalInfo, "jsonFiles"):
                datasetList = [ds.getID() for ds in self.dataset.origdatasets]
                # Get list of signal yields corresponding to the dataset order:
                srNsigs = [srNsigDict[dataID] if dataID in srNsigDict else 0.0
                       for dataID in datasetList]
                # Get computer
                computer = StatsComputer.forPyhf(dataset=self.dataset,
                                                       nsig=srNsigs,
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
                    element=self.avgElement, txnames=self.txnames, expected=expected
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
                r = (self.totalXsection() / upperLimit).asNumber()
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
        self.cachedObjs[expected]["nll"] = llhdDict["lbsm"]
        self.cachedObjs[expected]["nll_sm"] = llhdDict["lsm"]
        self.cachedObjs[expected]["nllmax"] = llhdDict["lmax"]
        self.cachedObjs[expected]["muhat"] = llhdDict["muhat"]
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
            from smodels.tools.runtime import _deltas_rel_default

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
            selectedTPs, key=lambda tp: (priority[tp.dataType()], tp.getRValue(expected=True))
        )
        # Now get a single TP for each result
        # (the highest ranking analyses come last and are kept in the dict)
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
                ret += tp.xsection.value
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
        if theoryPredictions and isinstance(theoryPredictions, list):

            if not maxCond:
                self._theoryPredictions = theoryPredictions
            else:
                newPredictions = []
                for theoPred in theoryPredictions:
                    mCond = theoPred.getmaxCondition()
                    if mCond != "N/A" and round(mCond / maxCond, 2) > 1.0:
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
        combinedRes = _getCombinedResultFor(dataSetResults, expResult )
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


def _getCombinedResultFor(dataSetResults, expResult ):
    """
    Compute the combined result for all datasets, if covariance
    matrices are available. Return a TheoryPrediction object
    with the signal cross-section summed over all the signal regions
    and the respective upper limit.

    :param datasetPredictions: List of TheoryPrediction objects for each signal region
    :param expResult: ExpResult object corresponding to the experimental result

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
    # Create a theory preidction object for the combined datasets:
    theoryPrediction = TheoryPrediction()
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
        expectedR = (xsec.value / dataset.getSRUpperLimit(expected=True)).asNumber()
        if expectedR > bestExpectedR or (expectedR == bestExpectedR and xsec.value > bestXsec):
            bestExpectedR = expectedR
            bestPred = pred
            bestXsec = xsec.value

    return bestPred


def _getDataSetPredictions(dataset, smsTopList, maxMassDist, deltas_rel=None):
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
        theoryPrediction = TheoryPrediction( deltas_rel)
        theoryPrediction.dataset = dataset
        theoryPrediction.txnames = cluster.txnames
        theoryPrediction.xsection = _evalConstraint(cluster)
        # Skip results with too small (invisible) cross-sections
        if theoryPrediction.xsection.value < 1e-6*fb:
            continue
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
