"""
.. module:: theoryPrediction
   :synopsis: Provides a class to encapsulate the results of the computation of
              reference cross sections and related functions.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

from smodels.base.physicsUnits import TeV, fb, UnitXSec
from smodels.experiment.datasetObj import CombinedDataSet
from smodels.experiment.databaseObj import Database
from smodels.matching.exceptions import SModelSMatcherError as SModelSError
from smodels.statistics.basicStats import observed, apriori, NllEvalType, \
         aposteriori
from smodels.matching import clusterTools
from smodels.base.smodelsLogging import logger
from smodels.tools.caching import roundCache,lru_cache
from typing import Union, Dict, Callable, Optional, List
import numpy as np

# number of digits for rounding the mu argument when computing likelihoods
mu_digits = 8

__all__ = [ "TheoryPrediction", "theoryPredictionsFor",
            "TheoryPredictionsCombiner" ]


class TheoryPrediction(object):
    """
    An instance of this class represents the results of the theory prediction
    for an analysis.
    """

    def __init__(self, deltas_rel : Union[None,float] = None):
        """
        Initialize the theory prediction object. deltas_rel is meant to be a
        constant.

        :param deltas_rel: relative uncertainty in signal (float).
                           Default value is 20%.
        """
        self.analysis = None
        self.xsection = None
        self.conditions = None
        self.mass = None
        self.totalwidth = None
        self.dataset = None
        self.txnames = []
        self.smsList = []
        if deltas_rel is None:
            from smodels.base.runtime import _deltas_rel_default
            deltas_rel = _deltas_rel_default
        self.deltas_rel = deltas_rel
        self._statsComputer = None

    def __str__(self):
        ret = f"{self.analysisId()}:{self.totalXsection()}"
        return ret

    def dataId(self) -> Union[None,str]:
        """
        Return ID of dataset
        """
        if self.dataset is None:
            return None
        
        return self.dataset.getID()

    def analysisId(self) -> Union[None,str]:
        """
        Return experimental analysis ID
        """

        if self.dataset is None:
            return None
        
        return self.dataset.globalInfo.id

    def dataType(self, short : bool = False ) -> Union[None,str]:
        """
        Return the type of dataset
        :param: short, if True, return abbreviation (ul,em,comb)

        :returns: data type, as a string, e.g. "efficiencyMap"
        """
        if self.dataset is None:
            return None

        if short:
            t = self.dataset.getType()
            D = {"upperLimit": "ul", "efficiencyMap": "em",
                 "combined": "comb"}
            if t in D.keys():
                return D[t]
            return "??"

        return self.dataset.getType()
    
    def type(self) -> str:
        """
        Return the type of theory prediction (single result)
        """

        return str(type(self).__name__)


    def computeXSection(self):
        if type(self.xsection) != type(None):
            return

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
        self.computeXSection()
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

    def getTxNamesWeights(self,sort : bool = True) -> dict:
        """
        Returns a dictionary with the txname objects as keys
        and their total weights as values.
        :param sort: If True, the dictionary is sorted according
        to the largest weights.
        """

        txnamesWeightsDict = {}
        for sms in self.smsList:
            if sms.txname not in txnamesWeightsDict:
                txnamesWeightsDict[sms.txname] = sms.weight.asNumber(fb)
            else:
                txnamesWeightsDict[sms.txname] += sms.weight.asNumber(fb)

        if sort:
            txnamesWeightsDict = dict(sorted(txnamesWeightsDict.items(),
                                        key=lambda x: (x[1],str(x[0])), reverse=True))

        return txnamesWeightsDict

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
        from smodels.statistics.statsTools import StatsComputer

        statsComputer = StatsComputer.forTheoryPrediction(self)

        if statsComputer is None:
            self._statsComputer = "N/A"
        else:
            self._statsComputer = statsComputer

    @lru_cache
    def getUpperLimit( self, evaluationType : NllEvalType = observed,
                       nSigma : int = 0, **kwargs ) -> UnitXSec:
        """
        Get the upper limit on sigma*eff.
        For UL-type results, use the UL map. For EM-Type returns
        the corresponding dataset (signal region) upper limit.
        For combined results, returns the upper limit on the
        total sigma*eff (for all signal regions/datasets).
        :param evaluationType: one of: observed, apriori, aposteriori
        :param nSigma: the upper limit for central value (0),
        + 1 sigma, - 1 sigma, etc.
        For error bands.
        :return: upper limit (Unum object)
        """
        
        if self.dataType() == "efficiencyMap":
            ul = self.dataset.getSRUpperLimit(evaluationType=evaluationType,
                    nSigma = nSigma )
        if self.dataType() == "upperLimit":
            ul = self.dataset.getUpperLimitFor(
                sms=self.avgSMS, txnames=self.txnames,
                evaluationType=evaluationType,
                nSigma = nSigma
            )
        if self.dataType() == "combined":
            ul = self.statsComputer.getUpperLimit(
                evaluationType = evaluationType,
                nSigma = nSigma,
                limit_on_xsec = True, **kwargs )
        return ul

    def getUpperLimitOnMu(self, evaluationType : NllEvalType = observed,
            nSigma : int = 0, **kwargs ) -> Union[None,float]:
        """
        Get upper limit on signal strength multiplier, using the
        theory prediction value and the corresponding upper limit
        (i.e. mu_UL = upper limit/theory xsec)
        :param evaluationType: one of: observed, apriori, aposteriori
        :param nSigma: the upper limit for central value (0),
        + 1 sigma, - 1 sigma, etc.
        For error bands.
        :returns: upper limit on signal strength multiplier mu
        """
        if "expected" in kwargs:
            import warnings
            warnings.warn ( "flag 'expected' in theoryPrediction.getRValue() renamed to evaluationType, please adapt!", DeprecationWarning, stacklevel=2 )
            evaluationType = kwargs["expected"]
            kwargs.pop ( "expected" )
        for k,v in kwargs.items():
            if k not in [ "pmSigma" ]:
                logger.error ( f"unknown argument {k} in theoryPrediction.getUpperLimitOnMu()" )
        upperLimit = self.getUpperLimit( evaluationType=evaluationType,
                                         nSigma=nSigma, **kwargs )
        xsec = self.totalXsection()
        if xsec is None or upperLimit is None:
            return None

        muUL = (upperLimit/xsec).asNumber()

        return muUL

    @lru_cache
    def getRValue( self, evaluationType : NllEvalType = observed,
                   nSigma : int = 0, **kwargs ) -> Union[float,None]:
        """
        Get the r value = theory prediction / experimental upper limit

        :param evaluationType: one of: observed, apriori, aposteriori
        """
        if "expected" in kwargs:
            import warnings
            warnings.warn ( "flag 'expected' in theoryPrediction.getRValue() renamed to evaluationType, please adapt!", DeprecationWarning, stacklevel=2 )
            evaluationType = kwargs["expected"]
            kwargs.pop ( "expected" )
        for k,v in kwargs.items():
            if k not in [ "pmSigma" ]:
                logger.error ( f"unknown argument {k} in theoryPrediction.getRValue()" )
        upperLimit = self.getUpperLimit(evaluationType, nSigma = nSigma, **kwargs )
        if upperLimit is None or upperLimit.asNumber(fb) == 0.0:
            r = None
            return r
        else:
            r = (self.totalXsection()/upperLimit).asNumber()
            return r

    def whenDefined( function : Callable ) -> Callable:
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
    def nllsm( self, evaluationType : NllEvalType = observed,
             return_nll : bool = False ) -> float:
        """likelihood at SM point, same as .def likelihood( ( mu = 0. )"""
        nllDict = self.computeStatistics(evaluationType)
        return nllDict["nllsm"]

    @whenDefined
    def nll_min( self, evaluationType : NllEvalType = observed,
                 return_dict : bool = False,
                 **kwargs ) -> Union[float,dict]:
        """ nll at mu_hat
        :param return_dict: if true, return also mu_hat and sigma_mu
        """
        nllDict = self.computeStatistics(evaluationType)
        if len(kwargs)>0:
            nllDict = self.statsComputer.nll_min (
                    evaluationType = evaluationType, **kwargs )
        if return_dict:
            return nllDict
        return nllDict["nll_min"]

    @whenDefined
    @roundCache(argname='mu',argpos=1,digits=mu_digits)
    def CLs( self, mu : float = 1., evaluationType : NllEvalType = observed,
             **kwargs ) -> Union[float,None]:
        """ obtain the CLs value of the combination for a given poi value "mu" """
        cls = self.statsComputer.CLs ( poi_test = mu,
                evaluationType = evaluationType, **kwargs )
        return cls

    @whenDefined
    def sigma_mu(self, evaluationType : NllEvalType = observed ):
        """sigma_mu of mu_hat"""
        llhDict = self.computeStatistics(evaluationType)
        return llhDict["sigma_mu"]

    @whenDefined
    def muhat(self, evaluationType : NllEvalType = observed ):
        """position of maximum likelihood"""
        llhDict = self.computeStatistics(evaluationType)
        return llhDict["muhat"]

    @whenDefined
    @roundCache(argname='mu',argpos=1,digits=mu_digits)
    def nll(self, mu=1.0, evaluationType : NllEvalType = observed,
            asimov : Optional[int] = None, **kwargs ) -> float:
        """
        get the nll for a signal strength modifier mu
        :param evaluationType: one of: observed, apriori, aposteriori
        :asimov: get for asimov data with this mu, or None
        """
        if asimov != None and abs(asimov)>1e-8 and abs(asimov-1)>1e-8:
            raise SModelSError (
               "currently we only handle asimov data for 0. or 1." )
        if not "writeYields" in kwargs or kwargs["writeYields"]==True:
            if abs(mu-5)<1e-5 and evaluationType == aposteriori:
                from smodels.statistics.nnInterface import writeOutYields
                writeOutYields ( self )
        if "writeYields" in kwargs:
            kwargs.pop ( "writeYields" )
        if "expected" in kwargs:
            import warnings
            warnings.warn ( "flag 'expected' in theoryPrediction.getRValue() renamed to evaluationType, please adapt!", DeprecationWarning, stacklevel=2 )
            evaluationType = kwargs["expected"]
        for k,v in kwargs.items():
            if k not in [ "pmSigma" ]:
                logger.error ( f"unknown argument {k} in theoryPrediction.likelihood()" )

        # for truncated gaussians the fits only work with negative signals!
        nll = self.statsComputer.nll(poi_test = mu,
                evaluationType = evaluationType, asimov = asimov, **kwargs )
        return nll

    @whenDefined
    def likelihood(self, mu : float = 1.0, evaluationType : NllEvalType = observed,
            return_nll : bool =False, asimov : Optional[int] = None,
            **kwargs ) -> float:
        """
        get the likelihood for a signal strength modifier mu.
        this is a method to prepare for a transition to dealing with nlls only

        :param evaluationType: one of: observed, apriori, aposteriori
        """
        mnll = self.nll ( mu=mu, evaluationType=evaluationType,
                         asimov=asimov, **kwargs )
        return self.nllToLikelihood ( mnll, return_nll )

    def nllToLikelihood ( self, nll : Union[None,float], return_nll : bool ):
        """ if not return_nll, then compute likelihood from nll """
        if return_nll:
            return nll
        return np.exp ( - nll ) if nll is not None else None

    @whenDefined
    @lru_cache
    def computeStatistics(self, evaluationType : NllEvalType = observed ):
        """
        Compute the likelihoods, and upper limit for this theory prediction.
        The resulting values are stored as the likelihood, lmax, and lsm
        attributes.
        :param evaluationType: one of: observed, apriori, aposteriori
        """
        # Compute likelihoods and related parameters:
        llhdDict = self.statsComputer.get_five_values(
                evaluationType = evaluationType, return_nll = True )
        return llhdDict


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
        ret = f"{len(self._theoryPredictions)} predictions: "
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


class TheoryPredictionsCombiner(TheoryPrediction):
    """
    Facility used to combine theory predictions from different analyes.
    If a list with a single TheoryPrediction is given, return the TheoryPrediction
    object.
    """

    def __new__(cls,theoryPredictions: list,
                slhafile : Union[None,str] = None,
                deltas_rel : Union[float,None] = None ):
        """
        If called with a list containing a single TheoryPrediction, return the TheoryPrediction object.
        Otherwise, create a TheoryPredictionsCombiner object.
        """

        if len(theoryPredictions) == 1:
            return theoryPredictions[0]
        else:
            tpCombiner = super(TheoryPredictionsCombiner, cls).__new__(cls)
            return tpCombiner

    def __init__( self, theoryPredictions: list,
                  slhafile : Union[None,str]=None,
                  deltas_rel : Union[float,None]=None ):
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
        self._statsComputer = None

    @classmethod
    def selectResultsFrom(cls, theoryPredictions : Union[List[TheoryPrediction],TheoryPredictionList], 
                          anaIDs : List[str] ):
        """
        Select the results from theoryPrediction list which match one
        of the IDs in anaIDs. If there are multiple predictions for the
        same ID for which a likelihood is available, it gives priority
        to the ones with largest evaluationType r-values.

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
                f"No likelihood available for {anaID}. This analysis will not be used in analysis combination." )
        # If no results are available, return None
        if len(selectedTPs) == 0:
            return None

        # Define a hierarchy for the results:
        priority = {"combined": 2, "efficiencyMap": 1, "upperLimit": 0}
        # Now sort by highest priority and then by highest evaluationType r-value:
        selectedTPs = sorted(
            selectedTPs, key=lambda tp: (priority[tp.dataType()],
                                 tp.getRValue(evaluationType=apriori) is not None,
                                 tp.getRValue(evaluationType=apriori))
        )
        # Now get a single TP for each result
        # (the highest ranking analyses with r != None come last
        # and are kept in the dict)
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
        Return a string with the IDs of all the experimental results
        used in the combination.
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

    def getTxNamesWeights(self, sort : bool = True ) -> dict:
        """
        Returns a dictionary with the txname objects as keys
        and their total weights as values.
        :param sort: If True, the dictionary is sorted according to
        the largest weights.
        """

        txnamesWeightsDict = {}
        for theoPred in self.theoryPredictions:
            txW = theoPred.getTxNamesWeights(sort=False)
            for tx,w in txW.items():
                if tx not in txnamesWeightsDict:
                    txnamesWeightsDict[tx] = 0.0
                txnamesWeightsDict[tx] += w

        if sort:
            txnamesWeightsDict = dict(sorted(txnamesWeightsDict.items(),
                                      key=lambda x: (x[1],x[0]), reverse=True))

        return txnamesWeightsDict

    def getLlhds(self, **kwargs ) -> dict:
        """
        Facility to access the likelihoods for the individual analyses and
        the combined likelihood.
        Returns a dictionary with the analysis IDs as keys and the likelihood
        values as values.  Mostly used for plotting the likelihoods.

        :param muvals: List with values for the signal strenth for which
        the likelihoods must be evaluated.
        :param evaluationType: returns the observed/priori expected/posteriori expected likelihood values.
        :param normalize: If True normalizes the likelihood by its integral
        over muvals.
        """
        return self.statsComputer.getLlhds( **kwargs )

    def describe(self):
        """returns a string containing a list of all analysisId and dataIds"""
        ids = []
        for tp in self.theoryPredictions:
            ids.append(f"{tp.analysisId()}:{tp.dataId()}")
        return f"SRs: {', '.join(ids)}"


def theoryPredictionsFor(database : Database, smsTopDict : Dict,
        maxMassDist : float = 0.2, useBestDataset : bool = True,
        combinedResults : bool = True,
        deltas_rel : Union[None,float] = None) -> list:
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
        errorMsg += f" database with the selected experimental results and not {type(database)}"
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
                if bestRes is not None:
                    expResults.append(bestRes) # Best result = combination if available

        for theoPred in expResults:
            theoPred.expResult = expResult
            theoPred.deltas_rel = deltas_rel
            srType = 'SR' # By default assume it is a SR
            # Check if SR corresponds to a CR or VR, in which case we do not report the upper limit
            if hasattr(expResult.globalInfo, "srMappings"):
                assert type ( expResult.globalInfo.srMappings) == dict, f"{expResult.globalInfo.id} srMappings are not a dict"
                srType = expResult.globalInfo.srMappings.get(theoPred.dataId(),None)
                if srType is not None:
                    srType = srType["type"]
                
            if srType == "SR":
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

def _isDatasetInCombination ( dataset, expResult ) -> Union[None,bool]:
    """ is a given dataset mentioned in the combination?
    we are allowing datasets in an expResult that is not mentioned
    in the combination of that result.
    :returns: true if its in, false if it is not in, None if there is no combination
    """
    assert hasattr ( dataset, "dataInfo" ), \
        "why does the dataset here not have a dataInfo?"
    dataId = dataset.getID()
    if not hasattr ( expResult.globalInfo, "srMappings" ):
        return False
    if not hasattr ( expResult.globalInfo, "srSets" ):
        raise SModelSError ( f"{expResult.globalInfo.id} has srMappings but no srSets" )
    # if dataId not in expResult.globalInfo.srMappings:  ## [!AL!] dataId refers to the smodels name, no? So it should not be in a key of srMappings, which corresponds to the label
        # return False
    
    for srSetName in expResult.globalInfo.statModels.keys():
        regionList = expResult.globalInfo.srSets[srSetName] # List of labels for the signal regions
        for region in regionList:
            regionDict = expResult.globalInfo.srMappings[region] # dictionary mapping the label to the smodels, pyhf,... names
            if dataId == regionDict['smodels']:
                return True
    return False

def _getCombinedResultFor(dataSetResults, expResult):
    """
    Compute the combined result for all datasets, if a statistical model is]
    available. Return a TheoryPrediction object
    with the signal cross-section summed over all the signal regions
    and the respective upper limit.

    :param datasetPredictions: List of TheoryPrediction objects
    for each signal region
    :param expResult: ExpResult object corresponding to the experimental result

    :return: TheoryPrediction object
    """

    if len(dataSetResults) == 1:
        return dataSetResults[0]
    elif not expResult.hasStatsModel():
        return None
    
    # Don't give combined result if all regions are CRs
    srTypeDict = {}
    globalInfo = expResult.globalInfo
    for predList in dataSetResults:
        for tpred in predList:
            dataId = tpred.dataId()
            srTypeDict[dataId] = 'SR'
            if hasattr ( globalInfo, "srMappings" ): ## [!AL!] dataId refers to the smodels name, no? So it should not be in a key of srMappings, which corresponds to the label
                    for regionDict in globalInfo.srMappings.values():
                        if dataId == regionDict['smodels']:
                            srTypeDict[dataId] = regionDict['type']
    
    if all(srType != 'SR' for srType in srTypeDict.values()):
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
        if not _isDatasetInCombination ( pred.dataset, expResult ):
            continue
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
    Returns the best result according to the evaluationType upper limit.
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
    # select the best one according to the evaluationType upper limit:
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
        if hasattr(globalInfo,"srMappings"):
            regionType = "SR" # Assume it is a SR by default
            for regionDict in globalInfo.srMappings.values():
                if dataset.dataInfo.dataId == regionDict["smodels"]:
                    regionType = regionDict["type"]
            
            if regionType != "SR":
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
        expectedR = (xsec/dataset.getSRUpperLimit(evaluationType=apriori)).asNumber()
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
        logger.error(f"Sqrt(s) defined with wrong units for {ID}")
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
        logger.warning(f"Unkown data type: {dataset.getType()}. Data will be ignored.")

    return clusters
