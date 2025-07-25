#!/usr/bin/env python3

"""
.. module:: statsTools
   :synopsis: a module that contains the class responsible for
              all statistical computations. Designed to
              eventually become simply a frontend for spey

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

__all__ = [ "StatsComputer", "getStatsComputerModule" ]

from typing import Union, Text, Dict, List
from smodels.statistics.exceptions import SModelSStatisticsError as SModelSError
from smodels.base.smodelsLogging import logger
from smodels.base.physicsUnits import fb
from smodels.statistics.simplifiedLikelihoods import LikelihoodComputer, UpperLimitComputer, Data
from smodels.statistics.pyhfInterface import PyhfData, PyhfUpperLimitComputer
from smodels.statistics.basicStats import observed, apriori, aposteriori, NllEvalType
from smodels.statistics.truncatedGaussians import TruncatedGaussians
from smodels.statistics.analysesCombinations import AnaCombLikelihoodComputer
from smodels.experiment.datasetObj import DataSet,CombinedDataSet
from typing import Union, Text

def getStatsComputerModule():
    """ very single convenience function to centralize
    switching between our stats code and spey. """
    from smodels.base import runtime
    if runtime._experimental["spey"]:
        from smodels.statistics.speyTools import SpeyComputer as StatsComputer
        return StatsComputer
    else:
        from smodels.statistics.statsTools import StatsComputer
        return StatsComputer

class StatsComputer:
    __slots__ = [ "nsig", "dataObject", "dataType", "likelihoodComputer", "data",
                  "upperLimitComputer", "deltas_sys", "allowNegativeSignals" ]

    def __init__ ( self, dataObject : Union['DataSet','CombinedDataSet', list], dataType : str,
                   nsig : Union[None,float,List,Dict] = None,
                   deltas_rel : Union[None,float] = None,
                   allowNegativeSignals : bool = False):
        """
         Initialise.

        :param dataObject: a smodels (combined)dataset or a list of theory predictions (for combination of analyses)
        :param nsig: signal yield, either as float or as list
        :param deltas_rel: relative error on signal. currently unused
        :allowNegativeSignals: if True, negative values for the signal (mu) are allowed.
        """

        if dataType not in [ "1bin", "SL", "pyhf", "truncGaussian", "analysesComb"]:
            logger.error ( f"I do not recognize the data type {dataType}" )
            raise SModelSError()

        self.dataObject = dataObject
        self.dataType = dataType
        self.nsig = nsig
        self.deltas_sys = deltas_rel
        if self.deltas_sys is None:
            self.deltas_sys = 0.
        self.allowNegativeSignals = allowNegativeSignals
        self.upperLimitComputer = None
        self.likelihoodComputer = None

    @classmethod
    def forSingleBin(cls, dataset, nsig, deltas_rel):
        """ get a statscomputer for an efficiency map (single bin).

        :param dataset: DataSet object
        :param nsig: Number of signal events for each SR
        :deltas_rel: Relative uncertainty for the signal

        :returns: a StatsComputer
        """
        computer = StatsComputer(dataObject=dataset,
                                 dataType="1bin",
                                 nsig=nsig, deltas_rel=deltas_rel)

        computer.getComputerSingleBin( )

        return computer

    @classmethod
    def forMultiBinSL(cls,dataset, nsig, deltas_rel):
        """ get a statscomputer for simplified likelihood combination.

        :param dataset: CombinedDataSet object
        :param nsig: Number of signal events for each SR
        :deltas_rel: Relative uncertainty for the signal

        :returns: a StatsComputer
        """

        computer = StatsComputer(dataObject=dataset,
                                 dataType="SL",
                                 nsig=nsig, deltas_rel=deltas_rel)

        computer.getComputerMultiBinSL( )

        return computer

    @classmethod
    def forPyhf(cls, dataset, nsig, deltas_rel):
        """ get a statscomputer for pyhf combination.

        :param dataset: CombinedDataSet object
        :param nsig: Number of signal events for each SR
        :deltas_rel: Relative uncertainty for the signal

        :returns: a StatsComputer
        """
        computer = StatsComputer(dataObject=dataset,
                                 dataType="pyhf",
                                 nsig=nsig, deltas_rel=deltas_rel)

        computer.getComputerPyhf( )

        return computer

    @classmethod
    def forTruncatedGaussian(cls,theorypred, corr : float =0.6 ):
        """ get a statscomputer for truncated gaussians
        :param theorypred: TheoryPrediction object
        :param corr: correction factor: \
                ULexp_mod = ULexp / (1. - corr*((ULobs-ULexp)/(ULobs+ULexp))) \
                a factor of corr = 0.6 is proposed.
        :returns: a StatsComputer
        """
        # marked as experimental feature
        if not hasattr(theorypred, "avgElement"):
            logger.error( f"theory prediction {theorypred.analysisId()} has no average element! why??" )
            return None

        eul = theorypred.dataset.getUpperLimitFor(
            element=theorypred.avgElement, txnames=theorypred.txnames, evaluationType=apriori
        )
        if eul is None:
            return None
        eul = eul / theorypred.xsection
        ul = theorypred.dataset.getUpperLimitFor(
            element=theorypred.avgElement, txnames=theorypred.txnames, evaluationType=observed
        ) / theorypred.xsection
        kwargs = { "upperLimitOnMu": float(ul), "expectedUpperLimitOnMu": float(eul),
                   "corr": corr }
        computer = StatsComputer(dataObject=theorypred.dataset,
                                 dataType="truncGaussian",
                                 nsig=0.,
                                 allowNegativeSignals=True)

        computer.getComputerTruncGaussian( **kwargs)

        return computer

    @classmethod
    def forAnalysesComb(cls,theoryPredictions, deltas_rel):
        """ get a statscomputer for combination of analyses
        :param theoryPredictions: list of TheoryPrediction objects
        :param deltas_rel: relative error for the signal
        :returns: a StatsComputer
        """

        # Only allow negative signal if all theory predictions allow for it
        allowNegativeSignals = all([tp.statsComputer.allowNegativeSignals
                                    for tp in theoryPredictions])

        computer = StatsComputer(dataObject=theoryPredictions,
                                 dataType="analysesComb",
                                 nsig=None, deltas_rel=deltas_rel,
                                 allowNegativeSignals=allowNegativeSignals)

        computer.getComputerAnalysesComb( )

        return computer

    def getComputerSingleBin(self):
        """
        Create computer from a single bin

        """

        dataset = self.dataObject
        data = Data( dataset.dataInfo.observedN, dataset.dataInfo.expectedBG,
                     dataset.dataInfo.bgError**2, deltas_rel = self.deltas_sys,
                     nsignal = self.nsig )
        # self.data = data
        self.likelihoodComputer = LikelihoodComputer ( data )
        self.upperLimitComputer = UpperLimitComputer ( self.likelihoodComputer )

    def getComputerMultiBinSL(self):
        """
        Create computer from a multi bin SL result
        """

        dataset = self.dataObject
        cov = dataset.globalInfo.covariance
        if type(cov) != list:
            raise SModelSError( f"covariance field has wrong type: {type(cov)}" )
        if len(cov) < 1:
            raise SModelSError( f"covariance matrix has length {len(cov)}." )

        nobs = [ x.dataInfo.observedN for x in dataset.origdatasets ]
        bg = [ x.dataInfo.expectedBG for x in dataset.origdatasets ]
        third_momenta = [ getattr ( x.dataInfo, "thirdMoment", None ) for x in dataset.origdatasets ]
        # FIXME are we sure thats right, its not the one below?
        # nobs = [ x.dataInfo.observedN for x in dataset._datasets ]
        # bg = [ x.dataInfo.expectedBG for x in dataset._datasets ]
        # third_momenta = [ getattr ( x.dataInfo, "thirdMoment", None ) for x in dataset._datasets ]
        c = third_momenta.count ( None )
        if c > 0:
            if c < len(third_momenta):
                logger.warning ( f"third momenta given for some but not all signal regions in {dataset.globalInfo.id}" )
            third_momenta = None

        data = Data( nobs, bg, cov, third_moment=third_momenta, nsignal = self.nsig,
                     deltas_rel = self.deltas_sys, lumi=dataset.getLumi() )
        self.data = data
        self.likelihoodComputer = LikelihoodComputer ( data )
        self.upperLimitComputer = UpperLimitComputer ( self.likelihoodComputer )

    def getComputerPyhf(self ):
        """
        Create computer for a pyhf result
        """

        globalInfo = self.dataObject.globalInfo
        jsonFiles = [js for js in globalInfo.jsonFiles]
        jsons = globalInfo.jsons.copy()
        # datasets = [ds.getID() for ds in dataset._datasets]
        datasets = [ds.getID() for ds in self.dataObject.origdatasets]
        # Filtering the json files by looking at the available datasets
        for jsName in globalInfo.jsonFiles:
            jsonSRs = []
            for ir,region in enumerate ( globalInfo.jsonFiles[jsName] ): 
                if not isinstance(region,dict):
                    raise SModelSError("The jsonFiles field should contain lists \
                                       of strings or dictionaries \
                                       (%s is not allowed)" %type(region))
                
                globalInfo.jsonFiles[jsName][ir] = region
                if region['type'] == 'SR':
                    jsonSRs.append(region['smodels'])
            if all([ds not in jsonSRs for ds in datasets]):
                # No datasets found for this json combination
                jsIndex = jsonFiles.index(jsName)
                jsonFiles.pop(jsIndex)
                jsons.pop(jsIndex)
                continue
            if not all([SR in datasets for SR in jsonSRs]):
                # Some SRs are missing for this json combination
                logger.error( f"Wrong json definition in globalInfo.jsonFiles for json : {jsName}" )

        jsonDictNames = {}
        for jsName in jsonFiles:
            jsonDictNames.update( { jsName: [ region['smodels'] for region in globalInfo.jsonFiles[jsName] if region is not None and "smodels" in region ] } )
        # jsonRegions = [ [region['smodels'] for region in globalInfo.jsonFiles[jsName]] for jsName in jsonFiles]
        jsonRegions = [ region for regions in jsonDictNames.values() for region in regions ]
        for ds in datasets:
            if not ds in jsonRegions:
                logger.info(f'Region {ds} does not appear in any json file for {globalInfo.id}')
        logger.debug(f"list of datasets: {datasets}")
        logger.debug(f"jsonFiles after filtering: {jsonFiles}")

        # Constructing the list of signals with subsignals matching each json
        nsignals = {}
        for jsName in jsonFiles:
            nsignals.update( { jsName: {} } )
        for name, nsig in self.nsig.items():
            for jsName in nsignals.keys():
                if name in jsonDictNames[jsName]:
                    nsignals[jsName].update( { name: nsig } )

        includeCRs = False
        if hasattr(globalInfo,'includeCRs'):
            includeCRs = globalInfo.includeCRs
        signalUncertainty = None
        if hasattr(globalInfo,"signalUncertainty"):
            signalUncertainty = globalInfo.signalUncertainty

        # Loading the jsonFiles, unless we already have them (because we pickled)
        data = PyhfData(nsignals, jsons, globalInfo.jsonFiles, includeCRs, signalUncertainty)
        if data.errorFlag:
            return None
        self.upperLimitComputer = PyhfUpperLimitComputer(data, lumi=self.dataObject.getLumi() )
        self.likelihoodComputer = self.upperLimitComputer # for pyhf its the same

    def getComputerTruncGaussian ( self, **kwargs ):
        """
        Create computer for truncated gaussians
        """
        computer = TruncatedGaussians ( **kwargs )
        self.data = None
        self.likelihoodComputer = computer
        self.upperLimitComputer = computer

    def getComputerAnalysesComb(self):
        """
        Create computer from a single bin
        :param nsig: signal yields.
        """

        self.upperLimitComputer = AnaCombLikelihoodComputer(theoryPredictions=self.dataObject,
                                                            deltas_rel=self.deltas_sys)
        self.likelihoodComputer = self.upperLimitComputer # for analyses combination its the same

    def get_five_values ( self, evaluationType : NllEvalType,
                      return_nll : bool = False,
                      check_for_maxima : bool = False )-> Dict:
        """
        Return the Five Values: l(bsm), l(sm), muhat, l(muhat), sigma(mu_hat)

        :param check_for_maxima: if true, then check lmax against l(sm) and l(bsm)
                                 correct, if necessary
        """
        ret = self.maximize_likelihood ( evaluationType = evaluationType, return_nll = return_nll  )
        if ret is None:
            return {}
        lmax = ret['lmax']

        lbsm = self.likelihood ( poi_test = 1., evaluationType=evaluationType, return_nll = return_nll )
        ret["lbsm"] = lbsm
        lsm = self.likelihood ( poi_test = 0., evaluationType=evaluationType, return_nll = return_nll )
        ret["lsm"] = lsm
        if check_for_maxima:
            if return_nll:
                if lsm < lmax: ## if return_nll is off, its the other way
                    muhat = ret["muhat"]
                    logger.debug(f"lsm={lsm:.2g} > lmax({muhat:.2g})={lmax:.2g}: will correct")
                    ret["lmax"] = lsm
                    ret["muhat"] = 0.0
                if lbsm < lmax:
                    muhat = ret["muhat"]
                    logger.debug(f"lbsm={lbsm:.2g} > lmax({muhat:.2g})={lmax:.2g}: will correct")
                    ret["lmax"] = lbsm
                    ret["muhat"] = 1.0
            else:
                if lsm > lmax:
                    muhat = ret["muhat"]
                    logger.debug(f"lsm={lsm:.2g} > lmax({muhat:.2g})={lmax:.2g}: will correct")
                    ret["lmax"] = lsm
                    ret["muhat"] = 0.0
                if lbsm > lmax:
                    muhat = ret["muhat"]
                    logger.debug(f"lbsm={lbsm:.2g} > lmax({muhat:.2g})={lmax:.2g}: will correct")
                    ret["lmax"] = lbsm
                    ret["muhat"] = 1.0

        return ret

    def likelihood ( self, poi_test : float, evaluationType : NllEvalType,
                  return_nll : bool, asimov : Union[None,float] = None ) -> float:
        """ simple frontend to individual computers """
        self.transform ( evaluationType )
        kwargs = { "evaluationType": evaluationType, "asimov": asimov }
        if self.dataType == "pyhf":
            if not "workspace_index" in kwargs:
                index = self.likelihoodComputer.getBestCombinationIndex()
                kwargs["workspace_index"] = index
            ret = self.likelihoodComputer.likelihood (
                    poi_test, return_nll = return_nll, **kwargs )
            return ret
        ret = self.likelihoodComputer.likelihood ( poi_test,
                                            return_nll = return_nll, **kwargs)
        return ret

    def CLs ( self, poi_test : float = 1., evaluationType : NllEvalType=observed ) -> Union[float,None]:
        """ compute CLs value for a given value of the poi """
        # self.transform ( evaluationType )
        if hasattr ( self.upperLimitComputer, "CLs" ):
            return self.upperLimitComputer.CLs ( poi_test, evaluationType )
        return None

    def transform ( self, evaluationType ):
        """ SL only. transform the data to evaluationType or observed """
        if self.dataType in [ "pyhf", "truncGaussian", "analysesComb" ]:
            return
        self.likelihoodComputer.transform ( evaluationType )

    def restore ( self, evaluationType ):
        """ SL only. transform the data to evaluationType or observed """
        if self.dataType in [ "pyhf", "truncGaussian", "analysesComb" ]:
            return
        if evaluationType != observed:
            return
        self.likelihoodComputer.model = self.likelihoodComputer.origModel

    def maximize_likelihood ( self, evaluationType : NllEvalType,
           return_nll : bool = False ) -> dict:
        """ simple frontend to the individual computers, later spey
        :param return_nll: if True, return negative log likelihood
        :returns: Dictionary of llhd (llhd at mu_hat), \
                  muhat, sigma_mu (sigma of mu_hat), \
                  optionally also theta_hat
        """
        self.transform ( evaluationType )
        kwargs = { }
        if self.dataType == "pyhf":
            if not "workspace_index" in kwargs:
                index = self.likelihoodComputer.getBestCombinationIndex()
                kwargs["workspace_index"] = index
                if evaluationType != observed:
                    kwargs["evaluationType"] = evaluationType
        elif self.dataType == "truncGaussian":
            kwargs["evaluationType"]=evaluationType
        elif self.dataType == "analysesComb":
            kwargs["evaluationType"]=evaluationType

        ret = self.likelihoodComputer.lmax ( return_nll = return_nll,
                                   allowNegativeSignals = self.allowNegativeSignals,
                                   **kwargs )
        return ret

    def poi_upper_limit ( self, evaluationType : NllEvalType,
           limit_on_xsec : bool = False ) -> float:
        """
        Simple frontend to the upperlimit computers, later to spey.poi_upper_limit

        :param limit_on_xsec: if True, then return the limit on the
                              cross section
        """
        if self.dataType == "pyhf":
            if all([s == 0 for s in self.nsig]):
                logger.warning("All signals are empty")
                return None
            index = self.likelihoodComputer.getBestCombinationIndex()
            if limit_on_xsec:
                ret = self.upperLimitComputer.getUpperLimitOnSigmaTimesEff(
                       evaluationType = evaluationType, workspace_index = index )
            else:
                ret = self.upperLimitComputer.getUpperLimitOnMu(
                       evaluationType = evaluationType, workspace_index = index )
        elif self.dataType in ["SL", "1bin", "truncGaussian"]:
            self.upperLimitComputer.likelihoodComputer.model = self.data
            if limit_on_xsec:
                ret = self.upperLimitComputer.getUpperLimitOnSigmaTimesEff(
                       evaluationType = evaluationType )
            else:
                ret = self.upperLimitComputer.getUpperLimitOnMu(
                       evaluationType = evaluationType )
        elif self.dataType in ["analysesComb"]:
            if limit_on_xsec:
                ret = self.upperLimitComputer.getUpperLimitOnSigmaTimesEff(evaluationType = evaluationType,
                                                                           allowNegativeSignals=self.allowNegativeSignals )
            else:
                ret = self.upperLimitComputer.getUpperLimitOnMu(evaluationType = evaluationType,
                                                                allowNegativeSignals=self.allowNegativeSignals )

        return ret

class SimpleStatsDataSet:
    """ a very simple data class that can replace a smodels.dataset,
    for 1d SL data only. used for testing and in dataPreparation """
    class SimpleInfo:
        def __init__ ( self, observedN : float, evaluationTypeBG : float,
                       bgError : float ):
            self.observedN = observedN
            self.expectedBG = expectedBG
            self.bgError = bgError

    class GlobalInfo:
        def __init__ ( self, lumi ):
            self.id = "SimpleStatsDataSet"
            self.lumi = lumi

    def __init__ ( self, observedN : float, evaluationTypeBG : float,
                   bgError : float, lumi : fb ):
        """ initialise the dataset with all relevant stats """
        self.dataInfo = self.SimpleInfo ( observedN, evaluationTypeBG, bgError )
        self.globalInfo = self.GlobalInfo( lumi )

    def getLumi ( self ):
        return self.globalInfo.lumi

    def getType ( self ):
        return "efficiencyMap"

if __name__ == "__main__":
    # nobs,bg,bgerr,lumi = 3., 4.1, 0.6533758489567854, 35.9/fb
    # nobs,bg,bgerr,lumi = 0, 1., 0.2, 35.9/fb
    # nobs,bg,bgerr,lumi = 0, 0.001, 0.01, 35.9/fb
    nobs,bg,bgerr,lumi = 3905,3658.3,238.767, 35.9/fb
    dataset = SimpleStatsDataSet ( nobs, bg, bgerr, lumi )
    computer = StatsComputer ( dataset, 1. )
    ul = computer.poi_upper_limit ( evaluationType = observed, limit_on_xsec = True )
    print ( "ul", ul )
    ule = computer.poi_upper_limit ( evaluationType = apriori, limit_on_xsec = True )
    print ( "ule", ule )
