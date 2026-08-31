#!/usr/bin/env python3

"""
.. module:: statsTools
   :synopsis: a module that contains the class responsible for
              all statistical computations. Designed to
              eventually become simply a frontend for spey

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

__all__ = [ "StatsComputer", "getCompRetrieverModule" ]

from typing import Union, Dict, Optional
from smodels.statistics.exceptions import SModelSStatisticsError as SModelSError
from smodels.base.smodelsLogging import logger
from smodels.base.physicsUnits import fb, UnitLumi
from smodels.statistics.simplifiedLikelihoods import SLLikelihoodComputer, SLUpperLimitComputer, SLData
from smodels.statistics.pyhfInterface import PyhfData, PyhfUpperLimitComputer
from smodels.statistics.nnInterface import NNData, NNUpperLimitComputer
from smodels.statistics.basicStats import observed, apriori, NllEvalType, exponentiateNLL
from smodels.statistics.truncatedGaussians import TruncatedGaussians
from smodels.statistics.analysesCombinations import AnaCombLikelihoodComputer
from smodels.base.physicsUnits import UnitXSec
from smodels.tools.caching import lru_cache
from smodels.experiment.datasetObj import DataSet, CombinedDataSet

def getCompRetrieverModule() -> type:
    """ very single convenience function to centralize
    switching between our stats code and spey. """
    from smodels.base import runtime
    if runtime._experimental["spey"]:
        from smodels.statistics.speyTools import SpeyRetriever as CompRetriever
        return CompRetriever
    else:
        from smodels.statistics.statsTools import CompRetriever
        return CompRetriever

class CompRetriever:
    """ simple class that retrieves and constructs the sub computers required by StatsComputer """

    @classmethod
    def forMultiBinSL(cls,regionSet : str, dataset: CombinedDataSet, nsigDict: dict,
            deltas_rel : Optional[float] = 0.0  ) -> SLUpperLimitComputer:
        """ get a subcomputer for simplified likelihood sr-combination.

        :param dataset: CombinedDataSet object
        :param regionSet: the regionSet defining the SRs for which to construct the computer
        :param nsigDict: Dictionary of signal events for each SR
        :param deltas_rel: Relative uncertainty for the signal

        :returns: a subcomputer
        """

        if regionSet not in dataset.globalInfo.statModels:
            raise SModelSError( f"{regionSet} not in statModels in {dataset.globalInfo.id}" )
        covs = dataset.globalInfo.cachedModels
        type_n_models = dataset.globalInfo.statModels[regionSet]        
        mtype,covname = type_n_models[0] # get first statistical model
        if mtype != "sl":
            raise SModelSError(f"expected sl but got {mtype} for type of stats model in {dataset.globalInfo.id}")
        
        srList = dataset.globalInfo.regionSets[regionSet]
        cov = covs[covname]
        if not isinstance(cov, list):
            logger.error(f"covariance field has wrong type: {type(cov)} in {dataset.globalInfo.id}")
            raise SModelSError(f"covariance field has wrong type: {type(cov)} in {dataset.globalInfo.id}")
        if len(cov) == 0:
            logger.error(f"covariance matrix has length {len(cov)} but regionSet {regionSet} has {len(srList)} signal regions in {dataset.globalInfo.id}")
            raise SModelSError(f"covariance matrix has length {len(cov)} but regionSet {regionSet} has {len(srList)} signal regions in {dataset.globalInfo.id}")
        
        # Collect relevant data:
        nobs = []
        bg = []
        nsig = []
        third_momenta = []
        for sr in srList:
            ds = dataset.getDataSet(sr)         
            if ds is None:
                raise SModelSError(f"SR {sr} defined in regionSet {regionSet} not found in dataset {dataset.globalInfo.id}")
            nobs.append(ds.dataInfo.observedN)
            bg.append(ds.dataInfo.expectedBG)
            nsig.append(nsigDict.get(sr, 0.0))
            third_momenta.append(getattr(ds.dataInfo, "thirdMoment", None ))

        c = third_momenta.count ( None )
        if c > 0:
            if c < len(third_momenta):
                logger.warning ( f"third momenta given for some but not all signal regions in {dataset.globalInfo.id}" )
            third_momenta = None

        data = SLData( nobs, bg, cov, third_moment=third_momenta,
                     nsignal = nsig,
                     deltas_rel = deltas_rel, lumi=dataset.getLumi(),
                     name = covname )
        likelihoodComputer = SLLikelihoodComputer ( data )
        computer = SLUpperLimitComputer ( likelihoodComputer )
        computer.dataType = "SL"
        computer.allowNegativeSignals = False
        return computer

    @classmethod
    def forSingleBin( cls, regionSet : str, dataset : DataSet, nsigDict : dict,
                      deltas_rel : float = 0.2,
                      lumi : Optional[UnitLumi]=None ) ->  SLUpperLimitComputer:
        """ get a sub computer for an efficiency map (single bin).

        :param dataset: DataSet object
        :param regionSet: the regionSet defining the SRs for which to construct the computer
        :param nsigDict: Dictionary of signal events for each SR
        :param deltas_rel: Relative uncertainty for the signal

        :returns: a sub computer
        """
        if regionSet != dataset.getID():
            logger.error ( f"regionSet {regionSet} does not match dataset id {dataset.getID()}" )
            raise SModelSError ( f"regionSet {regionSet} does not match dataset id {dataset.getID()}" )
        data = SLData( dataset.dataInfo.observedN, dataset.dataInfo.expectedBG,
                     dataset.dataInfo.bgError**2, deltas_rel = deltas_rel,
                     nsignal = list(nsigDict.values()), lumi = lumi,
                     name = regionSet )
        likelihoodComputer = SLLikelihoodComputer ( data )
        computer = SLUpperLimitComputer ( likelihoodComputer )
        computer.dataType = "1bin"
        computer.allowNegativeSignals = False
        return computer

    @classmethod
    def forNNs(cls, regionSet : str, dataset: CombinedDataSet, nsigDict: dict ) -> NNUpperLimitComputer:
        """ get a sub computer for an NN combination.

        :param regionSet: the regionSet defining the SRs for which to construct the computer
        :param dataset: CombinedDataSet object
        :param nsigDict: Dictionary of signal events for each SR

        :returns: a sub computer
        """
        globalInfo = dataset.globalInfo
        labelToONNX = {}
        regionMappings = globalInfo.regionMappings

        for sr in globalInfo.regionSets[regionSet]:
            if sr not in regionMappings:
                logger.error ( f"SR {sr} defined in regionSet {regionSet} not found in regionMappings for dataset {dataset.globalInfo.id}" )
                raise SModelSError ( f"SR {sr} defined in regionSet {regionSet} not found in regionMappings for dataset {dataset.globalInfo.id}" )
            labelToONNX[sr] = regionMappings[sr]["onnx"]

        if regionSet not in globalInfo.statModels:
            logger.error ( f"regionSet {regionSet} not found in statModels for dataset {dataset.globalInfo.id}" )
            raise SModelSError ( f"regionSet {regionSet} not found in statModels for dataset {dataset.globalInfo.id}" )

        modelList = globalInfo.statModels[regionSet]
        if len(modelList) == 0:
            logger.error ( f"no model defined for regionSet {regionSet} in dataset {dataset.globalInfo.id}" )
            raise SModelSError ( f"no model defined for regionSet {regionSet} in dataset {dataset.globalInfo.id}" )

        model_type, model_filename = modelList[0] # Always use first model
        if model_type != "onnx":
            logger.error ( f"model type {model_type} for regionSet {regionSet} in dataset {dataset.globalInfo.id} is not 'onnx'" )
            raise SModelSError ( f"model type {model_type} for regionSet {regionSet} in dataset {dataset.globalInfo.id} is not 'onnx'" )

        # Get dictionary for signal yields using the ONNX labels
        f_signals = {onnx_sr : nsigDict.get(label,0.0) for label,onnx_sr in labelToONNX.items()}

        data = NNData( f_signals, dataset )
        upperLimitComputer = NNUpperLimitComputer(data, lumi=dataset.getLumi(),
                                                  onnxfilename = model_filename )

        return upperLimitComputer

    @classmethod
    def forPyhf( cls, regionSet : str, dataset : CombinedDataSet,
                 nsigDict : dict ) -> PyhfUpperLimitComputer:
        """ get a sub computer for pyhf combination.

        :param regionSet: name of signal region set
        :param dataset: CombinedDataSet object
        :param nsigDict: Dictionary with signal yields for each SR

        :returns: a sub computer
        """

        globalInfo = dataset.globalInfo
        regionMappings = globalInfo.regionMappings
        labelToPyhf = {}
        for sr_label in globalInfo.regionSets[regionSet]:
            if sr_label not in regionMappings:
                logger.error ( f"SR {sr_label} defined in regionSet {regionSet} not found in regionMappings for dataset {dataset.globalInfo.id}" )
                raise SModelSError ( f"SR {sr_label} defined in regionSet {regionSet} not found in regionMappings for dataset {dataset.globalInfo.id}" )
            if regionMappings[sr_label]["smodels"] is not None:
                labelToPyhf[sr_label] = regionMappings[sr_label]["pyhf"]

        if regionSet not in globalInfo.statModels:
            logger.error ( f"regionSet {regionSet} not found in statModels for dataset {dataset.globalInfo.id}" )
            raise SModelSError ( f"regionSet {regionSet} not found in statModels for dataset {dataset.globalInfo.id}" )

        modelList = globalInfo.statModels[regionSet]
        if len(modelList) == 0:
            logger.error ( f"no model defined for regionSet {regionSet} in dataset {dataset.globalInfo.id}" )
            raise SModelSError ( f"no model defined for regionSet {regionSet} in dataset {dataset.globalInfo.id}" )

        model_type, model_filename = modelList[0] # Always use first model
        if model_type != "pyhf" and model_type != "full_pyhf":
            logger.error ( f"model type {model_type} for regionSet {regionSet} in dataset {dataset.globalInfo.id} is not 'pyhf/full_pyhf'" )
            raise SModelSError ( f"model type {model_type} for regionSet {regionSet} in dataset {dataset.globalInfo.id} is not 'pyhf/full_pyhf'" )

        # Get dictionary for signal yields using the pyhf labels
        nsignals = {label : nsigDict.get(label,0.0) for label in labelToPyhf.keys()}
        json = globalInfo.cachedModels[model_filename]
        regions = [regionMappings[label] for label in globalInfo.regionSets[regionSet]]
        logger.debug(f"list of datasets: {list(labelToPyhf.keys())}")


        includeCRs = False
        if hasattr(globalInfo,'includeCRs'):
            includeCRs = globalInfo.includeCRs
        signalUncertainty = None
        if hasattr(globalInfo,"signalUncertainty"):
            signalUncertainty = globalInfo.signalUncertainty

        # Loading the jsonFiles, unless we already have them (because we pickled)
        data = PyhfData(nsignals, json, regions,
                        includeCRs, signalUncertainty, globalInfo,
                        jsonFileName = model_filename )
        upperLimitComputer = PyhfUpperLimitComputer( data,
                                                    lumi=dataset.getLumi() )

        return upperLimitComputer

    @classmethod
    def forTruncatedGaussian(cls,theorypred: object, corr : float =0.6 ) -> Union[None,TruncatedGaussians]:
        """ get a sub computer for truncated gaussians
        :param theorypred: TheoryPrediction object
        :param corr: correction factor: \
                ULexp_mod = ULexp / (1. - corr*((ULobs-ULexp)/(ULobs+ULexp))) \
                a factor of corr = 0.6 is proposed.
        :returns: list of subComputers (with a single entry)
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
            element=theorypred.avgElement, txnames=theorypred.txnames,
            evaluationType=observed) / theorypred.xsection
        kwargs = { "upperLimitOnMu": float(ul),
                   "expectedUpperLimitOnMu": float(eul),
                   "corr": corr }
        computer = TruncatedGaussians ( **kwargs )
        return computer

    @classmethod
    def forAnalysesComb(cls,theoryPredictions: list, deltas_rel : Optional[float]) \
            -> AnaCombLikelihoodComputer:
        """ get a sub computer for combination of analyses
        :param theoryPredictions: list of TheoryPrediction objects
        :param deltas_rel: relative error for the signal
        :returns: a sub computer
        """
        # Only allow negative signal if all theory predictions allow for it
        allowNegativeSignals = all([tp.statsComputer.allowNegativeSignals
                                    for tp in theoryPredictions])

        computer = AnaCombLikelihoodComputer( theoryPredictions=theoryPredictions,
                                              deltas_rel=deltas_rel )
        computer.allowNegativeSignals = allowNegativeSignals
        return computer



class StatsComputer:
    """ this is the stats computer, it takes the subcomputers
    upon construction, and handles all the delegations to them,
    getting the most sensitive model, etc
    """

    def __init__ ( self, subComputers : list, allowNegativeSignals : bool = False ):
        """
         Initialise. allowNegativeSignals is true if its true for all
         subcomputers.
        """
        self.subComputers = subComputers
        self.allowNegativeSignals = allowNegativeSignals
        for computer in self.subComputers:
            computer.allowNegativeSignals = self.allowNegativeSignals

    @classmethod
    def forTheoryPrediction(cls, theoryPrediction: object) -> Union[None,'StatsComputer']:
        """Construct a StatsComputer from a TheoryPrediction.

        Inspects the data type and statistical model of the theory
        prediction and delegates to the appropriate sub-computer factory
        (SL, pyhf, NN, truncated Gaussian, or analyses combination).

        :param theoryPrediction: a TheoryPrediction or TheoryPredictionsCombiner object
        :returns: a StatsComputer instance, or None if no suitable computer can be built
        """

        CompRetriever = getCompRetrieverModule()

        dataType = theoryPrediction.dataType()
        tpType = theoryPrediction.type()
        computers = []

        if dataType == "upperLimit":
            from smodels.base.runtime import experimentalFeature
            if experimentalFeature( "truncatedGaussians" ):
                computer = CompRetriever.forTruncatedGaussian(theoryPrediction)
                if computer is not None:
                    computers.append(computer)

        elif dataType == "efficiencyMap":
            dataset = theoryPrediction.dataset
            nsigDict = {dataset.getID() : (theoryPrediction.xsection * dataset.getLumi()).asNumber()}
            computer = CompRetriever.forSingleBin(regionSet=theoryPrediction.dataset.getID(),dataset=theoryPrediction.dataset,
                                                nsigDict=nsigDict)
            computers.append(computer)

        elif dataType == "combined" and tpType == "TheoryPredictionsCombiner":
           # First make sure all theory predictions in the combiner have well-defined stats models
            if all(tp.statsComputer != 'N/A' for tp in theoryPrediction.theoryPredictions):
                computer = CompRetriever.forAnalysesComb(theoryPrediction.theoryPredictions,
                                                            theoryPrediction.deltas_rel)
                computers.append(computer)

        elif dataType == "combined" and tpType == "TheoryPrediction":
            dataset = theoryPrediction.dataset
            srNsigDictAll = { p.dataset.getID() : \
                (p.xsection*dataset.getLumi()).asNumber() \
                for p in theoryPrediction.datasetPredictions }

            for regionSet,modelList in dataset.globalInfo.statModels.items():
                if regionSet not in dataset.globalInfo.regionSets:
                    logger.error(f"A statistical model has been defined for {regionSet}, but it has not been found in regionSets")
                    raise ValueError(f"A statistical model has been defined for {regionSet}, but it has not been found in regionSets")

                if not modelList or len(modelList) == 0:
                    continue

                # Get the dict of signal yields for the given set of SRs:
                # (if the SR does not appear in theory predictions, set its signal yield to 0)
                srNsigDict = {sr: srNsigDictAll.get(sr, 0.0)
                            for sr in dataset.globalInfo.regionSets[regionSet]}

                # Always use the first model:
                model_type,_ = modelList[0]
                if model_type == "sl":
                    computers.append(CompRetriever.forMultiBinSL(regionSet=regionSet,dataset=dataset,
                                                            nsigDict=srNsigDict,
                                                            deltas_rel = theoryPrediction.deltas_rel ))
                elif model_type == "onnx":
                    computers.append(CompRetriever.forNNs(regionSet=regionSet,dataset=dataset,
                                                    nsigDict=srNsigDict))
                elif model_type == "pyhf":
                    computers.append(CompRetriever.forPyhf(regionSet=regionSet,dataset=dataset,
                                                    nsigDict=srNsigDict))
                else:
                    logger.error(f"Unknown statistical model type {model_type} for regionSet {regionSet} in dataset {dataset}")
                    raise SModelSError(f"Unknown statistical model type {model_type} for regionSet {regionSet} in dataset {dataset}")
        else:
            logger.error(f"Unknown data type {dataType} and type {tpType} for theory prediction {theoryPrediction}")
            raise SModelSError(f"Unknown data type {dataType} and type {tpType} for theory prediction {theoryPrediction}")


        if not isinstance(computers,list) or len(computers) == 0:
            return None
        else:
            allowNegativeSignals = all([comp.allowNegativeSignals for comp in computers])
            return cls(subComputers=computers, allowNegativeSignals=allowNegativeSignals)

    def get_five_values ( self, evaluationType : NllEvalType,
                      return_nll : bool = False,
                      check_for_maxima : bool = False )-> Dict:
        """
        Return the Five Values: l(bsm), l(sm), muhat, l(muhat), sigma(mu_hat)

        :param check_for_maxima: if true, then check lmax against l(sm) and l(bsm)
                                 correct, if necessary
        """
        assert return_nll == True, f"get_five_values return_nll {return_nll}"
        ret = self.nll_min ( evaluationType = evaluationType )
        if ret is None:
            return {}
        nll_min = ret[ "nll_min" ]

        nllbsm = self.nll ( poi_test = 1., evaluationType=evaluationType )
        nllsm = self.nll ( poi_test = 0., evaluationType=evaluationType )
        ret["nllbsm"] = nllbsm
        ret["nllsm"] = nllsm
        if check_for_maxima:
            if nllsm < nll_min: ## if return_nll is off, its the other way
                muhat = ret["muhat"]
                logger.debug(f"nllsm={nllsm:.2g} < nll_min({muhat:.2g})={nll_min:.2g}: will correct")
                ret[ "nll_min" ] = nllsm
                ret["muhat"] = 0.0
            if nllbsm < nll_min:
                muhat = ret["muhat"]
                logger.debug(f"nllbsm={nllbsm:.2g} < nll_min({muhat:.2g})={nll_min:.2g}: will correct")
                ret[ "nll_min" ] = nllbsm
                ret["muhat"] = 1.0
        return ret

    def nll ( self, poi_test : float, evaluationType : NllEvalType,
              asimov : Union[None,float] = None, **kwargs  ) -> Union[None,float]:
        """ simple frontend to individual computers """
        msm = self.getMostSensitiveModel()
        self.transform ( evaluationType )
        kwargs.update ( { "evaluationType": evaluationType, "asimov": asimov } )
        # kwargs = { "evaluationType": evaluationType, "asimov": asimov }
        if msm is None:
            return None
        ret = msm.nll ( poi_test, **kwargs)
        return ret

    def likelihood ( self, poi_test : float, evaluationType : NllEvalType,
                  return_nll : bool, asimov : Union[None,float] = None,
                  **kwargs ) -> Union[None,float]:
        """ convenience function, should become obsolete longterm """
        nll = self.nll ( poi_test, evaluationType, asimov, **kwargs )
        return exponentiateNLL ( nll, doIt = not return_nll )

    def CLs ( self, poi_test : float = 1.,
              evaluationType : NllEvalType=observed,
              **kwargs ) -> Union[float,None]:
        """ compute CLs value for a given value of the poi """
        msm = self.getMostSensitiveModel()
        # print ( f"@@ST0 getMostSensitiveModel eType {evaluationType} ret {ret}" )
        if msm is None:
            return None

        if hasattr ( msm , "CLs" ):
            return msm.CLs ( poi_test,
                    evaluationType = evaluationType, **kwargs )
        return None

    def transform ( self, evaluationType: NllEvalType ):
        """ SL only. transform the data to evaluationType or observed """
        for subComputer in self.subComputers:
            if subComputer is None:
                continue
            if subComputer.dataType in [ "pyhf", "truncGaussian", "analysesComb", "nn" ]:
                continue
            subComputer.likelihoodComputer.transform ( evaluationType )

    def restore ( self, evaluationType: NllEvalType ):
        """ SL only. Restore the data to the original observed values """
        if evaluationType != observed:
            return
        for subComputer in self.subComputers:
            if subComputer is None:
                continue
            if subComputer.dataType in [ "pyhf", "truncGaussian", "analysesComb" ]:
                continue
            subComputer.model = subComputer.origModel

    def getLlhds(self, **kwargs ) -> dict:
        """
        Facility to access the likelihoods for the individual analyses and
        the combined likelihood.
        Returns a dictionary with the analysis IDs as keys and the likelihood
        values as values.  Mostly used for plotting the likelihoods.

        :param muvals: List with values for the signal strenth for which
        the likelihoods must be evaluated.
        :param idx: index of subcomputer
        :param evaluationType: returns the observed/priori expected/posteriori expected likelihood values.
        :param normalize: If True normalizes the likelihood by its integral
        over muvals.
        """
        idx = kwargs.pop("idx",0)
        if idx >= len(self.subComputers):
             logger.error(f"only {len(self.subComputers)} computers for prediction but index {idx} was requested")
             raise SModelSError(f"only {len(self.subComputers)} computers for prediction index {idx} was requested")
        
        return self.subComputers[idx].getLlhds( **kwargs )

    def nll_min ( self, evaluationType : NllEvalType, ** kwargs ) -> Union[None,dict]:
        """
        :returns: dictionary with muhat, sigma_mu and nll_min as keys
        """
        msm = self.getMostSensitiveModel()
        if msm is None:
            return { "nll_min": float("nan" ), "mu_hat": float("nan"),
                     "sigma_mu": float("nan") }
        self.transform ( evaluationType )

        ret = msm.nll_min (
            evaluationType = evaluationType,
            allowNegativeSignals = msm.allowNegativeSignals, **kwargs )
        return ret

    @lru_cache
    def getMostSensitiveModel ( self ) -> Union[PyhfUpperLimitComputer,NNUpperLimitComputer,\
                                                SLUpperLimitComputer,AnaCombLikelihoodComputer,\
                                                TruncatedGaussians,None]:
        """ convenience function to get the most significant model

        :returns: dictionary with idx of the computer, ul_min,
        limit_on_xsecs
        and the name of the most sensitive model
        """
        ul_min = float("inf")
        most_sensitive_computer = None
        for i,computer in enumerate ( self.subComputers ):
            ul = computer.getUpperLimitOnMu ( evaluationType=apriori )
            if ul is not None and ul < ul_min:
                ul_min = ul
                most_sensitive_computer = computer
        return most_sensitive_computer

    def getTotalXSec ( self ) -> UnitXSec:
        """ get the total yield, summing over all computers """
        ret = 0.*fb
        for computer in self.subComputers:
            add = computer.getTotalXSec()
            ret += add
        return ret

    def getUpperLimit ( self, evaluationType : NllEvalType,
           limit_on_xsec : bool = False,
           nSigma : int = 0, **kwargs ) -> Union[float,UnitXSec,None]:
        """
        Simple frontend to the upperlimit computers, later to spey.poi_upper_limit

        :param limit_on_xsec: if True, then return the limit on the
        cross section
        :param nSigma: the upper limit for central value (0),
        + 1 sigma, - 1 sigma, etc. For error bands.
        :param kwargs: e.g. pmSigma
        """
        msm = self.getMostSensitiveModel()
        if msm is None:
            return None
        ulmu = msm.getUpperLimitOnMu(
                   evaluationType = evaluationType, nSigma = nSigma, **kwargs )
        if ulmu == None or not limit_on_xsec:
            return ulmu
        ret = ulmu * self.getTotalXSec()
        return ret

class SimpleStatsDataSet:
    """ a very simple data class that can replace a smodels.dataset,
    for 1d SL data only. used for testing and in dataPreparation """
    class SimpleInfo:
        def __init__ ( self, observedN : float, expectedBG : float,
                       bgError : float ):
            self.observedN = observedN
            self.expectedBG = expectedBG
            self.bgError = bgError

    class GlobalInfo:
        def __init__ ( self, lumi: UnitLumi ):
            self.id = "SimpleStatsDataSet"
            self.lumi = lumi

    def __init__ ( self, observedN : float, expectedBG : float,
                   bgError : float, lumi: UnitLumi = 1.0*fb ):
        """ initialise the dataset with all relevant stats """
        self.dataInfo = self.SimpleInfo ( observedN, expectedBG, bgError )
        self.globalInfo = self.GlobalInfo( lumi )

    def getLumi ( self ) -> UnitLumi:
        return self.globalInfo.lumi

    def getType ( self ) -> str:
        return "efficiencyMap"

if __name__ == "__main__":
    # nobs,bg,bgerr,lumi = 3., 4.1, 0.6533758489567854, 35.9/fb
    # nobs,bg,bgerr,lumi = 0, 1., 0.2, 35.9/fb
    # nobs,bg,bgerr,lumi = 0, 0.001, 0.01, 35.9/fb
    nobs,bg,bgerr,lumi = 3905,3658.3,238.767, 35.9/fb
    dataset = SimpleStatsDataSet ( nobs, bg, bgerr, lumi )
    computer = StatsComputer ( dataset, 1. )
    ul = computer.getUpperLimit ( evaluationType = observed,
                                    limit_on_xsec = True )
    print ( "ul", ul )
    ule = computer.getUpperLimit ( evaluationType = apriori,
                                     limit_on_xsec = True )
    print ( "ule", ule )
