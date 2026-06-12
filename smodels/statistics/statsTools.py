#!/usr/bin/env python3

"""
.. module:: statsTools
   :synopsis: a module that contains the class responsible for
              all statistical computations. Designed to
              eventually become simply a frontend for spey

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

__all__ = [ "StatsComputer", "getCompRetrieverModule" ]

from typing import Union, Dict, Optional, List
from smodels.statistics.exceptions import SModelSStatisticsError as SModelSError
from smodels.base.smodelsLogging import logger
from smodels.base.physicsUnits import fb, UnitLumi
from smodels.statistics.simplifiedLikelihoods import LikelihoodComputer, SLUpperLimitComputer, Data
from smodels.statistics.pyhfInterface import PyhfData, PyhfUpperLimitComputer
from smodels.statistics.nnInterface import NNData, NNUpperLimitComputer
from smodels.statistics.basicStats import observed, apriori, NllEvalType, exponentiateNLL
from smodels.statistics.truncatedGaussians import TruncatedGaussians
from smodels.statistics.analysesCombinations import AnaCombLikelihoodComputer
from smodels.base.physicsUnits import UnitXSec
from smodels.tools.caching import lru_cache

def getCompRetrieverModule():
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
    def forMultiBinSL(cls,srSet : str, dataset, nsigDict, 
            deltas_rel : Optional[float] = 0.0  ) -> SLUpperLimitComputer:
        """ get a subcomputer for simplified likelihood sr-combination.

        :param dataset: CombinedDataSet object
        :param srSet: the srSet defining the SRs for which to construct the computer
        :param nsigDict: Dictionary of signal events for each SR
        :deltas_rel: Relative uncertainty for the signal

        :returns: a subcomputer
        """
        assert srSet in dataset.globalInfo.statModels, f"{srSet} not in statModels in {dataset.globalInfo.id}"
        covs = dataset.globalInfo.cachedModels
        type_n_models = dataset.globalInfo.statModels[srSet]
        type_n_model = type_n_models[0] # get first one
        mtype = type_n_model[0]
        assert mtype == "sl", f"expected sl but got {mtype} for type of stats model"
        covname = type_n_model[1]
        #assert covname.endswith(".cov"), f"name of covariance matrix file {covname} does not end with .cov in {dataset.globalInfo.id}"
        offset = 0
        #subComputers = []
        srList = dataset.globalInfo.srSets[srSet]
        nsig = [nsigDict.get(sr, 0.0) for sr in srList ]
        cov = covs[covname]
        assert type(cov)==list, f"covariance field has wrong type: {type(cov)} in {dataset.globalInfo.id}"
        assert len(cov)>0, f"covariance matrix has length {len(cov)}."
        n = len(cov)

        nobs = [ x.dataInfo.observedN for x in dataset._datasets[offset:offset+n] ] ## [!AL!]: do we need _datasets? Can't we just loop over the SRs defined in srSets[seSet]? WW: we need e.g. observedN, this is the way to get it
        bg = [ x.dataInfo.expectedBG for x in dataset._datasets[offset:offset+n] ]
        nsig = nsig[offset:offset+n]
        third_momenta = [ getattr ( x.dataInfo, "thirdMoment", None ) for x in dataset._datasets[offset:offset+n] ]
        c = third_momenta.count ( None )
        if c > 0:
            if c < len(third_momenta):
                logger.warning ( f"third momenta given for some but not all signal regions in {dataset.globalInfo.id}" )
            third_momenta = None

        data = Data( nobs, bg, cov, third_moment=third_momenta,
                     nsignal = nsig,
                     deltas_rel = deltas_rel, lumi=dataset.getLumi(),
                     name = covname )
        likelihoodComputer = LikelihoodComputer ( data )
        computer = SLUpperLimitComputer ( likelihoodComputer )
        computer.dataType = "SL"
        computer.allowNegativeSignals = False
        return computer

    @classmethod
    def forSingleBin( cls, srSet :str, dataset, nsigDict, deltas_rel : float = 0.2,
                      lumi : Optional[UnitLumi]=None ) ->  SLUpperLimitComputer:
        """ get a sub computer for an efficiency map (single bin).

        :param dataset: DataSet object
        :param srSet: the srSet defining the SRs for which to construct the computer
        :param nsigDict: Dictionary of signal events for each SR
        :deltas_rel: Relative uncertainty for the signal

        :returns: a sub computer
        """
        if srSet != dataset.getID():
            logger.error ( f"srSet {srSet} does not match dataset id {dataset.getID()}" )
            raise SModelSError ( f"srSet {srSet} does not match dataset id {dataset.getID()}" )
        data = Data( dataset.dataInfo.observedN, dataset.dataInfo.expectedBG,
                     dataset.dataInfo.bgError**2, deltas_rel = deltas_rel,
                     nsignal = list(nsigDict.values()), lumi = lumi )
        likelihoodComputer = LikelihoodComputer ( data )
        computer = SLUpperLimitComputer ( likelihoodComputer )
        computer.dataType = "1bin"
        computer.allowNegativeSignals = False
        return computer

    @classmethod
    def forNNs(cls, srSet : str, dataset, nsigDict ) -> NNUpperLimitComputer:
        """ get a sub computer for an NN combination.

        :param srSet: the srSet defining the SRs for which to construct the computer
        :param dataset: CombinedDataSet object
        :param nsigDict: Dictionary of signal events for each SR
        :deltas_rel: Relative uncertainty for the signal

        :returns: a sub computer
        """
        globalInfo = dataset.globalInfo
        labelToONNX = {}
        srMappingsDict = globalInfo.srMappingsDict

        for sr in globalInfo.srSets[srSet]:
            if sr not in srMappingsDict:
                logger.error ( f"SR {sr} defined in srSet {srSet} not found in srMappings for dataset {dataset.globalInfo.id}" )
                raise SModelSError ( f"SR {sr} defined in srSet {srSet} not found in srMappings for dataset {dataset.globalInfo.id}" )
            labelToONNX[sr] = srMappingsDict[sr]["onnx"]

        if srSet not in globalInfo.statModels:
            logger.error ( f"srSet {srSet} not found in statModels for dataset {dataset.globalInfo.id}" )
            raise SModelSError ( f"srSet {srSet} not found in statModels for dataset {dataset.globalInfo.id}" )

        modelList = globalInfo.statModels[srSet]
        if len(modelList) == 0:
            logger.error ( f"no model defined for srSet {srSet} in dataset {dataset.globalInfo.id}" )
            raise SModelSError ( f"no model defined for srSet {srSet} in dataset {dataset.globalInfo.id}" )

        model_type, model_filename = modelList[0] # Always use first model
        if model_type != "onnx":
            logger.error ( f"model type {model_type} for srSet {srSet} in dataset {dataset.globalInfo.id} is not 'onnx'" )
            raise SModelSError ( f"model type {model_type} for srSet {srSet} in dataset {dataset.globalInfo.id} is not 'onnx'" )

        # Get dictionary for signal yields using the ONNX labels
        f_signals = {onnx_sr : nsigDict.get(label,0.0) for label,onnx_sr in labelToONNX.items()}

        data = NNData( f_signals, dataset )
        upperLimitComputer = NNUpperLimitComputer(data, lumi=dataset.getLumi(),
                                                  onnxfilename = model_filename )

        return upperLimitComputer

    @classmethod
    def forPyhf( cls, srSet : str, dataset : CombinedDataSet,
                 nsigDict : dict ) -> PyhfUpperLimitComputer:
        """ get a sub computer for pyhf combination.

        :param srSet: name of signal region set
        :param dataset: CombinedDataSet object
        :param nsigDict: Dictionary with signal yields for each SR

        :returns: a sub computer
        """

        globalInfo = dataset.globalInfo
        srMappingsDict = globalInfo.srMappingsDict
        labelToPyhf = {}
        for sr_label in globalInfo.srSets[srSet]:
            if sr_label not in srMappingsDict:
                logger.error ( f"SR {sr_label} defined in srSet {srSet} not found in srMappings for dataset {dataset.globalInfo.id}" )
                raise SModelSError ( f"SR {sr_label} defined in srSet {srSet} not found in srMappings for dataset {dataset.globalInfo.id}" )
            if srMappingsDict[sr_label]["smodels"] is not None:
                labelToPyhf[sr_label] = srMappingsDict[sr_label]["pyhf"]

        if srSet not in globalInfo.statModels:
            logger.error ( f"srSet {srSet} not found in statModels for dataset {dataset.globalInfo.id}" )
            raise SModelSError ( f"srSet {srSet} not found in statModels for dataset {dataset.globalInfo.id}" )

        modelList = globalInfo.statModels[srSet]
        if len(modelList) == 0:
            logger.error ( f"no model defined for srSet {srSet} in dataset {dataset.globalInfo.id}" )
            raise SModelSError ( f"no model defined for srSet {srSet} in dataset {dataset.globalInfo.id}" )

        model_type, model_filename = modelList[0] # Always use first model
        if model_type != "pyhf" and model_type != "full_pyhf":
            logger.error ( f"model type {model_type} for srSet {srSet} in dataset {dataset.globalInfo.id} is not 'pyhf/full_pyhf'" )
            raise SModelSError ( f"model type {model_type} for srSet {srSet} in dataset {dataset.globalInfo.id} is not 'pyhf/full_pyhf'" )

        # Get dictionary for signal yields using the pyhf labels
        nsignals = {label : nsigDict.get(label,0.0) for label in labelToPyhf.keys()}
        json = globalInfo.cachedModels[model_filename]
        regions = [srMappingsDict[label] for label in globalInfo.srSets[srSet]]
        logger.debug(f"list of datasets: {list(labelToPyhf.keys())}")


        includeCRs = False
        if hasattr(globalInfo,'includeCRs'):
            includeCRs = globalInfo.includeCRs
        signalUncertainty = None
        if hasattr(globalInfo,"signalUncertainty"):
            signalUncertainty = globalInfo.signalUncertainty

        # Loading the jsonFiles, unless we already have them (because we pickled)
        # ## [!AL!]: I've tried to simplify the code here, since we should return a single computer. But I am not sure if the input for PyhfData is correct
        data = PyhfData(nsignals, json, regions,
                        includeCRs, signalUncertainty, globalInfo,
                        jsonFileName = model_filename )
        upperLimitComputer = PyhfUpperLimitComputer( data,
                                                    lumi=dataset.getLumi() )

        return upperLimitComputer

    @classmethod
    def forTruncatedGaussian(cls,theorypred, corr : float =0.6 ) -> Union[None,TruncatedGaussians]:
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
    def forAnalysesComb(cls,theoryPredictions, deltas_rel : Optional[float]) \
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
    def forTheoryPrediction(cls, theoryPrediction) -> Union[None,'StatsComputer']:

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
            computer = CompRetriever.forSingleBin(srSet=theoryPrediction.dataset.getID(),dataset=theoryPrediction.dataset,
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
            # Get dictionary with dataset IDs and signal yields
            srNsigDictAll = {}
            for region in dataset.globalInfo.srMappings:
                srNsigDictAll[ region["smodels"] ] = 0.
            # Update with theory predictions
            srNsigDictAll.update({pred.dataset.getID() : (pred.xsection*dataset.getLumi()).asNumber()
                                    for pred in theoryPrediction.datasetPredictions})

            for srSet,modelList in dataset.globalInfo.statModels.items():
                if srSet not in dataset.globalInfo.srSets:
                    logger.error(f"A statistical model has been defined for {srSet}, but it has not been found in srSets")
                    raise ValueError(f"A statistical model has been defined for {srSet}, but it has not been found in srSets")

                if not modelList or len(modelList) == 0:
                    continue

                # Get the dict of signal yields for the given set of SRs:
                # (if the SR does not appear in theory predictions, set its signal yield to 0)
                srNsigDict = {sr: srNsigDictAll.get(sr, 0.0)
                            for sr in dataset.globalInfo.srSets[srSet]}

                # Always use the first model:
                model_type,_ = modelList[0]
                if model_type == "sl":
                    computers.append(CompRetriever.forMultiBinSL(srSet=srSet,dataset=dataset,
                                                            nsigDict=srNsigDict,
                                                            deltas_rel = theoryPrediction.deltas_rel ))
                elif model_type == "onnx":
                    computers.append(CompRetriever.forNNs(srSet=srSet,dataset=dataset,
                                                    nsigDict=srNsigDict))
                elif model_type == "pyhf":
                    computers.append(CompRetriever.forPyhf(srSet=srSet,dataset=dataset,
                                                    nsigDict=srNsigDict))
                else:
                    logger.error(f"Unknown statistical model type {model_type} for srSet {srSet} in dataset {dataset}")
                    raise SModelSError(f"Unknown statistical model type {model_type} for srSet {srSet} in dataset {dataset}")
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

    def transform ( self, evaluationType ):
        """ SL only. transform the data to evaluationType or observed """
        for subComputer in self.subComputers:
            if subComputer is None:
                continue
            if subComputer.dataType in [ "pyhf", "truncGaussian", "analysesComb", "nn" ]:
                continue
            subComputer.likelihoodComputer.transform ( evaluationType )

    def restore ( self, evaluationType ):
        """ SL only. transform the data to evaluationType or observed """
        if evaluationType != observed:
            return
        for subComputer in self.subComputers:
            if subComputer is None:
                continue
            if subComputer.dataType in [ "pyhf", "truncGaussian", "analysesComb" ]:
                continue
            subComputer.model = subComputer.origModel

    #def getLlhds(self,muvals,evaluationType : bool = False,
    #              normalize : bool = True,
    #              idx : int = 0 ) -> dict:
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
        idx = 0
        if "idx" in kwargs:
            idx = kwargs.pop("idx")
        assert idx < len(self.subComputers), f"only {len(self.subComputers)} computers for prediction but you wanted index {idx}"
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

    def getTotalXSec ( self ):
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
        def __init__ ( self, lumi ):
            self.id = "SimpleStatsDataSet"
            self.lumi = lumi

    def __init__ ( self, observedN : float, expectedBG : float,
                   bgError : float, lumi = 1.0*fb ):
        """ initialise the dataset with all relevant stats """
        self.dataInfo = self.SimpleInfo ( observedN, expectedBG, bgError )
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
    ul = computer.getUpperLimit ( evaluationType = observed,
                                    limit_on_xsec = True )
    print ( "ul", ul )
    ule = computer.getUpperLimit ( evaluationType = apriori,
                                     limit_on_xsec = True )
    print ( "ule", ule )
