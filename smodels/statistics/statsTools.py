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
    """ simple class that retrieves and constructs the sub computers,
    in StatsComputer """



    @classmethod
    def forMultiBinSL(cls,dataset, nsigDict, deltas_rel : float ) -> List[SLUpperLimitComputer]:
        """ get a subcomputer for simplified likelihood sr-combination.

        :param dataset: CombinedDataSet object
        :param nsigDict: Dictionary of signal events for each SR
        :deltas_rel: Relative uncertainty for the signal

        :returns: a subcomputer
        """
        covs = dataset.globalInfo.cachedModels
        offset = 0
        subComputers = []
        nsig = list(nsigDict.values())
        for covname,cov in covs.items():
            if not covname.endswith ( ".cov" ):
                continue
            if type(cov) != list:
                raise SModelSError( f"covariance field has wrong type: {type(cov)}" )
            if len(cov) < 1:
                raise SModelSError( f"covariance matrix has length {len(cov)}." )
            n = len(cov)

            # FIXME are we sure thats right, its not the one below?
            #nobs = [ x.dataInfo.observedN for x in dataset.origdatasets[offset:offset+n] ]
            #bg = [ x.dataInfo.expectedBG for x in dataset.origdatasets[offset:offset+n] ]
            #third_momenta = [ getattr ( x.dataInfo, "thirdMoment", None ) for x in dataset.origdatasets[offset:offset+n] ]
            nobs = [ x.dataInfo.observedN for x in dataset._datasets[offset:offset+n] ]
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
            subComputers.append ( computer )
            offset += n
        return subComputers

    @classmethod
    def forSingleBin( cls, dataset, nsigDict, deltas_rel : float = 0.2,
                      lumi : Optional[UnitLumi]=None ) ->  List[SLUpperLimitComputer]:
        """ get a sub computer for an efficiency map (single bin).

        :param dataset: DataSet object
        :param nsigDict: Dictionary of signal events for each SR
        :deltas_rel: Relative uncertainty for the signal

        :returns: a sub computer
        """
        data = Data( dataset.dataInfo.observedN, dataset.dataInfo.expectedBG,
                     dataset.dataInfo.bgError**2, deltas_rel = deltas_rel,
                     nsignal = list(nsigDict.values()), lumi = lumi )
        likelihoodComputer = LikelihoodComputer ( data )
        computer = SLUpperLimitComputer ( likelihoodComputer )
        computer.dataType = "1bin"
        computer.allowNegativeSignals = False
        return [ computer ]

    @classmethod
    def forNNs(cls, dataset, nsigDict ) -> List[NNUpperLimitComputer]:
        """ get a sub computer for an NN combination.

        :param dataset: CombinedDataSet object
        :param nsigDict: Dictionary of signal events for each SR
        :deltas_rel: Relative uncertainty for the signal

        :returns: a sub computer
        """
        globalInfo = dataset.globalInfo
        labelToONNX = {}
        labelToSModelS = {}

        for sr in globalInfo.srMappings:
           # nsignals[ sr["onnx"] ] = 0.
            if sr["label"] != None:
                labelToONNX [ sr["label"] ] = sr["onnx"]
                labelToSModelS [ sr["label"] ] = sr["smodels"]

        
        subComputers = []
        for srSetName, model_tuples in globalInfo.statModels.items():
            f_signals = {}
            for sr in globalInfo.srMappings:
                f_signals[ sr["onnx"] ] = 0.
            for label in globalInfo.srSets[srSetName]:
                smodelsName = labelToSModelS[label]
                if smodelsName in nsigDict:
                    f_signals[ labelToONNX[label] ] = \
                        nsigDict[ smodelsName ]
            model_tuple = model_tuples[0]
            modelfilename = model_tuple[1]
            model_type = model_tuple[0]
            if "pyhf" in model_type:
                continue
            data = NNData( f_signals, dataset )
            upperLimitComputer = NNUpperLimitComputer(data,
                    lumi=dataset.getLumi(),
                    onnxfilename = modelfilename )
            upperLimitComputer.allowNegativeSignals = False
            subComputers.append ( upperLimitComputer )
        return subComputers

    @classmethod
    def forPyhf(cls, dataset, nsigDict) -> List[PyhfUpperLimitComputer]:
        """ get a sub computer for pyhf combination.

        :param dataset: CombinedDataSet object
        :param nsigDict: Dictionary with signal yields for each SR

        :returns: a sub computer
        """
        globalInfo = dataset.globalInfo
        datasets = [ds.getID() for ds in dataset.origdatasets]
        jsonFiles, jsons = [], []
        jsonDictNames = {}
        for srSetName, model_tuples in globalInfo.statModels.items():
            model_tuple = model_tuples[0]
            model_type = model_tuple[0]
            model_name = model_tuple[1]
            if not "pyhf" in model_type:
                continue
            jsonFiles.append ( model_name )
            jsons.append ( globalInfo.cachedModels[model_name] )
            jsonDictNames[model_name]=cls.getAll ( \
                    globalInfo.srSets[srSetName], globalInfo.srMappingsDict )

        logger.debug(f"list of datasets: {datasets}")

        # Constructing the list of signals with subsignals matching each json
        nsignals = {}
        for jsName in jsonFiles:
            nsignals.update( { jsName: {} } )
        for name, nsig in nsigDict.items():
            for jsName in nsignals.keys():
                if name in jsonDictNames[jsName]:
                    nsignals[jsName].update( { name: nsig } )

        includeCRs = False
        if hasattr(globalInfo,'includeCRs'):
            includeCRs = globalInfo.includeCRs
        signalUncertainty = None
        if hasattr(globalInfo,"signalUncertainty"):
            signalUncertainty = globalInfo.signalUncertainty

        r_jsonFiles = {} ## here we reconstruct the old jsonFiles dict
        ## (for now, maybe we will see that there is a smarter way)
        for srSetName, model_tuples in globalInfo.statModels.items():
            model_tuple = model_tuples[0]
            if not model_tuple[0] in [ "pyhf", "full_pyhf" ]:
                continue
            region_names = globalInfo.srSets[srSetName]
            regions = cls.getRegions ( region_names,
                    globalInfo.srMappingsDict )
            if model_tuple[1] in r_jsonFiles:
                raise SModelSError ( f"model {model_tuple[1]} mentioned more than once in {globalInfo.path}" )
            r_jsonFiles[model_tuple[1]]=regions

        subComputers = []
        for nsignal,json,(jsonFileName,r_jsonFile) in zip(nsignals.values(),jsons,r_jsonFiles.items() ):
            # Loading the jsonFiles, unless we already have them (because we pickled)
            data = PyhfData(nsignal, json, r_jsonFile,
                            includeCRs, signalUncertainty, globalInfo,
                            jsonFileName = jsonFileName )
            upperLimitComputer = PyhfUpperLimitComputer( data,
                                                         lumi=dataset.getLumi() )
            upperLimitComputer.allowNegativeSignals = False
            subComputers.append ( upperLimitComputer )
        return subComputers

    @classmethod
    def forTruncatedGaussian(cls,theorypred, corr : float =0.6 ) -> List[TruncatedGaussians]:
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
            return []

        eul = theorypred.dataset.getUpperLimitFor(
            element=theorypred.avgElement, txnames=theorypred.txnames, evaluationType=apriori
        )
        if eul is None:
            return []
        eul = eul / theorypred.xsection
        ul = theorypred.dataset.getUpperLimitFor(
            element=theorypred.avgElement, txnames=theorypred.txnames, 
            evaluationType=observed) / theorypred.xsection
        kwargs = { "upperLimitOnMu": float(ul), 
                   "expectedUpperLimitOnMu": float(eul),
                   "corr": corr }
        computer = TruncatedGaussians ( **kwargs )
        return [ computer ]

    @classmethod
    def forAnalysesComb(cls,theoryPredictions, deltas_rel : Optional[float]) \
            -> List[AnaCombLikelihoodComputer]:
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
        return [ computer ]

    @classmethod
    def getAll ( cls, srSet : list, srMappingsDict : dict ) -> list:
        """ for the srSet, return list only of signal regions,
        drop others """
        ret = []
        for r in srSet:
            if r in srMappingsDict:
                m = srMappingsDict[r]
                if True: # m["type"] == "SR":
                    ret.append ( r )
        return ret

    @classmethod
    def getRegions ( cls, region_labels : list, srMappingsDict : dict ) -> list:
        """ given a list of the region_labels, return the
        corresponding list of the region dictionaries """
        ret = []
        for label in region_labels:
            if not label in srMappingsDict:
                continue
            ret.append ( srMappingsDict[label] )
        return ret

class StatsComputer:
    """ this is the stats computer, it takes the subcomputers
    upon construction, and handles all the delegations to them,
    getting the most sensitive model, etc
    """

    def __init__ ( self, subComputers : list ):
        """
         Initialise. allowNegativeSignals is true if its true for all
         subcomputers.
        """
        self.subComputers = subComputers
        self.allowNegativeSignals = True
        for computer in subComputers:
            if not self.allowNegativeSignals:
                self.allowNegativeSignals = False
                break

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
