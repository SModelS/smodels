#!/usr/bin/env python3

"""
.. module:: speyTools
   :synopsis: a module that contains tools and convenience methods
              that we use in connection with spey.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

__all__ = [ "SpeyRetriever", "SpeyAnalysesCombosComputer" ]

from typing import Union, Text, Tuple, Dict, List, Optional
import sys
from spey import ExpectationType, StatisticalModel, get_backend
from smodels.statistics.basicStats import exponentiateNLL
from smodels.base.physicsUnits import UnitXSec
import spey
# spey.set_optimiser( "iminuit" )

try:
    from spey.system.exceptions import AsimovTestStatZero
except ImportError: # comes only with newer versions of spey
    AsimovTestStatZero = Exception # a dummy so we can still try
from smodels.base.smodelsLogging import logger
from smodels.base.physicsUnits import fb, UnitLumi
from smodels.experiment.datasetObj import DataSet
from smodels.statistics.basicStats import observed, apriori, aposteriori, NllEvalType
from smodels.base.crossSection import XSection
from smodels.base.physicsUnits import fb, UnitXSec
import numpy as np

_debug = { "writePoint": False } # for debugging only

class SpeyRetriever:
    @classmethod
    def forMultiBinSL(cls,dataset, nsig, deltas_rel : Optional[float] = 0.0 ) -> list:
        """ get a subcomputer for simplified likelihood sr-combination.

        :param dataset: CombinedDataSet object
        :param nsig: number of signal events. For simplified likelihood backend
        this input can contain `np.array` or `List[float]` which contains
        signal yields per region.
        For `pyhf` backend this input evaluationType to be a JSON-patch
        i.e. `List[Dict]`, see `pyhf` documentation for details on
        JSON-patch format.
        :deltas_rel: Relative uncertainty for the signal

        :returns: a subcomputer
        """
        obsN = [ x.dataInfo.observedN for x in dataset._datasets ]
        bg = [ x.dataInfo.expectedBG for x in dataset._datasets ]
        # cov = dataset.globalInfo.covariance
        covs = dataset.globalInfo.cachedModels
        subComputers = []
        lumi = float ( dataset.getLumi().asNumber(1./fb) )
        thirdmomenta=[]
        for ds in dataset._datasets:
            if hasattr ( ds.dataInfo, "thirdMoment" ):
                thirdmomenta.append ( ds.dataInfo.thirdMoment )

        for covname,cov in covs.items():
            if not covname.endswith ( ".cov" ):
                continue
            if type(cov) != list:
                raise SModelSError( f"covariance field has wrong type: {type(cov)}" )
            if len(cov) < 1:
                raise SModelSError( f"covariance matrix has length {len(cov)}." )
            n = len(cov)
            bg = [ x.dataInfo.expectedBG for x in dataset._datasets[offset:offset+n] ]
            nsig = nsig[offset:offset+n]
            third_momenta = [ getattr ( x.dataInfo, "thirdMoment", None ) for x in dataset._datasets[offset:offset+n] ]
            c = third_momenta.count ( None )
            if c > 0:
                if c < len(third_momenta):
                    logger.warning ( f"third momenta given for some but not all signal regions in {dataset.globalInfo.id}" )
                third_momenta = None
            xsec = sum ( nsig ) / dataset.getLumi()
            if third_momenta == None:
                try:
                    stat_wrapper = get_backend("default.correlated_background")
                except spey.PluginError as e: ## older spey?
                    stat_wrapper = get_backend("default_pdf.correlated_background")
                if _debug["writePoint"]:
                    f=open ( "data.txt","wt" )
                    f.write ( f"obsN={obsN}\n" )
                    f.write ( f"bg={bg}\n" )
                    f.write ( f"cov={cov}\n" )
                    f.write ( f"nsig={nsig}\n" )
                    f.write ( f"analysis='{dataset.globalInfo.id}'\n" )
                    f.write ( f"lumi={lumi}\n" )
                    f.close()
    # import sys; sys.exit()

                speyModel = stat_wrapper( data = obsN,
                                background_yields = bg, covariance_matrix = cov,
                                signal_yields = nsig,
                                xsection = [ x / lumi for x in nsig ],
                                analysis = dataset.globalInfo.id,
                )
                facade = SpeyModelFacade ( speyModel, "SL", covname, xsec )
                subComputers.append ( facade )
            else:
                # SLv2
                try:
                    stat_wrapper = get_backend("default.third_moment_expansion")
                except ImportError as e:
                    stat_wrapper = get_backend("default_pdf.third_moment_expansion")
                speyModel = stat_wrapper( data = obsN,
                                background_yields = bg, covariance_matrix = cov,
                                signal_yields = nsig,
                                xsection = [ x / lumi for x in nsig ],
                                third_moment = thirdmomenta,
                                analysis = dataset.globalInfo.id,
                )
                if _debug["writePoint"]:
                    f=open ( "data.txt","wt" )
                    f.write ( f"obsN={obsN}\n" )
                    f.write ( f"bg={bg}\n" )
                    f.write ( f"cov={cov}\n" )
                    f.write ( f"nsig={nsig}\n" )
                    f.write ( f"analysis='{dataset.globalInfo.id}'\n" )
                    f.write ( f"lumi={lumi}\n" )
                    f.close()
                facade = SpeyModelFacade ( speyModel, "SL", covname, xsec )
                subComputers.append ( facade )
        return subComputers

    @classmethod
    def forSingleBin( cls, dataset, nsig ) -> list:
        """ get a sub computer for an efficiency map (single bin).

        :param dataset: DataSet object
        :param nsig: Number of signal events for each SR
        :deltas_rel: Relative uncertainty for the signal

        :returns: a sub computer
        :raises NotImplementedError: If requested backend has not been recognised.
        """
        try:
            stat_wrapper = get_backend("default.uncorrelated_background")
        except spey.PluginError as e: ## older spey?
            stat_wrapper = get_backend("default_pdf.uncorrelated_background")
        id = f"{dataset.globalInfo.id}:{dataset.dataInfo.dataId}"

        speyModel = stat_wrapper(
                        data = [float(dataset.dataInfo.observedN)],
                        background_yields = [float(dataset.dataInfo.expectedBG)],
                        absolute_uncertainties = [float(dataset.dataInfo.bgError)],
                        signal_yields = [nsig],
                        xsection = float (nsig/dataset.getLumi().asNumber(1./fb)),
                        analysis = id,
#                        backend = 'simplified_likelihoods'
        )
        name = dataset.dataInfo.dataId
        xsec = nsig / dataset.globalInfo.lumi
        facade = SpeyModelFacade ( speyModel, "1bin", name, xsec )
        return [ facade ]

    @classmethod
    def forNNs(cls, dataset, nsig ) -> list:
        """ get a sub computer for an NN combination.

        :param dataset: CombinedDataSet object
        :param nsig: Number of signal events for each SR
        :deltas_rel: Relative uncertainty for the signal

        :returns: a sub computer
        """
        import spey
        stat_wrapper = spey.get_backend('ml.likelihoods')
        globalInfo = dataset.globalInfo
        statModels = globalInfo.statModels
        labelToONNX = {}
        labelToSModelS = {}
        hasJsonsWithoutMLModels = False

        for sr in globalInfo.srMappings:
           # nsignals[ sr["onnx"] ] = 0.
            if sr["label"] != None:
                labelToONNX [ sr["label"] ] = sr["onnx"]
                labelToSModelS [ sr["label"] ] = sr["smodels"]

        subComputers = cls.forPyhf ( dataset, nsig )

        for srSetName, model_tuples in globalInfo.statModels.items():
            f_signals = {}
            for sr in globalInfo.srMappings:
                f_signals[ sr["onnx"] ] = 0.
            for label in globalInfo.srSets[srSetName]:
                smodelsName = labelToSModelS[label]
                if smodelsName in nsig:
                    f_signals[ labelToONNX[label] ] = \
                        nsig[ smodelsName ]
            model_tuple = model_tuples[0]
            modelfilename = model_tuple[1]
            if "pyhf" in model_tuple[0]:
                continue
            onnxBlob=dataset.globalInfo.cachedModels[srSetName]
            # self.speyModel = stat_wrapper(nsig,onnxBlob) # this is how i want it long run
            ## the following code is just for now to see if it works in principle
            import tempfile
            tempf = tempfile.mktemp ( prefix="/tmp", postfix=".onnx" )
            # tempf = "/tmp/my.onnx"
            f = open ( tempf, "wb" )
            import onnx
            onnx.save ( onnxBlob, f )
            f.close()
            speyModel = stat_wrapper(nsig,tempf)
            if os.path.exists ( tempf ):
                os.unlink ( tempf )
                xsec = sum(nsig) / dataset.globalInfo.lumi
                facade = speyModelFacade ( upperLimitComputer, "nn", 
                                           modelfilename, xsec )
                subComputers.append ( facade )
        return subComputers

    @classmethod
    def forPyhf(cls, dataset, nsig) -> list:
        """ get a sub computer for pyhf combination.

        :param dataset: CombinedDataSet object
        :param nsig: Number of signal events for each SR

        :returns: a sub computer
        """
        stat_wrapper = get_backend("pyhf")
        from smodels.statistics.speyPyhf import SpeyPyhfData
        data = SpeyPyhfData.createDataObject ( dataset, nsig )
        models = []
        patches = data.patchMaker()
        for i in range( len(data.inputJsons ) ):
            # idx, _ = self.getBestCombinationIndex( data )
            inputJson = data.inputJsons[i]
            signal_patch = patches[i]
            #print ( "inputJsons", inputJson )
            # import IPython; IPython.embed( colors = "neutral" ); sys.exit()
            analysis = dataset.globalInfo.id

            speyModel = stat_wrapper( analysis = analysis,
                            signal_patch = signal_patch,
                            background_only_model = inputJson )
            xsec = sum (nsig ) / dataset.globalInfo.lumi
            facade = SpeyModelFacade ( speyModel, "pyhf", mname, xsec )
            models.append ( facade )
        return models

    @classmethod
    def forTruncatedGaussian(cls,theorypred, corr : float =0.6 ) -> list:
        """ get a sub computer for truncated gaussians
        :param theorypred: TheoryPrediction object
        :param corr: correction factor: \
                ULexp_mod = ULexp / (1. - corr*((ULobs-ULexp)/(ULobs+ULexp))) \
                a factor of corr = 0.6 is proposed.
        :returns: list of subComputers (with a single entry)
        """
        logger.error ( f"speyTools no truncated Gaussian backend exists" )
        return None

    @classmethod
    def forAnalysesComb(cls,theoryPredictions, deltas_rel : Optional[float]) \
            -> list:
        """ get a sub computer for combination of analyses
        :param theoryPredictions: list of TheoryPrediction objects
        :param deltas_rel: relative error for the signal
        :returns: a sub computer
        """
        from smodels.statistics.analysesCombinations import AnaCombLikelihoodComputer
        computer = AnaCombLikelihoodComputer( theoryPredictions=theoryPredictions,
                                              deltas_rel=deltas_rel )
        computer.dataType = "analysesComb"
        computer.allowNegativeSignals = False
        #computer.dataType = "analysesComb"
        #computer.allowNegativeSignals = allowNegativeSignals
        return [ computer ]

class SimpleSpeyDataSet:
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
            self.id = "SimpleSpeyDataSet"
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

class SpeyModelFacade:
    """ very simple container that wraps around spey models,
    and adapts the API to SModelS """

    def __init__ ( self, speyModel, dataType : str, name : str,
                   totalXsec: UnitXSec ):
        self.speyModel = speyModel
        self.dataType = dataType
        self.name = name
        self.likelihoodComputer = self
        self.allowNegativeSignals = False
        self.totalXsec = totalXsec

    def nll_min ( self, evaluationType : NllEvalType = observed,
                  allowNegativeSignals : bool = False,
                  **kwargs ) -> dict:
        kwargs["return_nll"]=True
        kwargs["allow_negative_signal"] = allowNegativeSignals
        speyret = self.speyModel.maximize_likelihood ( **kwargs )
        muhat = float(speyret[0])
        ret = { "muhat": muhat, "nll_min": float(speyret[1]) }
        ## if we need sigma_mu we need to add it here
        test_statistic = "q" if self.allowNegativeSignals else "qmutilde"
        exp = self.translateExpectationType ( evaluationType )
        #sigma_mu = self.speyModel.sigma_mu( poi_test=muhat,
        #        expected=exp, test_statistics=test_statistic )
        sigma_mu = 1.
        ret[ "sigma_mu" ] = sigma_mu
        return ret

    def nll ( self, mu : float, evaluationType : NllEvalType = observed,
              asimov : Optional[float] = None,
              **kwargs ):
        exp = self.translateExpectationType ( evaluationType )
        kwargs["return_nll"]=True
        ret = self.speyModel.likelihood ( poi_test=mu,
               expected = exp, **kwargs )
        return float(ret)

    def getTotalXSec ( self ):
        return self.totalXsec

    @classmethod
    def transform ( cls, evaluationType : NllEvalType ):
        return

    @classmethod
    def translateExpectationType ( cls, evaluationType : NllEvalType ) -> ExpectationType:
        """ translate the specification for evaluationType values from smodels
            lingo to spey convention """
        if type(evaluationType)==ExpectationType:
            return expected
        expectedDict = { observed: ExpectationType.observed,
                         apriori: ExpectationType.apriori,
                         aposteriori: ExpectationType.aposteriori}
        if evaluationType in expectedDict:
            return expectedDict[evaluationType]
        logger.error( f'{expected} is not a valid expectation type. Possible expectation types are True (observed), False (apriori) and "posteriori".' )
        return None

    def getUpperLimitOnMu ( self, evaluationType : NllEvalType = observed,
            nSigma : int = 0, **kwargs ) -> float:
        exp = self.translateExpectationType ( evaluationType )
        expected_pvalue = "nominal"
        if nSigma != 0:
            expected_pvalue = "1sigma"

        ret = self.speyModel.poi_upper_limit ( expected = exp,
               expected_pvalue = expected_pvalue )
        if nSigma == 0:
            ret = float ( ret )
        elif nSigma == 1:
            ret = float(ret[0])
        elif nSigma == -1:
            ret = float(ret[-1])
        return ret

if __name__ == "__main__":
    # nobs,bg,bgerr,lumi = 3., 4.1, 0.6533758489567854, 35.9/fb
    # nobs,bg,bgerr,lumi = 0, 1., 0.2, 35.9/fb
    # nobs,bg,bgerr,lumi = 0, 0.001, 0.01, 35.9/fb
    nobs,bg,bgerr,lumi = 3905,3658.3,238.767, 35.9/fb
    dataset = SimpleSpeyDataSet ( nobs, bg, bgerr, lumi )
    computer = SpeyRetriever ( dataset, "1bin", 1. )
    ul = computer.getUpperLimit ( evaluationType = observed,
                                  limit_on_xsec = True )
    print ( "ul", ul )
    ule = computer.getUpperLimit ( evaluationType = apriori,
                                   limit_on_xsec = True )
    print ( "ule", ule )
