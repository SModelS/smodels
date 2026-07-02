#!/usr/bin/env python3

"""
.. module:: speyTools
   :synopsis: a module that contains tools and convenience methods
              that we use in connection with spey.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

__all__ = [ "SpeyRetriever"]

from typing import Optional
from spey import ExpectationType, get_backend
from smodels.base.physicsUnits import UnitXSec
from smodels.statistics.exceptions import SModelSStatisticsError as SModelSError
from smodels.experiment.datasetObj import DataSet, CombinedDataSet
# spey.set_optimiser( "iminuit" )

from smodels.base.smodelsLogging import logger
from smodels.base.physicsUnits import fb, UnitLumi
from smodels.statistics.basicStats import observed, apriori, aposteriori, NllEvalType
from smodels.statistics.analysesCombinations import AnaCombLikelihoodComputer

_debug = { "writePoint": False } # for debugging only


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

        # optimiser_arguments = None
        # optimiser_arguments = { "method": "SLSQP" }
        optimiser_arguments = { "method": "L-BFGS-B" }
        optimiser_arguments.update ( kwargs )
        ret = self.speyModel.poi_upper_limit ( expected = exp,
               expected_pvalue = expected_pvalue,
               optimiser_arguments = optimiser_arguments )
        if nSigma == 0:
            ret = float ( ret )
        elif nSigma == 1:
            ret = float(ret[0])
        elif nSigma == -1:
            ret = float(ret[-1])
        return ret


class SpeyRetriever:
    """ simple class that retrieves and constructs the sub computers
    using the Spey interface."""
    

    @classmethod
    def forMultiBinSL(cls, srSet: str, dataset : CombinedDataSet, nsigDict : dict,
            deltas_rel : Optional[float] = 0.0 ) -> SpeyModelFacade:
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

        import spey
        covs = dataset.globalInfo.cachedModels
        lumi = float ( dataset.getLumi().asNumber(1./fb) )
        thirdmomenta=[]
        for ds in dataset._datasets:
            if hasattr ( ds.dataInfo, "thirdMoment" ):
                thirdmomenta.append ( ds.dataInfo.thirdMoment )
        type_n_models = dataset.globalInfo.statModels[srSet]
        type_n_model = type_n_models[0] # get first one
        mtype = type_n_model[0]
        assert mtype == "sl", f"expected sl but got {mtype} for type of stats model"
        covname = type_n_model[1]
        nsig = list(nsigDict.values())
        cov = covs[covname]
        offset=0
        assert type(cov)==list, f"covariance field has wrong type: {type(cov)} in {dataset.globalInfo.id}"
        assert len(cov)>0, f"covariance matrix has length {len(cov)}."

        n = len(cov)
        obsN = [ x.dataInfo.observedN for x in dataset._datasets[offset:offset+n] ]
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
            except spey.PluginError: ## older spey?
                stat_wrapper = get_backend("default_pdf.correlated_background")
            if _debug["writePoint"]:
                obsN = [ x.dataInfo.observedN for x in dataset._datasets[:n] ]
                f=open ( "data.txt","wt" )
                f.write ( f"obsN={obsN}\n" )
                f.write ( f"bg={bg}\n" )
                f.write ( f"cov={cov}\n" )
                f.write ( f"nsig={nsig}\n" )
                f.write ( f"analysis='{dataset.globalInfo.id}'\n" )
                f.write ( f"lumi={lumi}\n" )
                f.close()

            speyModel = stat_wrapper( data = obsN,
                            background_yields = bg, covariance_matrix = cov,
                            signal_yields = nsig,
            #                xsection = [ x / lumi for x in nsig ],
                            analysis = dataset.globalInfo.id,
            )
            facade = SpeyModelFacade ( speyModel, "SL", covname, xsec )
            return facade
        # SLv2
        try:
            stat_wrapper = get_backend("default.third_moment_expansion")
        except ImportError:
            stat_wrapper = get_backend("default_pdf.third_moment_expansion")
        speyModel = stat_wrapper( data = obsN,
                        background_yields = bg, covariance_matrix = cov,
                        signal_yields = nsig,
        #                xsection = [ x / lumi for x in nsig ],
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
        return facade

    @classmethod
    def forSingleBin( cls, srSet: str, dataset, nsigDict, deltas_rel : float = 0.2,
                      lumi : Optional[UnitLumi]=None ) -> SpeyModelFacade:
        """ get a sub computer for an efficiency map (single bin).

        :param dataset: DataSet object
        :param nsigDict: Dictionary with signal yields for each SR
        :deltas_rel: Relative uncertainty for the signal

        :returns: a sub computer
        :raises NotImplementedError: If requested backend has not been recognised.
        """

        import spey

        try:
            stat_wrapper = get_backend("default.uncorrelated_background")
        except spey.PluginError: ## older spey?
            stat_wrapper = get_backend("default_pdf.uncorrelated_background")
        id = f"{dataset.globalInfo.id}:{dataset.dataInfo.dataId}"
        nsig = nsigDict.get(srSet, 0)
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
        return facade

    @classmethod
    def forNNs(cls, srSet: str, dataset, nsigDict ) -> SpeyModelFacade:
        """ get a sub computer for an NN combination.

        :param dataset: CombinedDataSet object
        :param nsigDict: Dictionary with signal yields for each SR
        :deltas_rel: Relative uncertainty for the signal

        :returns: a sub computer
        """

        logger.error ("speyTools backend to NN is not implemented yet!" )
        raise SModelSError ("speyTools backend to NN is not implemented yet!" )
    
        globalInfo = dataset.globalInfo
        labelToONNX = {}
        srMappings = globalInfo.srMappings

        for sr in globalInfo.regionSets[srSet]:
            if sr not in srMappings:
                logger.error ( f"SR {sr} defined in srSet {srSet} not found in srMappings for dataset {dataset.globalInfo.id}" )
                raise SModelSError ( f"SR {sr} defined in srSet {srSet} not found in srMappings for dataset {dataset.globalInfo.id}" )
            labelToONNX[sr] = srMappings[sr]["onnx"]

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

        onnxBlob=dataset.globalInfo.cachedModels[model_filename]
        # self.speyModel = stat_wrapper(nsig,onnxBlob) # this is how i want it long run
        ## the following code is just for now to see if it works in principle
        import tempfile
        tempf = tempfile.mktemp ( prefix="/tmp/", suffix=".onnx" )
        # tempf = "/tmp/my.onnx"
        f = open ( tempf, "wb" )
        import onnx
        onnx.save ( onnxBlob, f )
        f.close()
        xsec = sum(list(nsigDict.values())) / dataset.globalInfo.lumi
        stat_wrapper = get_backend("onnx")
        speyModel = stat_wrapper(nsigDict,tempf)
        facade = SpeyModelFacade ( speyModel, "nn",
                                    model_filename, xsec )
        # os.unlink ( tempf )
        return facade

    @classmethod
    def forPyhf(cls, srSet: str, dataset, nsigDict) -> SpeyModelFacade:
        """ get a sub computer for pyhf combination.

        :param dataset: CombinedDataSet object
        :param nsigDict: Dictionary with signal yields for each SR

        :returns: a sub computer
        """
        stat_wrapper = get_backend("pyhf")
        from smodels.statistics.speyPyhf import SpeyPyhfData
        data = SpeyPyhfData.createDataObject ( dataset, nsigDict, srSet )
        models = []
        patches = data.patchMaker()
        # for i in range( len(data.inputJsons ) ):
            # idx, _ = self.getBestCombinationIndex( data )
        inputJson = data.inputJson # s[i]
        signal_patch = patches # [i]
        #print ( "inputJsons", inputJson )
        # import IPython; IPython.embed( colors = "neutral" ); sys.exit()
        analysis = dataset.globalInfo.id

        speyModel = stat_wrapper( analysis = analysis,
                        signal_patch = signal_patch,
                        background_only_model = inputJson )
        xsec = sum (list(nsigDict.values())) / dataset.globalInfo.lumi
        mname = data.jsonFile# s[i]
        facade = SpeyModelFacade ( speyModel, "pyhf", mname, xsec )
        return facade

    @classmethod
    def forTruncatedGaussian(cls, theorypred, corr : float =0.6 ) -> None:
        """ get a sub computer for truncated gaussians
        :param theorypred: TheoryPrediction object
        :param corr: correction factor: \
                ULexp_mod = ULexp / (1. - corr*((ULobs-ULexp)/(ULobs+ULexp))) \
                a factor of corr = 0.6 is proposed.
        :returns: list of subComputers (with a single entry)
        """
        logger.error ( "speyTools no truncated Gaussian backend exists" )
        return None

    @classmethod
    def forAnalysesComb(cls,theoryPredictions, deltas_rel : float = 0.0) -> AnaCombLikelihoodComputer:
        """ get a sub computer for combination of analyses
        :param theoryPredictions: list of TheoryPrediction objects
        :param deltas_rel: relative error for the signal
        :returns: a sub computer
        """

        computer = AnaCombLikelihoodComputer( theoryPredictions=theoryPredictions,
                                              deltas_rel=deltas_rel )
        return computer

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
