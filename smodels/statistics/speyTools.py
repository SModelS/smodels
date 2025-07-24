#!/usr/bin/env python3

"""
.. module:: speyTools
   :synopsis: a module that contains tools and convenience methods
              that we use in connection with spey.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

__all__ = [ "SpeyComputer", "SpeyAnalysesCombosComputer" ]

from typing import Union, Text, Tuple, Dict, List
import sys
from spey import ExpectationType, StatisticalModel, get_backend
import spey
try:
    from spey.system.exceptions import AsimovTestStatZero
except ImportError: # comes only with newer versions of spey
    AsimovTestStatZero = Exception # a dummy so we can still try
from smodels.base.smodelsLogging import logger
from smodels.base.physicsUnits import fb
from smodels.experiment.datasetObj import DataSet
from smodels.statistics.basicStats import observed, apriori, aposteriori, NllEvalType
from smodels.base.crossSection import XSection
import numpy as np

_debug = { "writePoint": False } # for debugging only

class SpeyComputer:
    """ the facade that delegates all statistical computations to spey.
    takes care of all interactions with spey except for analysis combinations. """

    __slots__ = [ "speyModels", "dataset", "backendType", "nsig",
                  "model_index" ]

    def __init__ ( self, dataset, backendType : str, nsig : Union[float,list],
                   deltas_rel : Union[None,float] = None ):
        """ initialise with dataset.
        :param dataset: a smodels (combined)dataset, or a onnx file (for now, will change)
        :param backendType: the type of backend to use ( "1bin", "SL", "pyhf", "ML" )
        :param nsig: signal yield, either as float or as list
        :param deltas_rel: relative error on signal. currently unused
        """
        if type(dataset) not in [ list] and dataset.getType() not in \
                [ "efficiencyMap", "combined" ]:
            logger.error(f"I do not recognize the dataset type {dataset.getType()}")

        self.backendType = backendType
        self.dataset = dataset
        self.nsig = nsig
        self.speyModels = self.getStatModels ( nsig )
        self.model_index = 0 # we might have several models and need to choose

    def getStatModels(self, nsig ) -> StatisticalModel:
        """ retrieve the statistical model """
        assert self.backendType in [ "1bin", "SL", "ML", "pyhf" ], f"unknown backend type {self.backendType}"
        if hasattr ( self.dataset.globalInfo, "onnxFile" ):
            # an onnxfile is defined, we use it!
            return self.getNNModel ( nsig )
        if self.backendType == "pyhf":
            return self.getStatModelsPyhf ( nsig )
        if self.dataset.getType() == "efficiencyMap":
            return self.getStatModelSingleBin ( nsig )
        return self.getStatModelMultiBin ( nsig )

    @classmethod
    def create ( cls, dataset : DataSet, xsection : XSection,
                 predictions, deltas_rel : float ):
        """ generic entry point for theoryPrediction. given the
            data, create the appropriate SpeyComputer

        :param dataset: the dataset to create the stats model for
        :param xsection: the cross section, needed just as a scale factor
        :param predictions: the individual predictions, in case of multiple SRs
            for e.g. upper limits on them
        """
        dataType = dataset.getType()
        if dataType == "upperLimit":
            from smodels.tools.runtime import experimentalFeatures
            if not experimentalFeatures():
                computer = 'N/A'
            else:
                computer = SpeyComputer.forTruncatedGaussian(self)
                if computer is None: # No evaluationType UL available
                    computer = 'N/A'

        elif dataType == "efficiencyMap":
            nsig = (xsection.value * dataset.getLumi()).asNumber()
            computer = SpeyComputer.forSingleBin(dataset=dataset,
                                                  nsig=nsig,deltas_rel=deltas_rel)

        elif dataType == "combined":
            ## FIXME this should be removable
            # Get dictionary with dataset IDs and signal yields
            srNsigDict = {pred.dataset.getID() :
                          (pred.xsection.value*pred.dataset.getLumi()).asNumber()
                          for pred in predictions}

            # Get ordered list of datasets:
            if hasattr(dataset.globalInfo, "covariance"):
                datasetList = dataset.globalInfo.datasetOrder[:]
                # Get list of signal yields corresponding to the dataset order:
                srNsigs = [srNsigDict[dataID] if dataID in srNsigDict else 0.0
                       for dataID in datasetList]
                # Get computer
                computer = SpeyComputer.forMultiBinSL(dataset=dataset,
                                        nsig=srNsigs,deltas_rel = deltas_rel)

            elif hasattr(dataset.globalInfo, "jsonFiles"):
                datasetList = [ds.getID() for ds in dataset.origdatasets]
                # Get list of signal yields corresponding to the dataset order:
                srNsigs = [srNsigDict[dataID] if dataID in srNsigDict else 0.0
                       for dataID in datasetList]
                # Get computer
                computer = SpeyComputer.forPyhf(dataset=dataset,
                                         nsig=srNsigs, deltas_rel = deltas_rel)
        return computer


    @classmethod
    def forSingleBin(cls, dataset, nsig, deltas_rel):
        """ get a speycomputer for an efficiency map (single bin).

        :param dataset: DataSet object
        :param nsig: Number of signal events for each SR
        :deltas_rel: Relative uncertainty for the signal

        :returns: a SpeyComputer
        """
        computer = SpeyComputer(dataset=dataset, backendType="1bin",
                                nsig=nsig, deltas_rel=deltas_rel)

        return computer

    @classmethod
    def forMultiBinSL(cls,dataset, nsig, deltas_rel):
        """ get a statscomputer for simplified likelihood combination.

        :param dataset: CombinedDataSet object
        :param nsig: Number of signal events for each SR
        :deltas_rel: Relative uncertainty for the signal

        :returns: a StatsComputer
        """

        computer = SpeyComputer(dataset=dataset,
                                backendType="SL",
                                 nsig=nsig, deltas_rel=deltas_rel)
        return computer

    def getNNModel ( self, nsig ):
        """ here is the code for how we create a spey speyModel that uses an NN
            as its backend """
        #For the moment, the ml-likelihood backned has to be dowanloaded and installed manually. In the future, the backend could be directly embedded and installed in SmodelS.
#        MLlikePath='/Users/humberto/Documents/work/learn_pyhf_smodels/ML_LHClikelihoods'
#        MLlikePath='/home/walten/git/ML_LHClikelihoods/'

#        sys.path.append(MLlikePath)
        import spey
        stat_wrapper = spey.get_backend('ml.likelihoods')

        #The current ML-likelihood backend takes a local path. This has to be
        #updated, I think the idea is that the backend takes 'dataset' as input,
#        where 'dataset' is a loaded onnx model. Thus, the onnx has to be loaded
#        by SmodelS beforehand. It seems to be the more consistent wrt the other
#        backends.

        #network_path='/Users/humberto/Documents/work/learn_pyhf_smodels/ML_LHClikelihoods/ML_models/ATLAS-SUSY-2018-04/ensemble_model.onnx'
        onnxBlob=self.dataset.globalInfo.onnx
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
        #self.speyModel = get_ml_model ( ... )
        return [ speyModel ]

    def getStatModelMultiBin(self,
        nsig: Union[np.ndarray, List[Dict[Text, List]], List[float]] ):
        """
        Create a statistical model from multibin data.

        :param nsig: number of signal events. For simplified likelihood backend
            this input can contain `np.array` or `List[float]` which contains
            signal yields per region.
            For `pyhf` backend this input evaluationType to be a JSON-patch
            i.e. `List[Dict]`, see `pyhf` documentation for details on
            JSON-patch format.
        :param delta_sys: systematic uncertainty on signal. Currently unused.
        :param allow_negative_signal: if True, the evaluationType upper limit on mu,
            used to find the best statistical model, can be negative.

        :returns: spey StatisticalModel object.

        """
        dataset = self.dataset
        obsN = [ x.dataInfo.observedN for x in dataset._datasets ]
        bg = [ x.dataInfo.expectedBG for x in dataset._datasets ]
        cov = dataset.globalInfo.covariance
        lumi = float ( dataset.getLumi().asNumber(1./fb) )
        thirdmomenta=[]
        for ds in dataset._datasets:
            if hasattr ( ds.dataInfo, "thirdMoment" ):
                thirdmomenta.append ( ds.dataInfo.thirdMoment )
        if len(thirdmomenta)==0: # SLv1
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
            return [ speyModel ]
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
        return [ speyModel ]

    def getStatModelsPyhf(self, nsig: Union[float, np.ndarray] ):
        """
        Create statistical model from a pyhf json file

        :param nsig: signal yields.
        :return: spey StatisticalModel object.
        """
        dataset = self.dataset
        stat_wrapper = get_backend("pyhf")
        from smodels.tools.speyPyhf import SpeyPyhfData
        data = SpeyPyhfData.createDataObject ( dataset, self.nsig )
        models = []
        patches = data.patchMaker()
        for i in range( len(data.inputJsons ) ):
            # idx = self.getBestCombinationIndex( data )
            inputJson = data.inputJsons[i]
            signal_patch = patches[i]
            #print ( "inputJsons", inputJson )
            # import IPython; IPython.embed( colors = "neutral" ); sys.exit()
            analysis = dataset.globalInfo.id

            speyModel = stat_wrapper( analysis = analysis,
                            signal_patch = signal_patch,
                            background_only_model = inputJson )
            models.append ( speyModel )
        self.speyModels = models
        self.model_index = self.getBestCombinationIndex( data )
        return models

    def getStatModelSingleBin(self, nsig: Union[float, np.ndarray],
            delta_sys : Union[None,float] = None,
            allow_negative_signal : bool = False
    ):
        """
        Create statistical model from a single bin or multiple uncorrelated regions.

        :param nsig: signal yields.
        :return: spey StatisticalModel object.

        :raises NotImplementedError: If requested backend has not been recognised.
        """
        dataset = self.dataset
        stat_wrapper = get_backend("default_pdf.uncorrelated_background")
        id = f"{dataset.globalInfo.id}:{dataset.dataInfo.dataId}"

        speyModel = stat_wrapper(
                        data = [float(dataset.dataInfo.observedN)],
                        background_yields = [float(dataset.dataInfo.expectedBG)],
                        absolute_uncertainties = [float(dataset.dataInfo.bgError)],
                        signal_yields = nsig,
                        xsection = float (nsig/dataset.getLumi().asNumber(1./fb)),
                        analysis = id,
#                        backend = 'simplified_likelihoods'
        )
        return [ speyModel ]

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
            return expectedDict[expected]
        logger.error( f'{expected} is not a valid expectation type. Possible expectation types are True (observed), False (apriori) and "posteriori".' )
        return None

    def get_five_values ( self, evaluationType : NllEvalType,
                      return_nll : bool = False,
                      check_for_maxima : bool = False )-> Dict:
        """ return the Five Values: l(bsm), l(sm), muhat, l(muhat), sigma(mu_hat)

        :param check_for_maxima: if true, then check lmax against l(sm) and l(bsm)
             correct, if necessary
        """
        ret = self.maximize_likelihood ( evaluationType = evaluationType, return_nll = return_nll  )
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


    def poi_upper_limit ( self, evaluationType : NllEvalType,
           limit_on_xsec : bool = False, model_index : Union [int,None] = None ) -> float:
        """ simple frontend to spey::poi_upper_limit

        :param limit_on_xsec: if True, then return the limit on the
            cross section
        :param model_index: if None, then get upper limit for most sensitive model,
            if integer, get UL for that model
        """
        if model_index == None:
            model_index = self.model_index
        exp = self.translateExpectationType ( evaluationType )

        try:
            ret = self.speyModels[model_index].poi_upper_limit ( evaluationType = exp )
        except ValueError as e:
            logger.warning ( f"when computing upper limit for SL: {e}. Will try with other method" )
            sys.exit(-1)
        except AsimovTestStatZero as e:
            logger.debug ( f"spey returned: {e}. will interpret as ul=inf" )
            ret = float("inf")
        ret = float(ret) # cast for the printers
        if limit_on_xsec:
            totsig = self.nsig
            if type ( self.nsig ) in [ list ]:
                totsig = sum ( self.nsig )
            xsec = totsig / self.dataset.globalInfo.lumi
            ret = ret * xsec
        return ret

    def getBestCombinationIndex(self, data ):
        """find the index of the best evaluationType combination"""
        if len(data.inputJsons) == 1:
            return 0
        logger.debug( f"Finding best evaluationType combination among {len(data.inputJsons)} workspace(s)" )
        ulMin = float("+inf")
        i_best = None
        for i_ws in range(len(data.inputJsons)):
            if data.totalYield == 0.:
                continue
            if data.zeroSignalsFlag[i_ws] == True:
                logger.debug( f"Workspace number {i_ws} has zero signals" )
                continue
            else:
                ul = self.poi_upper_limit(evaluationType=apriori, model_index=i_ws)
            if ul == None:
                continue
            if ul < ulMin:
                ulMin = ul
                i_best = i_ws
        return i_best


    def asimov_likelihood ( self, poi_test : float, evaluationType : NllEvalType,
                            return_nll : bool ) -> float:
        """ simple frontend to spey functionality """
        self.checkMinimumPoi( poi_test )
        evaluationType = self.translateExpectationType ( evaluationType )
        return self.speyModel[self.model_index].asimov_likelihood ( poi_test = poi_test,
            evaluationType = evaluationType, return_nll = return_nll )

    @classmethod
    def forPyhf(cls, dataset, nsig, deltas_rel):
        """ get a statscomputer for pyhf combination.

        :param dataset: CombinedDataSet object
        :param nsig: Number of signal events for each SR
        :deltas_rel: Relative uncertainty for the signal

        :returns: a StatsComputer
        """
        computer = SpeyComputer(dataset=dataset, backendType="pyhf",
                                nsig=nsig, deltas_rel=deltas_rel)
        return computer

    def checkMinimumPoi ( self, poi_test : float ):
        """ check if poi is below minimum_poi """
        config = self.speyModels[self.model_index].backend.config()
        if poi_test < config.minimum_poi:
            logger.error ( f'Calling likelihood for {self.dataset.globalInfo.id} (using combination of SRs) for a mu giving a negative total yield. mu = {mu} and minimum_mu = {config.minimum_poi}.' )

    def likelihood ( self, poi_test : float, evaluationType : NllEvalType,
                            return_nll : bool ) -> float:
        """ simple frontend to spey functionality """
        self.checkMinimumPoi ( poi_test )
        evaluationType = self.translateExpectationType ( evaluationType )
        ret = self.speyModels[self.model_index].likelihood ( poi_test = poi_test,
            evaluationType = evaluationType, return_nll = return_nll )
        return float(ret)

    def maximize_likelihood ( self, evaluationType : NllEvalType,
           allow_negative_signal : bool = True,
           return_nll : bool = False  ) -> Tuple[float,float]:
        """ simple frontend to spey functionality
        :param return_nll: if True, return negative log likelihood
        :param allow_negative_signal: allow also negative muhats
        :returns: tuple of muhat,lmax
        """
        evaluationType = self.translateExpectationType ( evaluationType )
        speyret = self.speyModels[self.model_index].maximize_likelihood ( evaluationType = evaluationType,
                allow_negative_signal = allow_negative_signal,
                return_nll = return_nll )
        ret = { "muhat": float(speyret[0]), "lmax": float(speyret[1]) }
        ## not clear if bounds will be hard bounds
        if not allow_negative_signal and speyret[0]< 0.:
            llhd0 = self.likelihood ( 0., evaluationType = evaluationType, return_nll = return_nll )
            ret = { "muhat": 0., "lmax": float(llhd0) }
        return ret

    def sigma_mu ( self, poi_test : float, evaluationType : NllEvalType, 
                   allow_negative_signal : bool = False ) -> float:
        """ determine sigma at poi_test.
        :param: FIXME allow_negative_signal should not be needed!
        """
        test_statistic = "q" if allow_negative_signal else "qmutilde"
        exp = SpeyComputer.translateExpectationType ( evaluationType )
        sigma_mu = self.speyModels[self.model_index].sigma_mu( poi_test=poi_test,evaluationType=exp,
                                            test_statistics=test_statistic )
        return float(sigma_mu)

    def maximize_asimov_likelihood ( self, evaluationType : NllEvalType,
           return_nll : bool = False ) -> Tuple[float,float]:
        """ simple frontend to spey functionality
        :param return_nll: if True, return negative log likelihood
        :param allow_negative_signal: allow also negative muhats
        :returns: tuple of muhat,lmax
        """
        evaluationType = self.translateExpectationType ( evaluationType )
        # init = self.getSpeyInitialisation ( True )
        # opt = init["optimiser"]
        # opt["test_statistics"]="qmutilde"
        ret = self.speyModels[self.model_index].maximize_asimov_likelihood ( evaluationType = evaluationType,
            return_nll = return_nll ) # , **opt )
        assert ret[0]>=0., "maximum of asimov likelihood should not be below zero"
        for k,v in ret.items():
            ret[k]=float(v)
        #if not allow_negative_signal and ret[0]< 0.:
        #    ret = ( 0., self.likelihood ( 0., evaluationType = evaluationType, return_nll = return_nll ) )
        return ret

    @property
    def xsection(self):
        return self.speyModels[self.model_index].xsection

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


class SpeyAnalysesCombosComputer:
    """ the facade class to spey for analysis combinations """
    def __init__  ( self, theorypreds, deltas_rel : float ):
        self.theorypreds = theorypreds
        self.deltas_rel = deltas_rel
        # self.pprint()
        from spey import UnCorrStatisticsCombiner
        models = []
        self.xsecs = []
        for pred in theorypreds:
            computer = SpeyComputer.create ( pred.dataset, pred.xsection, None,
                    deltas_rel )
            models.append ( computer.speyModels[0] )
            xsec = pred.xsection.value
            self.xsecs.append ( xsec )
        self.speyModel = UnCorrStatisticsCombiner ( *models )

    def pprint ( self ):
        """ print some info about the analyses combination """
        print ( "[SpeyAnalysesCombo] combining:" )
        for tp in self.theorypreds:
            print ( f"    {tp.analysisId()}:{tp.dataId()}" )

    def poi_upper_limit ( self, evaluationType : NllEvalType,
           limit_on_xsec : bool = False ) -> float:
        """ simple frontend, to spey::poi_upper_limit

        :param limit_on_xsec: if True, then return the limit on the
                              cross section
        """
        exp = SpeyComputer.translateExpectationType ( evaluationType )
        ret = self.speyModel.poi_upper_limit ( evaluationType = exp )
        ret = float(ret) # cast for the printers
        if limit_on_xsec:
            totxsec = sum(self.xsecs,0.*fb)
            ret = ret * totxsec
        return ret

    def likelihood ( self, poi_test : float, evaluationType : NllEvalType,
                            return_nll : bool ) -> float:
        """ simple frontend to spey functionality """
        evaluationType = SpeyComputer.translateExpectationType ( evaluationType )
        ret = self.speyModel.likelihood ( poi_test = poi_test,
            evaluationType = evaluationType, return_nll = return_nll )
        return float(ret)

    def maximize_likelihood ( self, evaluationType : NllEvalType,
           allow_negative_signal : bool = True,
           return_nll : bool = False  ) -> Tuple[float,float]:
        """ simple frontend to spey functionality
        :param return_nll: if True, return negative log likelihood
        :param allow_negative_signal: allow also negative muhats
        :returns: tuple of muhat,lmax
        """
        evaluationType = SpeyComputer.translateExpectationType ( evaluationType )
        speyret = self.speyModel.maximize_likelihood ( evaluationType = evaluationType,
                allow_negative_signal = allow_negative_signal,
                return_nll = return_nll )
        ret = { "muhat": float(speyret[0]), "lmax": float(speyret[1]) }
        ## not clear if bounds will be hard bounds
        if not allow_negative_signal and speyret[0]< 0.:
            l0 = self.likelihood ( 0., evaluationType = evaluationType, return_nll = return_nll )
            ret = { "muhat": 0., "lmax": float(l0) }
        return ret

    def get_five_values ( self, evaluationType : NllEvalType,
                      return_nll : bool = False,
                      check_for_maxima : bool = False )-> Dict:
        """ method returning the Five Values:
        l(bsm), l(sm), muhat, l(muhat), sigma(mu_hat)

        :param check_for_maxima: if true, then check lmax against l(sm) and l(bsm)
            correct, if necessary
        """
        ret = self.maximize_likelihood ( evaluationType = evaluationType, return_nll = return_nll  )
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

# print ( "five values", ret, [ type(v) for k,v in ret.items() ] )
        return ret

    def sigma_mu ( self, poi_test : float, evaluationType : NllEvalType, 
                   allow_negative_signal : bool = False ) -> float:
        """ determine sigma at poi_test.
        :param: FIXME allow_negative_signal should not be needed!
        """
        test_statistic = "q" if allow_negative_signal else "qmutilde"
        exp = SpeyComputer.translateExpectationType ( evaluationType )
        sigma_mu = self.speyModel.sigma_mu( poi_test=poi_test,evaluationType=exp,
                                            test_statistics=test_statistic )
        return float(sigma_mu)


if __name__ == "__main__":
    # nobs,bg,bgerr,lumi = 3., 4.1, 0.6533758489567854, 35.9/fb
    # nobs,bg,bgerr,lumi = 0, 1., 0.2, 35.9/fb
    # nobs,bg,bgerr,lumi = 0, 0.001, 0.01, 35.9/fb
    nobs,bg,bgerr,lumi = 3905,3658.3,238.767, 35.9/fb
    dataset = SimpleSpeyDataSet ( nobs, bg, bgerr, lumi )
    computer = SpeyComputer ( dataset, "1bin", 1. )
    ul = computer.poi_upper_limit ( evaluationType = observed, limit_on_xsec = True )
    print ( "ul", ul )
    ule = computer.poi_upper_limit ( evaluationType = apriori, limit_on_xsec = True )
    print ( "ule", ule )
