#!/usr/bin/env python3

"""
.. module:: speyTools
   :synopsis: a module that contains tools and convenience methods 
              that we use in connection with spey.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

__all__ = [ "SpeyComputer" ]

from typing import Union, Text, Tuple, Dict, List
from spey import ExpectationType, StatisticalModel, get_backend
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import fb
import numpy as np

class SpeyComputer:
    __slots__ = [ "statModel", "dataset", "weight", "backendType", "nsig" ]

    def __init__ ( self, dataset, backendType : str, nsig : Union[float,list], 
                   deltas_rel : Union[None,float] = None ):
        """ initialise with dataset.
        :param dataset: a smodels (combined)dataset, or a onnx file (for now, will change)
        :param backendType: the type of backend to use ( "1bin", "SL", "pyhf", "ML" )
        :param nsig: signal yield, either as float or as list
        :param deltas_rel: relative error on signal. currently unused
        """
        if type(dataset) not in [ list] and dataset.getType() not in [ "efficiencyMap", "combined" ]:
            logger.error ( f"I do not recognize the dataset type {dataset.getType()}" )

        self.backendType = backendType
        self.dataset = dataset
        self.nsig = nsig
        self.statModel = self.getStatModel ( nsig )

    def getStatModel(self, nsig ) -> StatisticalModel:
        """ retrieve the statistical model """
        if type(self.dataset)==list: # ok, we have a combined dataset!
            return self.getAnalysisCombinationModel ( nsig )
        if hasattr ( self.dataset.globalInfo, "onnxFile" ):
            # an onnxfile is defined, we use it!  
            return self.getNNModel ( nsig ) 
        if self.backendType == "pyhf":
            return self.getStatModelPyhf ( nsig )
        if self.dataset.getType() == "efficiencyMap":
            return self.getStatModelSingleBin ( nsig )
        return self.getStatModelMultiBin ( nsig )

    def getAnalysisCombinationModel ( self, nsig ):
        computer = SpeyComputer(dataset=dataset, backendType="combo", 
                                nsig=nsig, deltas_rel=deltas_rel)
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
        """ here is the code for how we create a spey statModel that uses an NN 
            as its backend """
        #from spey import get_ml_model # import the spey method for it
               
        #For the moment, the ml-likelihood backned has to be dowanloaded and installed manually. In the future, the backend could be directly embedded and installed in SmodelS.
        import sys
        MLlikePath='/Users/humberto/Documents/work/learn_pyhf_smodels/ML_LHClikelihoods'
        sys.path.append(MLlikePath)
        import spey
        stat_wrapper = spey.get_backend('ml.likelihoods')
        
        #The current ML-likelihood backend takes a local path. This has to be updated, I think the idea is that the backend takes 'dataset' as input, where 'dataset' is a loaded onnx model. Thus, the onnx has to be loaded by SmodelS beforehand. It seems to be the more consistent wrt the other backends.
        
#network_path='/Users/humberto/Documents/work/learn_pyhf_smodels/ML_LHClikelihoods/ML_models/ATLAS-SUSY-2018-04/ensemble_model.onnx'
        network_path=self.dataset.globalInfo.onnxFile
        self.statModel = stat_wrapper(nsig,network_path)
        
        #self.statModel = get_ml_model ( ... )
        return self.statModel

    def getStatModelMultiBin(self,
        nsig: Union[np.ndarray, List[Dict[Text, List]], List[float]] ):
        """
        Create a statistical model from multibin data.

        :param nsig: number of signal events. For simplified likelihood backend this input can
                       contain `np.array` or `List[float]` which contains signal yields per region.
                       For `pyhf` backend this input expected to be a JSON-patch i.e. `List[Dict]`,
                       see `pyhf` documentation for details on JSON-patch format.
        :param delta_sys: systematic uncertainty on signal. Only used for simplified likelihood backend.
        :param allow_negative_signal: if True, the expected upper limit on mu, used to find the best statistical model, can be negative.

        :return: spey StatisticalModel object.

        :raises NotImplementedError: if input patter does not match to any backend specific input option.
        """
        dataset = self.dataset
        stat_wrapper = get_backend("default_pdf.correlated_background")
        obsN = [ x.dataInfo.observedN for x in dataset._datasets ]
        bg = [ x.dataInfo.expectedBG for x in dataset._datasets ]
        cov = dataset.globalInfo.covariance
        lumi = float ( dataset.getLumi().asNumber(1./fb) )
        # print ( "input", len(obsN), len(bg), len(cov), len(nsig) )
        self.statModel = stat_wrapper( data = obsN,
                        background_yields = bg, covariance_matrix = cov,
                        signal_yields = nsig,
                        xsection = [ x / lumi for x in nsig ],
                        analysis = dataset.globalInfo.id,
#                        backend = 'simplified_likelihoods'
        )
        return self.statModel

    def getStatModelPyhf(self, nsig: Union[float, np.ndarray] ):
        """
        Create statistical model from a single bin or multiple uncorrelated regions.

        :param nsig: signal yields.
        :return: spey StatisticalModel object.

        :raises NotImplementedError: If requested backend has not been recognised.
        """
        dataset = self.dataset
        stat_wrapper = get_backend("pyhf")
        from smodels.tools.speyPyhf import SpeyPyhfData
        data = SpeyPyhfData.createDataObject ( dataset, self.nsig )
        # import IPython; IPython.embed( colors = "neutral" ); sys.exit()
        analysis = dataset.globalInfo.id

        self.statModel = stat_wrapper( analysis = analysis,
                        signal_patch = data.patchMaker(),
                        background_only_model = data.inputJsons )
        return self.statModel

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

        self.statModel = stat_wrapper(
                        data = [float(dataset.dataInfo.observedN)],
                        background_yields = [float(dataset.dataInfo.expectedBG)],
                        absolute_uncertainties = [float(dataset.dataInfo.bgError)],
                        signal_yields = nsig,
                        xsection = float (nsig/dataset.getLumi().asNumber(1./fb)),
                        analysis = dataset.globalInfo.id,
#                        backend = 'simplified_likelihoods'
        )
        return self.statModel

    def translateExpectationType ( self, expected : Union [ bool, Text ] ) -> ExpectationType:
        """ translate the specification for expected values from smodels
            lingo to spey convention """
        if type(expected)==ExpectationType:
            return expected
        expectedDict = {False: ExpectationType.observed,
                        True: ExpectationType.apriori,
                        "apriori": ExpectationType.apriori,
                        "posteriori": ExpectationType.aposteriori}
        if expected in expectedDict:
            return expectedDict[expected]
        logger.error('%s is not a valid expectation type. Possible expectation types are True (observed), False (apriori) and "posteriori".' %expected)
        return None

    def get_five_values ( self, expected : Union [ bool, Text ],
                      return_nll : bool = False,
                      check_for_maxima : bool = False )-> Dict:
        """ return the Five Values: l(bsm), l(sm), muhat, l(muhat), sigma(mu_hat) 
        :param check_for_maxima: if true, then check lmax against l(sm) and l(bsm)
             correct, if necessary
        """
        ret = self.maximize_likelihood ( expected = expected, return_nll = return_nll  )
        lmax = ret['lmax']

        lbsm = self.likelihood ( poi_test = 1., expected=expected, return_nll = return_nll )
        ret["lbsm"] = lbsm
        lsm = self.likelihood ( poi_test = 0., expected=expected, return_nll = return_nll )
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


    def maximize_likelihood_timothee ( self, expected : Union[bool,Text],
            allow_negative_signal : bool = True,
            return_nll : bool = False  ) -> Tuple[float,float]:
        """ maximize likelihood, timothee style. i.e. if expected is
            posteriori then maximize asimov likelihood but with expected="apriori"
            will worry later about the correctness of this. """
        if expected == "posteriori":
            return self.maximize_asimov_likelihood(expected="apriori",
#                   allow_negative_signal = allow_negative_signal, # 
                   return_nll = return_nll )
        return self.maximize_likelihood( allow_negative_signal = allow_negative_signal,
               return_nll = return_nll, expected = expected )

    def poi_upper_limit ( self, expected : Union [ bool, Text ],
           limit_on_xsec : bool = False ) -> float:
        """ simple frontend, to spey::poi_upper_limit 
        :param limit_on_xsec: if True, then return the limit on the
                              cross section
        """
        from spey import ExpectationType
        exp = ExpectationType.aposteriori
        if expected == False:
            exp = ExpectationType.observed
        if expected in [ True, "prior" ]:
            exp = ExpectationType.apriori

        try:
            ret = self.statModel.poi_upper_limit ( expected = exp )
        except ValueError as e:
            logger.warning ( f"when computing upper limit for SL: {e}. Will try with other method" )
            sys.exit(-1)
        if limit_on_xsec:
            xsec = sum ( self.nsig ) / self.dataset.globalInfo.lumi    
            ret = ret * xsec
        return ret

    def _getBestStatModel(self, nsig, allow_negative_signal=False):
        """
        find the index of the best expected combination.

        :param nsig: list of signal yields.
        :param allow_negative_signal: if True, the expected upper limit on mu, used to find the best statistical model, can be negative.

        :return: the spey StatisticalModel object that computed the minimal apriori-expected poi upper limit.
        """
        dataset = self.dataset
        # Get the list of the names of the signal regions used in the json files
        listOfSRInJson=[]
        for SRnames in dataset.globalInfo.jsonFiles.values():
            listOfSRInJson += SRnames

        patches, listOfSignals = dataset._getPatches(nsig)
        mu_ul_exp_min = np.inf
        # Find best combination of signal regions
        for index, (patch, json) in enumerate(zip(patches,dataset.globalInfo.jsons)):
            # If the expected signal is 0 for each SR in the combined set of SRs, skip
            if all([sig==0. for sig in listOfSignals[index]]):
                continue
            # The x-section is at the level of the TheoryPrediction
            # if there are multiple sets of SRs, set a xsec_UL for the whole analysis, i.e. that uses all the SRs,
            # so that the resulting R value Is for the whole analysis
            xsec = float ( sum(nsig)/dataset.getLumi().asNumber(1./fb ) )
            from spey import get_correlated_nbin_statistical_model
            # It is possible to do differently and to set a xsec_UL on each set of SRs but that is not how it done in SModelS so far
            # xsec = sum(listOfSignals[index])/self.getLumi()
            statModel = get_correlated_nbin_statistical_model(analysis=dataset.globalInfo.id,
                                                            signal_yields=patch,
                                                            data=json,
                                                            xsection=xsec
                                                            )
            self.statModel = statModel
            # If all the SRs are used in the json files and there is only one json files, there is only one statModel.
            # No need to compute mu_ul_exp if not needed.
            if all([ds.dataInfo.dataId in listOfSRInJson for ds in dataset._datasets]) and len(dataset.globalInfo.jsons) == 1:
                return statModel

            config = statModel.backend.config()
            inits = self.getSpeyInitialisation ( True )
            lim = inits["limits"]
            bounds = config.bounds()
            try:
                mu_ul_exp = statModel.poi_upper_limit(
                        expected=ExpectationType.apriori, 
                        low_init = lim["low_ini"], hig_init = lim["hig_init"], 
                        maxiter = lim["maxiter"], 
                        optimiser_arguments = inits["optimiser"] )
            except ValueError as e:
                # if we dont get an answer, might just be this super region
                # is not a good choice
                logger.warn ( f"when trying to find best super region: {e}. will skip this one." )
                continue
            while abs(mu_ul_exp - bounds[config.poi_index][1]) / ( mu_ul_exp + bounds[config.poi_index][1]) <= 0.1:
                logger.debug('Expected upper limit on poi reached the upper bound. Will try again after increasing the upper bound.')
                bounds[config.poi_index] = (bounds[config.poi_index][1], bounds[config.poi_index][1]*10)
                mu_ul_exp = statModel.poi_upper_limit(
                        expected=ExpectationType.apriori, 
                        low_init = lim["low_ini"], hig_init = lim["hig_init"], 
                        maxiter = lim["maxiter"], 
                        optimiser_arguments = inits["optimiser"] )

            if mu_ul_exp == 0 and not allow_negative_signal:
                logger.warning(f'Expected upper limit on poi is negative when searching for the best statistical model for analysis {self.globalInfo.id}.')
            if mu_ul_exp == None:
                continue
            elif mu_ul_exp < mu_ul_exp_min:
                mu_ul_exp_min = mu_ul_exp
                bestStatModel = statModel

        # Check if a non-combined (uncorrelated) signal region is more contraining than the best combination obtained above
        # Check if a signal region is not in the list of SR names used in the json files
        for sig,ds in zip(nsig,self.dataset._datasets):
            dI = ds.dataInfo
            if dI.dataId not in listOfSRInJson:
                xsec = sig/self.dataset.getLumi()
                # Don't bother to compute eUL again (one could do it again if needed)
                if xsec.asNumber(fb) == 0.:
                    ## we have no values
                    continue
                mu_ul_exp = dI.expectedUpperLimit/xsec
                if mu_ul_exp < mu_ul_exp_min:
                    logger.info("Best constraining model is a single uncorrelated model.")
                    mu_ul_exp_min = mu_ul_exp
                    bestStatModel = ds.getStatModel(sig)

        if mu_ul_exp_min == np.inf:
            logger.error(f'No minimal upper limit on POI found for {self.globalInfo.id}')
            return None
        return bestStatModel

    def asimov_likelihood ( self, poi_test : float, expected : Union[bool,Text], 
                            return_nll : bool ) -> float:
        """ simple frontend to spey functionality """
        self.checkMinimumPoi( poi_test )
        expected = self.translateExpectationType ( expected )
        return self.statModel.asimov_likelihood ( poi_test = poi_test,
            expected = expected, return_nll = return_nll )

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

    @classmethod
    def forAnalysesComb(cls,theoryPredictions, deltas_rel):
        """ get a statscomputer for combination of analyses
        :param theoryPredictions: list of TheoryPrediction objects
        :param deltas_rel: relative error for the signal
        :returns: a StatsComputer
        """
        
        # Only allow negative signal if all theory predictions allow for it
        #allow_negative_signal = all([tp.statsComputer.allow_negative_signal
        #                            for tp in theoryPredictions])

        computer = SpeyComputer(dataset=theoryPredictions,  
                                backendType="analysesComb", 
                                nsig=None, deltas_rel=deltas_rel)
        
        return computer

    def checkMinimumPoi ( self, poi_test : float ):
        """ check if poi is below minimum_poi """
        config = self.statModel.backend.config()
        if poi_test < config.minimum_poi:
            logger.error ( f'Calling likelihood for {self.dataset.globalInfo.id} (using combination of SRs) for a mu giving a negative total yield. mu = {mu} and minimum_mu = {config.minimum_poi}.' )

    def likelihood ( self, poi_test : float, expected : Union[bool,Text], 
                            return_nll : bool ) -> float:
        """ simple frontend to spey functionality """
        self.checkMinimumPoi ( poi_test )
        expected = self.translateExpectationType ( expected )
        # init = self.getSpeyInitialisation ( True )
        return self.statModel.likelihood ( poi_test = poi_test,
            expected = expected, return_nll = return_nll )

    def maximize_likelihood ( self, expected : Union[bool,Text],
           allow_negative_signal : bool = True,
           return_nll : bool = False  ) -> Tuple[float,float]:
        """ simple frontend to spey functionality
        :param return_nll: if True, return negative log likelihood
        :param allow_negative_signal: allow also negative muhats
        :returns: tuple of muhat,lmax
        """
        expected = self.translateExpectationType ( expected )
        speyret = self.statModel.maximize_likelihood ( expected = expected, 
                allow_negative_signal = allow_negative_signal,
                return_nll = return_nll )
        ret = { "muhat": speyret[0], "lmax": speyret[1] }
        ## not clear if bounds will be hard bounds
        if not allow_negative_signal and speyret[0]< 0.:
            ret = { "muhat": 0., "lmax": self.likelihood ( 0., expected = expected, return_nll = return_nll ) }
        return ret

    def sigma_mu ( self, poi_test : float, expected : Union[bool,Text], allow_negative_signal : bool = False ):
        """ determine sigma at poi_test.
        :param: FIXME allow_negative_signal should not be needed!
        """
        test_statistic = "q" if allow_negative_signal else "qmutilde"
        sigma_mu = self.statModel.sigma_mu( poi_test=poi_test,expected=expected,
                                            test_statistics=test_statistic )
        return sigma_mu

    def maximize_asimov_likelihood ( self, expected : Union[bool,Text],
           return_nll : bool = False ) -> Tuple[float,float]:
        """ simple frontend to spey functionality 
        :param return_nll: if True, return negative log likelihood
        :param allow_negative_signal: allow also negative muhats
        :returns: tuple of muhat,lmax
        """
        expected = self.translateExpectationType ( expected )
        # init = self.getSpeyInitialisation ( True )
        # opt = init["optimiser"]
        # opt["test_statistics"]="qmutilde"
        ret = self.statModel.maximize_asimov_likelihood ( expected = expected, 
            return_nll = return_nll ) # , **opt )
        assert ret[0]>=0., "maximum of asimov likelihood should not be below zero"
        #if not allow_negative_signal and ret[0]< 0.:
        #    ret = ( 0., self.likelihood ( 0., expected = expected, return_nll = return_nll ) )
        return ret

    @property
    def xsection(self):
        return self.statModel.xsection

class SimpleSpeyDataSet:
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
            self.id = "SimpleSpeyDataSet"
            self.lumi = lumi

    def __init__ ( self, observedN : float, expectedBG : float,
                   bgError : float, lumi : fb ):
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
    dataset = SimpleSpeyDataSet ( nobs, bg, bgerr, lumi )
    computer = SpeyComputer ( dataset, 1. )
    ul = computer.poi_upper_limit ( expected = False, limit_on_xsec = True )
    print ( "ul", ul )
    ule = computer.poi_upper_limit ( expected = True, limit_on_xsec = True )
    print ( "ule", ule )
