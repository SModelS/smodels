#!/usr/bin/env python3

"""
.. module:: speyTools
   :synopsis: a module that contains tools and convenience methods 
              that we use in connection with spey.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

__all__ = [ "SpeyComputer" ]

from typing import Union, Text, Tuple, Dict, List
from spey import ExpectationType, StatisticalModel
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import fb
import numpy as np

class SpeyComputer:
    __slots__ = [ "statModel", "dataset" ]

    def __init__ ( self, dataset, nsig : Union[float,list], 
                   deltas_rel : Union[None,float] = None ):
        """ initialise with dataset.
        :param dataset: a smodels (combined)dataset
        :param nsig: signal yield, either as float or as list
        :param deltas_rel: relative error on signal. currently unused
        """
        if deltas_rel != None:
            logger.warning("Relative uncertainty on signal not supported by spey for a single region.")
        if dataset.getType() not in [ "efficiencyMap", "combined" ]:
            logger.error ( f"I do not recognize the dataset type {dataset.getType()}" )

        self.dataset = dataset
        self.statModel = self.getStatModel ( nsig )

    def getStatModel(self, nsig ) -> StatisticalModel:
        """ retrieve the statistical model """
        if self.dataset.getType() == "efficiencyMap":
            return self.getStatModelSingleBin ( nsig )
        return self.getStatModelMultiBin ( nsig )

    def getStatModelMultiBin(self,
        nsig: Union[np.ndarray, List[Dict[Text, List]], List[float]],
        delta_sys: float = 0.0,
        allow_negative_signal = False,

    ):
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

        if dataset.type == "simplified":
            nobs = [x.dataInfo.observedN for x in dataset.origdatasets]
            cov = dataset.globalInfo.covariance
            if type(cov) != list:
                raise SModelSError("covariance field has wrong type: %s" % type(cov))
            if len(cov) < 1:
                raise SModelSError("covariance matrix has length %d." % len(cov))
            bg = [x.dataInfo.expectedBG for x in dataset.origdatasets]
            third_moment = dataset.globalInfo.third_moment if hasattr(dataset.globalInfo, "third_moment") else None
            xsec = float ( sum(nsig)/dataset.getLumi().asNumber(1./fb) )
            from spey import get_correlated_nbin_statistical_model

            self.statModel = get_correlated_nbin_statistical_model(analysis = dataset.globalInfo.id,
                                                                signal_yields = nsig,
                                                                data = nobs,
                                                                covariance_matrix = cov,
                                                                backgrounds = bg,
                                                                third_moment = third_moment,
                                                                delta_sys = delta_sys,
                                                                xsection = xsec
                                                                )
        elif dataset.type == "pyhf":
            self.statModel = dataset._getBestStatModel(nsig, allow_negative_signal=allow_negative_signal)
        else:
            logger.error(f'Dataset of type "{dataset.type}" for analysis {dataset.globalInfo.id} is not of type "simplified" or "pyhf".')

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
        from spey import get_uncorrelated_nbin_statistical_model

        self.statModel = get_uncorrelated_nbin_statistical_model(
                            data = float(dataset.dataInfo.observedN),
                            backgrounds = float(dataset.dataInfo.expectedBG),
                            background_uncertainty = float(dataset.dataInfo.bgError),
                            signal_yields = nsig,
                            xsection = float (nsig/dataset.getLumi().asNumber(1./fb)),
                            analysis = dataset.globalInfo.id,
                            backend = 'simplified_likelihoods'
        )
        return self.statModel

    def translateExpectationType ( self, expected : Union [ bool, Text ] ) -> ExpectationType:
        """ translate the specification for expected values from smodels
            lingo to spey convention """
        expectedDict = {False:ExpectationType.observed,
                        True:ExpectationType.apriori,
                        "apriori": ExpectationType.apriori,
                        "posteriori":ExpectationType.aposteriori}
        if expected in expectedDict:
            return expectedDict[expected]
        logger.error('%s is not a valid expectation type. Possible expectation types are True (observed), False (apriori) and "posteriori".' %expected)
        return None

    def maximize_likelihood_timothee ( self, expected : Union[bool,Text],
            allowNegativeSignals : bool = True,
            return_nll : bool = False  ) -> Tuple[float,float]:
        """ maximize likelihood, timothee style. i.e. if expected is
            posteriori then maximize asimov likelihood but with expected="apriori"
            will worry later about the correctness of this. """
        if expected == "posteriori":
            return self.maximize_asimov_likelihood(expected="apriori",
#                   allowNegativeSignals = allowNegativeSignals, # 
                   return_nll = return_nll )
        return self.maximize_likelihood( allowNegativeSignals = allowNegativeSignals,
               return_nll = return_nll, expected = expected )

    def poi_upper_limit ( self, expected : Union [ bool, Text ],
           limit_on_xsec : bool = False ) -> float:
        """ simple frontend, to spey::poi_upper_limit 
        :param limit_on_xsec: if True, then return the limit on the
                              cross section
        """

        init, bounds, args = self.getSpeyInitialisation ( False, False )
        expected = self.translateExpectationType ( expected )
        try:
            ret = self.statModel.poi_upper_limit ( expected, par_bounds = bounds,
                init_pars = init, **args )
        except ValueError as e:
            logger.warning ( f"when computing upper limit for SL: {e}. Will try with other method" )
            self.alternateMethod ( args )
            ret = statModel.poi_upper_limit ( expected=expected, par_bounds=bounds, 
                init_pars = init, **args )
        if limit_on_xsec and type(ret) not in [ None ]:
            ret = ret * self.xsection * fb
        return ret

    def asimov_likelihood ( self, poi_test : float, expected : Union[bool,Text], 
                            return_nll : bool ) -> float:
        """ simple frontend to spey functionality """
        self.checkMinimumPoi( poi_test )
        expected = self.translateExpectationType ( expected )
        init, bounds, args = self.getSpeyInitialisation ( True )
        return self.statModel.asimov_likelihood ( poi_test = poi_test,
            expected = expected, return_nll = return_nll, par_bounds = bounds,
            init_pars = init, **args )

    def checkMinimumPoi ( self, poi_test : float ):
        """ check if poi is below minimum_poi """
        config = self.statModel.backend.model.config()
        if poi_test < config.minimum_poi:
            logger.error ( f'Calling likelihood for {dataset.globalInfo.id} (using combination of SRs) for a mu giving a negative total yield. mu = {mu} and minimum_mu = {config.minimum_poi}.' )

    def likelihood ( self, poi_test : float, expected : Union[bool,Text], 
                            return_nll : bool ) -> float:
        """ simple frontend to spey functionality """
        self.checkMinimumPoi ( poi_test )
        expected = self.translateExpectationType ( expected )
        init, bounds, args = self.getSpeyInitialisation ( True )
        return self.statModel.likelihood ( poi_test = poi_test,
            expected = expected, return_nll = return_nll, par_bounds = bounds,
            init_pars = init, **args )

    def maximize_likelihood ( self, expected : Union[bool,Text],
           allowNegativeSignals : bool = True,
           return_nll : bool = False  ) -> Tuple[float,float]:
        """ simple frontend to spey functionality
        :param return_nll: if True, return negative log likelihood
        :param allowNegativeSignals: allow also negative muhats
        :returns: tuple of muhat,lmax
        """
        expected = self.translateExpectationType ( expected )
        init, bounds, args = self.getSpeyInitialisation ( allowNegativeSignals = allowNegativeSignals )
        ret = self.statModel.maximize_likelihood ( expected = expected, 
            par_bounds = bounds, return_nll = return_nll, init_pars = init, **args )
        return ret

    def likelihood_timothee ( self, poi_test : float, expected : Union[bool,Text], 
                            return_nll : bool ) -> float:
        """ simple frontend to spey functionality
        likelihood, timothee style. i.e. if expected is
        posteriori then compute asimov likelihood but with expected="apriori"
        will worry later about the correctness of this. """
        if expected == "posteriori":
            return self.asimov_likelihood ( poi_test = poi_test, expected="apriori",
                return_nll = return_nll )
        return self.likelihood ( poi_test = poi_test, expected = expected, return_nll = return_nll )

    def sigma_mu ( self, poi_test : float, expected : Union[bool,Text], allowNegativeSignals : bool = False ):
        """ determine sigma at poi_test.
        :param: FIXME allowNegativeSignals should not be needed!
        """
        test_statistic = "q" if allowNegativeSignals else "qmutilde"
        sigma_mu = self.statModel.sigma_mu( poi_test=poi_test,expected=expected,
                                            test_statistics=test_statistic )
        return sigma_mu

    def maximize_asimov_likelihood ( self, expected : Union[bool,Text],
           return_nll : bool = False ) -> Tuple[float,float]:
        """ simple frontend to spey functionality 
        :param return_nll: if True, return negative log likelihood
        :param allowNegativeSignals: allow also negative muhats
        :returns: tuple of muhat,lmax
        """
        expected = self.translateExpectationType ( expected )
        init, bounds, args = self.getSpeyInitialisation ( True )
        return self.statModel.maximize_asimov_likelihood ( expected = expected, 
            par_bounds = bounds, init_pars = init, return_nll = return_nll,
            test_statistics="qmutilde", **args )

    def getSpeyInitialisation ( self, allowNegativeSignals : bool = False,
                                initial_bracket : bool = True ):
        """ get decent initial bounds and initial values for a statModel
        :param allowNegativeSignals: if true, then bound the poi to positive values
        :param initial_bracket: also supply an initial bracket for CLs-alpha root
                                finding?
        """
        if self.dataset.getType()=="efficiencyMap":
            ini = self.getInitialisationForSingleRegions ( allowNegativeSignals )
            return self.filterInitialBracket ( ini, initial_bracket )
            
        if self.dataset.type=="simplified":
            ini = self.getInitialisationForSL ( allowNegativeSignals )
            return self.filterInitialBracket ( ini, initial_bracket )
        if dataset.type!="pyhf":
            raise Exception ( f"dont know that dataset type {dataset.type}" )
        ini = self.getInitialisationForPyhf ( dataset, allowNegativeSignals )
        return self.filterInitialBracket ( ini, initial_bracket )

    def filterInitialBracket ( self, args : tuple, initial_bracket : bool ):
        """ possibly filter out the suggestion for an initial bracket,
            i.e. remove low_init and hig_init from args[2] "
        :param initial_bracket: if True, then actually leave them in, do nothing
        """
        if initial_bracket:
            return args
        args[2].pop ( "low_init", None )
        args[2].pop ( "hig_init", None )
        return args

    def getInitialisationForSingleRegions ( self, allowNegativeSignals : bool = False ) -> Tuple[List,List,Dict]:
        """ here we globally steer the initialisation for the case
            of single region models """
        model = self.statModel.backend.model
        config = model.config()
        init = config.suggested_init
        bounds = config.suggested_bounds
        args = {}
        return init,bounds,args

    def getInitialisationForPyhf ( self, allowNegativeSignals : bool = False ):
        """ get decent initial bounds and initial values for a statModel
        :param allowNegativeSignals: if true, then bound the poi to positive values
        """
        statModel = dataset.statModel
        model = statModel.backend.model
        config = model.config()
        init = config.suggested_init
        bounds = config.suggested_bounds
        if False: ## print it
            print ( "datasets", [ x.dataInfo.dataId for x in dataset._datasets ] )
            bg = [ x.dataInfo.expectedBG for x in dataset._datasets ]
            signal = model.signal[0]["value"]["data"]
            print ( "signals", model.signal[0]["value"]["data"] )
            print ( "nobs", model.background["channels"][0]["samples"][0]["data"] )
            print ( "suggested bounds", bounds )
            if allowNegativeSignals:
                bounds[config.poi_index] = (config.minimum_poi, 100)
            else:
                bounds[config.poi_index] = (0, 100)
            mui = []
            for i in range(len(nobs)):
                mui.append ( (nobs[i]-bg[i])/signal[i] )
            print ( "mui", mui )
        # args = { "maxiter": 500, "method": "SLSQP" } ## extra args for the optimizers
        # method BFGS, SLSQP
        #args = { "maxiter": 500, "method": "BFGS", "ntrials": 3,
        #         "xrtol": 1e-6 "low_init": bounds[0][0], 
    #                "hig_init": bounds[0][1] 
        #}
        # args = { "maxiter": 500, "ntrials": 1, "method": "SLSQP" }
        # args = { "maxiter": 500, "ntrials": 1, "method": None }
        args = {}
        return init,bounds,args

    def alternateMethod ( self, args : dict ):
        """ try a different method. i hope we wont need this in the long run """
        if not "method" in args or args["method"] is None:
            args["method"]="BFGS"
            return
        if args["method"]=="SLSQP":
            args["method"]="BFGS"
            return
        if args["method"]=="BFGS":
            args["method"]="SLSQP"
            return
        args["method"]="L-BFGS-B"

    @property
    def xsection(self):
        return self.statModel.xsection

    def getInitialisationForSL ( self, allowNegativeSignals : bool = False ):
        """ get decent initial bounds and initial values for an SL statModel
        :param allowNegativeSignals: if true, then bound the poi to positive values
        """
        import numpy as np
        statModel = self.statModel
        config = statModel.backend.model.config()
        init = config.suggested_init
        bounds = config.suggested_bounds
        args = {} ## extra args for the optimizers
        if False: ## these were the old values!
            bounds = [(suggested[0]-800,suggested[1]+800) for suggested in config.suggested_bounds]
            bounds[config.poi_index] = (-520, 800) # for now!
            return bounds,init,args
        assert config.poi_index == 0, f"Error: I assume the poi index to be zero, not {config.poi_index}"
        init[1:] = statModel.backend.model.observed - statModel.backend.model.background - statModel.backend.model.signal
        numerator, denominator = [], []
        for i in range ( len ( statModel.backend.model.observed ) ):
            cov = statModel.backend.model.covariance[i][i]
            x = np.sqrt ( cov )
            bounds[i+1]=(-5*x+init[i+1],5*x+init[i+1])
            sig = statModel.backend.model.signal[i]
            obs = statModel.backend.model.observed[i]
            bg = statModel.backend.model.background[i]
            if sig > 0. and cov > 0.:
                # for the given region, mui would be the best bet for muhat
                mui = ( obs - bg ) / sig
                # its variance is this
                cov_mui =  ( obs + cov ) / ( sig**2 )
                # the inverse of which will be our weights
                wi = 1. / cov_mui
                numerator.append ( wi * mui )
                denominator.append ( wi )
        ## ok so the bounds should be -5*x,5*x with x being np.sqrt(statModel.backend.model.covariance[i][i], the initial values just the diff between observation and expectation
        init_muhat = np.sum ( numerator ) / np.sum ( denominator)
        if not allowNegativeSignals and init_muhat<0.:
            init_muhat = 0.
        err_muhat = np.sqrt ( len(denominator) / np.sum ( denominator ) )
        # err_muhat = np.sqrt ( np.sum ( totweight )**(-2) * np.sum ( covest ) )

        init[config.poi_index] = init_muhat
        minmu, maxmu = -5*err_muhat + init_muhat, 5*err_muhat + init_muhat
        if not allowNegativeSignals and minmu < 0.:
            minmu = 0.
        bounds [ config.poi_index ] = ( minmu, maxmu )
        if False:
            print ( f"residual ({len(init)}) is", init_muhat, "+-", err_muhat )
            print ( "errors are", covest )
            print ( "residuals are", residuals )
            print ( "init pars in srCombinations are", init[:3] )
            print ( "signals   in srCombinations are", statModel.backend.model.signal[:] )
            print ( "deltas  in srCombinations are", statModel.backend.model.observed - statModel.backend.model.background )
            print ( "bounds    in srCombinations are", bounds[:3] )
        args = { "maxiter": 500, "method": None, "ntrials": 1,
                 "low_init": bounds[0][0], 
                    "hig_init": bounds[0][1] }
        # args["method"]="BFGS"
        args["tol"]=1e-3
        # args["xrtol"]=1e-6
        # print ( f"speyTools: initbracket is", args["low_init"], args["hig_init"] )
        args = { "maxiter": 500, "ntrials": 1, "method": None }
        return init,bounds,args

