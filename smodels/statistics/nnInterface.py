#!/usr/bin/env python3

"""
.. module:: nnInterface
   :synopsis: Code that delegates the computation of limits and likelihoods to machine learned models

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from typing import Union, Text, Tuple, Callable, Dict
import copy
import numpy as np
import sys
import onnxruntime
from smodels.base.smodelsLogging import logger
from smodels.statistics.basicStats import determineBrentBracket, CLsfromNLL, CLsfromLsb, getNumericalVariance
from scipy import optimize, differentiate

nninfo = {
    "hasgreeted": False,
}

class NNData:
    """
    Holds data for use in the machine learned models
    :ivar nsignals: signal predictions list divided into sublists, one for each json file
    """

    def __init__(self, nsignals, dataObject ):
        self.nsignals = nsignals  # fb
        self.getTotalYield()
        self.dataObject = dataObject
        self.globalInfo = dataObject.globalInfo
        datasets = [ds.getID() for ds in self.dataObject.origdatasets]
        self.origDataSetOrder = datasets

    def getTotalYield ( self ):
        """ the total yield in all signal regions """
        S = 0
        for signal in self.nsignals.values():
            if isinstance(signal, list):
                for sig in signal:
                    S += sig
            else:
                S += signal
        self.totalYield = S

class NNUpperLimitComputer:
    """
    Class that computes the upper limit using the jsons files and signal
    informations in the 'data' instance of 'NNData'
    """

    def __init__(self, data, cl=0.95, lumi=None ):
        """

        :param data: instance of 'NNData' holding the signals information
        :param cl: confdence level at which the upper limit is desired to be computed
        :ivar data: created from data
        :ivar nsignals: signal predictions list divided into sublists, one for each json file
        :ivar zeroSignalsFlag: list boolean flags in case all signals are zero for a specific json
        :ivar cl: created from cl
        :ivar alreadyBeenThere: boolean flag that identifies when nsignals accidentally passes twice at two identical values
        """

        self.data = data
        import onnxruntime
        self.regressors = {}
        for jsonfilename,onnx in self.data.globalInfo.onnxes.items():
            sess = onnxruntime.InferenceSession (  onnx )
            self.regressors[jsonfilename]={"session": sess,
                "dim": sess.get_inputs()[0].shape[1] }
        #self.regressor = onnxruntime.InferenceSession ( self.data.globalInfo.onnx )
        # store the dimensionality of the input vector that the model
        # asks us to. this may be different from the number of our
        # signal regions, as control regions may have been added.
        # we will pad with zeroes
        # self.regressor_dim = self.regressor.get_inputs()[0].shape[1]
        self.lumi = lumi
        self.nsignals = copy.deepcopy ( self.data.nsignals )
        logger.debug("Signals : {}".format(self.nsignals))
        # self.zeroSignalsFlag = self.data.zeroSignalsFlag
        self.cl = cl
        self.sigma_mu = None
        self.alreadyBeenThere = (
            False  # boolean to detect wether self.signals has returned to an older value
        )
        self.welcome()

    def getMostSensitiveModel ( self ):
        if hasattr ( self, "mostSensitiveModel" ):
            return self.mostSensitiveModel
        jsonfiles = list(self.data.globalInfo.onnxMeta.keys())
        if len(jsonfiles)==1:
            self.mostSensitiveModel = jsonfiles[0]
            return self.mostSensitiveModel
        mumin,modelToUse=float("inf"),None
        for model in jsonfiles:
            mu = model.getUpperLimitOnMu ( expected=True, modelToUse = model )
            if mu < mumin:
                modelToUse = model
                mumin = mu
        self.mostSensitiveModel = modelToUse
        self.mumin = mumin
        # print ( f"@@0 for now use first model {modelToUse}" )
        return self.mostSensitiveModel

    def isControlRegion ( self, srname : str, modelToUse : str ) -> bool:
        """ check if srname is control region
        :returns: true if srname is control region
        """
        def isCR ( dct : Dict )-> bool:
            """ tiny helper """
            if not "type" in dct:
                return False
            return dct["type"]=="CR"

        for dct in self.data.globalInfo.mlModels[modelToUse]:
            # print ( f"@@ {dct}" )
            name = dct["onnx"]
            if name == srname:
                return isCR ( dct )
            pname = name+"-0"
            if pname == srname:
                return isCR ( dct )
        return False

    def negative_log_likelihood(self, poi_test,
        modelToUse : Union[None,str] = None,
        outputType : str = "extended" ):
        """ the method that really wraps around the llhd computation.
        :param modelToUse: if given, compute the nll for that model.
        If None compute for most sensitive analysis.
        :param outputType: if 'extended' return dictionary with all
        values, if 'observed' return nll_obs_1, if 'expected' return
        nll_exp_1

        :returns: dictionary with nlls, obs and exp, mu=0 and 1
        """
        try:
            poi_test = poi_test[0]
        except (TypeError,IndexError) as e:
            pass

        if modelToUse == None:
            modelToUse = self.getMostSensitiveModel ( )

        syields = []
        for srname,smyield in self.data.globalInfo.onnxMeta[modelToUse]["smYields"].items():
            p1 = srname.rfind("-")
            realname = srname[:p1]
            if not realname in self.nsignals:
                realname = f"{realname}[{srname[p1+1:]}]"
                assert realname in self.nsignals, \
                  f"cannot find sr name {realname} in {self.nsignals}"
            # smodelsname = self.data.globalInfo
            signal = float ( self.nsignals[realname]*poi_test )
            if self.isControlRegion ( srname, modelToUse ):
                if hasattr ( self.data.globalInfo, "includeCRs" ) and self.data.globalInfo.includeCRs == False:
                    continue
                obsyield = self.data.globalInfo.onnxMeta[modelToUse]["obsYields"][srname]
                ## seems like a CR! replaced bkgexpected with observed (postfit)
                smyield = obsyield
            tot = smyield + signal
            syields.append ( tot )
            # print ( f"@@0 the smyield of {srname} is {smyield} poi_test is {poi_test} signal {signal}" )

        #for i in range(4):
        #    syields.append(0.)

        scaled_signal_yields = np.array( [syields], dtype=np.float32 )

        #if poi_test == 0.:
        #    print ( f"inputMeans {self.data.globalInfo.inputMeans[modelToUse]}" )
        #    print ( f"inputErrors {self.data.globalInfo.inputErrors[modelToUse]}" )
        for i,x in enumerate(scaled_signal_yields[0]):
            t = 0. # x
            err = self.data.globalInfo.onnxMeta[modelToUse]["inputErrors"][i]
            if err > 1e-20:
                t = (x - self.data.globalInfo.onnxMeta[modelToUse]["inputMeans"][i])/err
            #else:
            #    t = # - self.data.globalInfo.inputMeans[i]
            scaled_signal_yields[0][i]=t

        #if poi_test == 0.:
        #    print ( f"@@X we evaluate at {scaled_signal_yields}" )
        if len(scaled_signal_yields[0])!=self.regressors[modelToUse]["dim"]:
            dim_nn = self.regressors[modelToUse]["dim"]
            dim_input = len(scaled_signal_yields[0])
            line=f"the network wants {dim_nn} input dimensions, but we supply {dim_input}. fix it!"
            logger.error ( "[nnInterface]", line )
            print ( line )
            sys.exit()
        arr = self.regressors[modelToUse]["session"].run(None,
                {"input_1":scaled_signal_yields})
        # print ( f"@@arr {arr}" )
        arr = arr[0][0]
        nll0obs =  self.data.globalInfo.onnxMeta[modelToUse]["nLL_obs_mu0"]
        nll0exp =  self.data.globalInfo.onnxMeta[modelToUse]["nLL_exp_mu0"]
        nllA0obs =  self.data.globalInfo.onnxMeta[modelToUse]["nLLA_obs_mu0"]
        nllA0exp =  self.data.globalInfo.onnxMeta[modelToUse]["nLLA_exp_mu0"]
        i_exp, i_obs, i_expA, i_obsA = -4, -3, -2, -1 # the indices
        expDelta = self.data.globalInfo.onnxMeta[modelToUse]["inputMeans"][i_exp]
        obsDelta = self.data.globalInfo.onnxMeta[modelToUse]["inputMeans"][i_obs]
        expDeltaA = self.data.globalInfo.onnxMeta[modelToUse]["inputMeans"][i_expA]
        obsDeltaA = self.data.globalInfo.onnxMeta[modelToUse]["inputMeans"][i_obsA]
        expErr = self.data.globalInfo.onnxMeta[modelToUse]["inputErrors"][i_exp]
        obsErr = self.data.globalInfo.onnxMeta[modelToUse]["inputErrors"][i_obs]
        expErrA = self.data.globalInfo.onnxMeta[modelToUse]["inputErrors"][i_expA]
        obsErrA = self.data.globalInfo.onnxMeta[modelToUse]["inputErrors"][i_obsA]
        nll1exp = nll0exp + arr[i_exp]*expErr + expDelta
        nll1obs = nll0obs + arr[i_obs]*obsErr + obsDelta
        #print ( f"@@5 onnxMeta", self.data.globalInfo.onnxMeta )
        nllA1exp = nllA0exp + arr[i_expA]*expErrA + expDeltaA
        nllA1obs = nllA0obs + arr[i_obsA]*obsErrA + obsDeltaA

        if False and poi_test == 0.:
            print ( f"@@5 nll0obs {nll0obs} {nll0exp}" )
            print ( f"@@5 nllA0obs {nllA0obs} {nllA0exp}" )
            print ( f"@@5 poi_test {poi_test}" )
            print ( f"@@5 nll1obs {float(nll1obs)} {float(nll1exp)}" )
            print ( f"@@5 nllA1obs {float(nllA1obs)} {float(nllA1exp)}" )

        ret = { "nll_exp_0": nll0exp, "nll_exp_1": float(nll1exp),
                "nll_obs_0": nll0obs, "nll_obs_1": float(nll1obs),
                "nllA_exp_0": nllA0exp, "nllA_exp_1": float(nllA1exp),
                "nllA_obs_0": nllA0obs, "nllA_obs_1": float(nllA1obs) }
        if abs(poi_test)<1e-10:
            if abs(nll0obs-nll1obs)>1e-1:
                logger.error ( f"mu={poi_test:.2f} but nll0obs {nll0obs:.4f}!= nll1obs {nll1obs:.4f}. obsDelta {obsDelta} obsErr {obsErr} arr {arr}" )
                # ret["nll_obs_1"]=nll0obs
            if abs(nll0exp-nll1exp)>1e-1:
                logger.error ( f"mu={poi_test:.2f} but nll0exp {nll0exp:.4f}!= nll1exp {nll1exp:.4f}." )
                # ret["nll_exp_1"]=nll0exp
            if False:
                ret["nll_exp_1"]=ret["nll_exp_0"]
                ret["nll_obs_1"]=ret["nll_obs_0"]
                ret["nllA_obs_1"]=ret["nllA_obs_0"]
                ret["nllA_exp_1"]=ret["nllA_exp_0"]
        if outputType == "observed":
            return ret["nll_obs_1"]
        if outputType == "expected":
            return ret["nll_exp_1"]
        if outputType != "extended":
            logger.error ( f"outputType {outputType} unknown. should be one of 'observed', 'expected', 'extended'." )
            sys.exit(-1)
        return ret

    def welcome(self):
        """
        greet the world
        """

        if nninfo["hasgreeted"]:
            return
        logger.info( f"NN interface, we are using onnxruntime v{onnxruntime.__version__}" )
        nninfo["hasgreeted"] = True

    def likelihood( self, mu=1.0, return_nll=False, expected=False,
              modelToUse : Union[None,str] = None, asimov : bool = False ):
        """
        Returns the value of the likelihood. \
        Inspired by the 'pyhf.infer.mle' module but for non-log likelihood

        :param return_nll: if true, return nll, not llhd
        :param expected: if False, compute expected values, if True, \
            compute a priori expected, if "posteriori" compute posteriori \
            expected
        :param modelToUse: if given, compute likelihood for that model.
        :param asimov: if true, compute for asimov data
        If None compute for most sensitive analysis.
        """
        ret = self.negative_log_likelihood(mu,modelToUse=modelToUse)
        if expected:
            nll = ret['nll_exp_1']
            if asimov:
                nll = ret['nllA_exp_1']
        else:
            nll = ret['nll_obs_1']
            if asimov:
                nll = ret['nllA_obs_1']

        logger.debug( f"Calling likelihood")
        return self.exponentiateNLL ( nll, not return_nll )

    def exponentiateNLL(self, nll, doIt = True ):
        """if doIt, then compute likelihood from nll,
        else return nll"""
        if nll == None:
            return None
            #if doIt:
            #    return 0.0
            #return 9000.0
        if doIt:
            return np.exp(-nll )
        return nll

    def lmax( self, return_nll=False, expected=False,
              allowNegativeSignals=True,
              modelToUse : Union[None,str] = None,
              asimov : bool = False ):
        """
        Returns the (negative log) max likelihood

        :param return_nll: if true, return nll, not llhd
        :param expected: if False, compute expected values, if True,
        compute a priori expected, if "posteriori" compute posteriori
        expected
        :param allowNegativeSignals: if False, then negative nsigs are
        replaced with 0.
        :param modelToUse: if given, compute lmax for that model.
        :param asimov: if true, compute for asimov data
        If None compute for most sensitive analysis.
        """
        if modelToUse == None:
            modelToUse = self.getMostSensitiveModel ( )
        # FIXME maximize this one function
        muhat,nllmin = self.data.globalInfo.onnxMeta[modelToUse]["nLL_obs_max"]
        if asimov:
            muhat,nllmin = self.data.globalInfo.onnxMeta[modelToUse]["nLLA_obs_max"]
            if expected:
                muhat,nllmin = self.data.globalInfo.onnxMeta[modelToUse]["nLLA_exp_max"]
        elif expected:
            muhat,nllmin = self.data.globalInfo.onnxMeta[modelToUse]["nLL_exp_max"]

        outputType = "observed"
        if expected:
            outputType = "expected"
        #print ( f"@@ outputType {outputType} " )
        #bounds=[(-1,10)]
        #def callme ( args ):
        #    print ( f"@@ callme {args}" )
        options = { "disp": False, "maxiter": 200 }

        ## FIXME compute sigma_mu, compute via nllA

        def myNLL ( x ):
            # print ( f"@@0 myNLL x={x}" )
            if type(x) in [ list, np.array, np.ndarray ]:
                ret = []
                for xi in x:
                    ret.append ( myNLL ( xi ) )
                return np.array ( ret )
            ret = self.negative_log_likelihood ( x, modelToUse=modelToUse,
                                                 outputType=outputType )
            # print ( f"@@0 myNLL ret={x}" )
            return ret

        method = "Nelder-Mead"
        initx0s = [ 0., .1, -.1, .3, -.3, 1., -1., 3., -3., 10., -10., 100,-100 ]
        bounds=[(-100,100)]
        if not allowNegativeSignals:
            bounds=[(0,100)]
        for x0 in initx0s:
            o = optimize.minimize ( self.negative_log_likelihood, x0=x0,
                    args=(modelToUse,outputType), tol=1e-8, options = options,
                    method = method, bounds=bounds )
            #if o.fun < 0:
            #    print ( f"@@ o {o}" )
            if o.success == True and o.fun>0:
                muhat, nllmin = o.x[0], o.fun
                # import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()
                o = differentiate.hessian ( myNLL, np.array ( [ muhat ] ) )
                hessian = o.ddf[0][0][0]
                sigma_mu = 0.
                if hessian > 0.:
                    sigma_mu = np.sqrt ( 1. / hessian )

                ret = { "nll_min": nllmin, "muhat": muhat, "sigma_mu": sigma_mu }
                # print ( f"@@ muhat {muhat} nllmin {nllmin} hessian {hessian} allowNegativeSignals {allowNegativeSignals}" )
                # print ( f"@@ negative_log_likelihood {self.negative_log_likelihood(muhat)}" )
                return ret
            if x0 == initx0s:
                method = "L-BFGS-B"
        #print ( f"could not find nll_min" )
        #print ( f"nll(0) {self.negative_log_likelihood(0.)}" )
        #print ( f"nll(.1) {self.negative_log_likelihood(.1)}" )
        #print ( f"nll(-.1) {self.negative_log_likelihood(-.1)}" )
        #sys.exit()
        logger.warning ( f"could not find nll_min!" )
        return None

    def getUpperLimitOnSigmaTimesEff(self, expected=False,
            modelToUse : Union[None,str] = None ):
        """
        Compute the upper limit on the fiducial cross section sigma times efficiency:

        :param expected:  - if set to 'True': uses expected SM backgrounds as signals
                          - else: uses 'self.nsignals'
        :param modelToUse: if given, compute the nll for that model.
        If None compute for most sensitive analysis.
        :return: the upper limit on sigma times eff at 'self.cl' level (0.95 by default)
        """
        if self.data.totalYield == 0.:
            return None
        else:
            ul = self.getUpperLimitOnMu( expected=expected, modelToUse=modelToUse )
            if ul == None:
                return ul
            if self.lumi is None:
                logger.error(f"asked for upper limit on fiducial xsec, but no lumi given with the data")
                return ul
            xsec = self.data.totalYield / self.lumi
            return ul * xsec

    def getUpperLimitOnMu(self, expected=False, allowNegativeSignals = False,
            modelToUse : Union[None,str] = None ) -> float:
        """
        Compute the upper limit on the signal strength modifier with:
            - by default, the combination of the workspaces contained into self.workspaces

        :param expected:  - if set to 'True': uses expected SM backgrounds as signals
                          - else: uses 'self.nsignals'
        :param modelToUse: if given, compute the nll for that model.
        If None compute for most sensitive analysis.
        :return: the upper limit at 'self.cl' level (0.95 by default)
        """
        if modelToUse == None:
            if hasattr ( self, "mumin" ):
                return self.mumin
            modelToUse = self.getMostSensitiveModel()
        mu_hat, sigma_mu, clsRoot = self.getCLsRootFunc(expected=expected,
                              allowNegativeSignals=allowNegativeSignals,
                              modelToUse = modelToUse)
        a, b = determineBrentBracket(mu_hat, sigma_mu, clsRoot,
                allowNegative = allowNegativeSignals )
        mu_lim = optimize.brentq(clsRoot, a, b, rtol=1e-03, xtol=1e-06)
        if expected == "posteriori":
            print ( f"@@nnInterface.getUpperLimitOnMu r={1./mu_lim:.3f} expected {expected}" )
        return mu_lim

    def getCLsRootFunc(self, expected: bool = False,
            allowNegativeSignals : bool = True,
            modelToUse : Union[None,str] = None ) -> Tuple[float, float, Callable]:
        """
        Obtain the function "CLs-alpha[0.05]" whose root defines the upper limit,
        plus mu_hat and sigma_mu

        :param expected: if True, compute expected likelihood, else observed
        :param modelToUse: if given, compute the nll for that model.
        If None compute for most sensitive analysis.
        """
        # a posteriori expected is needed here
        # mu_hat is mu_hat for signal_rel
        fmh = self.lmax(expected="posteriori", allowNegativeSignals=allowNegativeSignals,
                             return_nll=True, modelToUse = modelToUse )
        mu_hat, sigma_mu, nll0A = fmh["muhat"], fmh["sigma_mu"], fmh["nll_min"]

        nll0 = nll0A
        
        if expected != "posteriori":
            fmh = self.lmax(expected=expected, allowNegativeSignals=allowNegativeSignals, modelToUse = modelToUse )
            mu_hat, sigma_mu, _ = fmh["muhat"], fmh["sigma_mu"], fmh["nll_min"]
            mu_hat = mu_hat if mu_hat is not None else 0.0
            nll0 = self.likelihood(mu_hat, expected=expected, return_nll=True,
                    modelToUse = modelToUse )

        # logger.error ( f"COMB nll0A {nll0A:.3f} mu_hatA {mu_hatA:.3f}" )
        # return 1.
        # print ( f"@@X getCLsRootFunc expected {expected} nll0 {nll0:.1f} nll0A {nll0A:.1f}" )

        def clsRootTevatron( mu: float, return_type: Text = "CLs-alpha",
                     modelToUse : Union[None,str] = None ) -> float:
            # at - infinity this should be .95,
            # at + infinity it should -.05
            # Make sure to always compute the correct llhd value (from theoryPrediction)
            # and not used the cached value (which is constant for mu~=1 an mu~=0)
            nll_b = self.likelihood(0., return_nll=True, expected=expected, modelToUse = modelToUse )
            if nll_b is None:
                return None
            sigma2_b = getNumericalVariance ( self, 0. )
            nll_sb = self.likelihood(mu, return_nll=True, expected=expected, modelToUse = modelToUse )
            if nll_sb is None:
                return None
            sigma2_sb = getNumericalVariance ( self, 1. )
            CLs = CLsfromLsb(nll_sb, nll_b, sigma2_sb, sigma2_b, return_type=return_type)
            return CLs

        def clsRootAsimov( mu: float, return_type: Text = "CLs-alpha",
                     modelToUse : Union[None,str] = None ) -> float:
            # at - infinity this should be .95,
            # at + infinity it should -.05
            # Make sure to always compute the correct llhd value (from theoryPrediction)
            # and not used the cached value (which is constant for mu~=1 an mu~=0)
            nllA = self.likelihood(mu, return_nll=True, modelToUse = modelToUse, asimov = True )
            nll = nllA
            if expected != "posteriori":
                nll = self.likelihood(mu, return_nll=True, expected=expected, modelToUse = modelToUse )
            return CLsfromNLL(nllA, nll0A, nll, nll0, (mu_hat > mu), return_type=return_type) if nll is not None else None

        from smodels.base import runtime
        useTevatron = runtime.experimentalFeature ( "tevatroncls" )
        if useTevatron:
            return mu_hat, sigma_mu, clsRootTevatron
        return mu_hat, sigma_mu, clsRootAsimov
