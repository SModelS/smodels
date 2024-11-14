#!/usr/bin/env python3

"""
.. module:: nnInterface
   :synopsis: Code that delegates the computation of limits and likelihoods to machine learned models

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from typing import Union, Text, Tuple, Callable
import copy
import numpy as np
import sys
import onnxruntime
from smodels.base.smodelsLogging import logger
from smodels.statistics.basicStats import determineBrentBracket, CLsfromNLL, CLsfromLsb
from scipy import optimize

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
        jsonfiles = list(self.data.globalInfo.smYields.keys())
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

    def negative_log_likelihood(self, poi_test,
        modelToUse : Union[None,str] = None):
        """ the method that really wraps around the llhd computation.
        :param modelToUse: if given, compute the nll for that model.
        If None compute for most sensitive analysis.

        :returns: dictionary with nlls, obs and exp, mu=0 and 1
        """

        #for modelname,model in self.data.globalInfo.smYields.items():
        #    compute_upper_limit(model)
        #choose_most_sensitive_model
        if modelToUse == None:
            modelToUse = self.getMostSensitiveModel ( )

        syields = []
        for srname,smyield in self.data.globalInfo.smYields[modelToUse].items():
            p1 = srname.rfind("-")
            realname = srname[:p1]
            # ic ( realname )
            if not realname in self.nsignals:
                realname = f"{realname}[{srname[p1+1:]}]"
                assert realname in self.nsignals, \
                  f"cannot find sr name {realname} in {self.nsignals}"
            smodelsname = self.data.globalInfo
            signal = float ( self.nsignals[realname]*poi_test )
            tot = smyield + signal
            syields.append ( tot )

        scaled_signal_yields = np.array( [syields], dtype=np.float32 )

        for i,x in enumerate(scaled_signal_yields[0]):
            t = 0. # x
            err = self.data.globalInfo.inputErrors[modelToUse][i]
            if err > 1e-20:
                t = (x - self.data.globalInfo.inputMeans[modelToUse][i])/err
            #else:
            #    t = # - self.data.globalInfo.inputMeans[i]
            scaled_signal_yields[0][i]=t

        arr = self.regressors[modelToUse]["session"].run(None, {"input_1":scaled_signal_yields})
        arr = arr[0][0]
        nll0obs =  self.data.globalInfo.nll_obs_mu0[modelToUse]
        nll0exp =  self.data.globalInfo.nll_exp_mu0[modelToUse]
        expDelta = self.data.globalInfo.inputMeans[modelToUse][-2]
        obsDelta = self.data.globalInfo.inputMeans[modelToUse][-1]
        expErr = self.data.globalInfo.inputErrors[modelToUse][-2]
        obsErr = self.data.globalInfo.inputErrors[modelToUse][-1]
        nll1exp = nll0exp + arr[0]*expErr + expDelta
        nll1obs = nll0obs + arr[1]*obsErr + obsDelta
        ret = { "nll_exp_0": nll0exp, "nll_exp_1": nll1exp,
                "nll_obs_0": nll0obs, "nll_obs_1": nll1obs }
        return ret

    def welcome(self):
        """
        greet the world
        """

        if nninfo["hasgreeted"]:
            return
        logger.info( f"NN interface, we are using xxx" )
        nninfo["hasgreeted"] = True

    def likelihood( self, mu=1.0, return_nll=False, expected=False,
              modelToUse : Union[None,str] = None ):
        """
        Returns the value of the likelihood. \
        Inspired by the 'pyhf.infer.mle' module but for non-log likelihood

        :param return_nll: if true, return nll, not llhd
        :param expected: if False, compute expected values, if True, \
            compute a priori expected, if "posteriori" compute posteriori \
            expected
        :param modelToUse: if given, compute likelihood for that model.
        If None compute for most sensitive analysis.
        """
        ret = self.negative_log_likelihood(mu,modelToUse=modelToUse)
        if mu==0.0 and expected:
            nll = ret['nll_exp_0']
        elif mu==0.0 and not expected:
            nll = ret['nll_obs_0']
        elif mu!=0.0 and expected:
            nll = ret['nll_exp_1']
        elif mu!=0.0 and not expected:
            nll = ret['nll_obs_1']
        else:
            raise NotImplementedError(f'Request for {mu} received but only m=0 and mu=1 are implemented.')

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
              modelToUse : Union[None,str] = None ):
        """
        Returns the (negative log) max likelihood

        :param return_nll: if true, return nll, not llhd
        :param expected: if False, compute expected values, if True,
        compute a priori expected, if "posteriori" compute posteriori
        expected
        :param allowNegativeSignals: if False, then negative nsigs are
        replaced with 0.
        :param modelToUse: if given, compute lmax for that model.
        If None compute for most sensitive analysis.
        """
        if modelToUse == None:
            modelToUse = self.getMostSensitiveModel ( )
        nll0 = self.likelihood ( 0., expected = expected, return_nll = True )
        def getNLL ( mu ):
            d = self.negative_log_likelihood (mu)
            if mu==0.0 and expected:
                nll = d['nll_exp_0']
            elif mu==0.0 and not expected:
                nll = d['nll_obs_0']
            elif mu!=0.0 and expected:
                nll = d['nll_exp_1']
            elif mu!=0.0 and not expected:
                nll = d['nll_obs_1']

            return nll

        ret = { "nll_min": nll0, "muhat": 0., "sigma_mu": 0. }

        muini = [ 0.0, 1.0, 3.0, -1.0, 10.0, 0.1, -0.1 ]
        if not allowNegativeSignals:
            muini = [ 0.0, 1.0, 3.0, 10.0, 0.1 ]
        for mu0 in muini:
            o = optimize.minimize ( getNLL, mu0, method="Nelder-Mead" )
            # print ( f"@@6 o={o} nll0 {nll0}" )
            # sys.exit()
            if o.success == True:
                muhat = o.x[0]
                if not allowNegativeSignals and muhat < 0.:
                    continue
                nll_min = o.fun
                #    muhat = 0.
                #    nll_min = nll0
                if 0. < nll_min < ret["nll_min"]:
                    ret = { "nll_min": nll_min, "muhat": muhat }
        ret["nll0"]=nll0
        ret["lmax"]=ret["nll_min"] # FIXME the lmax name!!!
        ret["expected"]=expected
        o = optimize.minimize ( getNLL, ret["muhat"], method="BFGS" )
        sigma_mu = np.sqrt ( o.hess_inv[0][0] )
        ret["sigma_mu"]=sigma_mu

        return ret

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
        # print ( f"@@5 getUpperLimitOnMu mu_hat {mu_hat} {sigma_mu} mu_lim {mu_lim}" )
        return mu_lim


    def getVarForMu ( self, mu : float ) -> float:
        """ get the variance around mu """
        # np.diff(np.diff([x*x for x in range(0,10)]))
        # print ( f"@@1 getVarForMu {mu}" )
        # print ( f"@@1 nll={self.negative_log_likelihood ( mu )}" )
        dx = 1e-3
        hessian = np.diff ( np.diff ( [ self.negative_log_likelihood(x)["nll_obs_1"] for x in np.arange ( mu - 3e-3, mu + 3e-3, dx ) ] ) )
        h = hessian[hessian>0.] # if only some are negative, remove them
        if len(h)==0: # if all are negative, invert them
            h = np.abs ( hessian )
        h = np.mean ( h ) / dx /dx
        #if h < 0.:
        #    print ( f"@@X hessians {h} {hessian}" )
        return 1./h

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
        fmh = self.lmax(expected=expected, allowNegativeSignals=allowNegativeSignals, modelToUse = modelToUse )
        mu_hat, sigma_mu, _ = fmh["muhat"], fmh["sigma_mu"], fmh["lmax"]
        mu_hat = mu_hat if mu_hat is not None else 0.0
        nll0 = self.likelihood(mu_hat, expected=expected, return_nll=True,
                modelToUse = modelToUse )
        # a posteriori expected is needed here
        # mu_hat is mu_hat for signal_rel
        fmh = self.lmax(expected="posteriori", allowNegativeSignals=allowNegativeSignals,
                             return_nll=True, modelToUse = modelToUse )
        _, _, nll0A = fmh["muhat"], fmh["sigma_mu"], fmh["lmax"]

        # logger.error ( f"COMB nll0A {nll0A:.3f} mu_hatA {mu_hatA:.3f}" )
        # return 1.

        def clsRootTevatron( mu: float, return_type: Text = "CLs-alpha",
                     modelToUse : Union[None,str] = None ) -> float:
            # at - infinity this should be .95,
            # at + infinity it should -.05
            # Make sure to always compute the correct llhd value (from theoryPrediction)
            # and not used the cached value (which is constant for mu~=1 an mu~=0)
            nll_b = self.likelihood(0., return_nll=True, expected=expected, modelToUse = modelToUse )
            if nll_b is None:
                return None
            sigma2_b = self.getVarForMu ( 0. )
            nll_sb = self.likelihood(mu, return_nll=True, expected=expected, modelToUse = modelToUse )
            if nll_sb is None:
                return None
            sigma2_sb = self.getVarForMu ( 1. )
            CLs = CLsfromLsb(nll_sb, nll_b, sigma2_sb, sigma2_b, return_type=return_type)
            return CLs

        def clsRootAsimov( mu: float, return_type: Text = "CLs-alpha",
                     modelToUse : Union[None,str] = None ) -> float:
            # at - infinity this should be .95,
            # at + infinity it should -.05
            # Make sure to always compute the correct llhd value (from theoryPrediction)
            # and not used the cached value (which is constant for mu~=1 an mu~=0)
            nll = self.likelihood(mu, return_nll=True, expected=expected, modelToUse = modelToUse )
            nllA = self.likelihood(mu, expected="posteriori", return_nll=True, modelToUse = modelToUse )
            return CLsfromNLL(nllA, nll0A, nll, nll0, return_type=return_type) if nll is not None else None


        # return mu_hat, sigma_mu, clsRootTevatron
        return mu_hat, sigma_mu, clsRootAsimov


    def transform ( self, expected : Union [ Text, bool ] ):
        """ replace the actual observations with backgrounds,
            if expected is True or "posteriori" """
        # always start from scratch
        # self.model = copy.deepcopy ( self.origModel )
        if expected == False:
            return
        self.data.observed = self.model.backgrounds
