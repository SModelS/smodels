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
from smodels.statistics.basicStats import determineBrentBracket, CLsfromNLL
from scipy import optimize

nninfo = {
    "hasgreeted": False,
}

class NNData:
    """
    Holds data for use in the machine learned models
    :ivar nsignals: signal predictions list divided into sublists, one for each json file
    :ivar modelFile: path to onnx model file
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
        S = sum ( self.nsignals )
        self.totalYield = S
        self.zeroSignalsFlag = self.totalYield == 0.

class NNUpperLimitComputer:
    """
    Class that computes the upper limit using the jsons files and signal informations in the 'data' instance of 'NNData'
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
        self.regressor = onnxruntime.InferenceSession ( self.data.globalInfo.onnx )
        # store the dimensionality of the input vector that the model
        # asks us to. this may be different from the number of our 
        # signal regions, as control regions may have been added.
        # we will pad with zeroes
        self.regressor_dim = self.regressor.get_inputs()[0].shape[1]
        self.lumi = lumi
        self.nsignals = copy.deepcopy ( self.data.nsignals )
        logger.debug("Signals : {}".format(self.nsignals))
        self.zeroSignalsFlag = self.data.zeroSignalsFlag
        self.cl = cl
        self.sigma_mu = None
        self.alreadyBeenThere = (
            False  # boolean to detect wether self.signals has returned to an older value
        )
        self.welcome()

    def negative_log_likelihood(self, poi_test: float):
        """ the method that really wraps around the llhd computation.

        :returns: dictionary with nlls, obs and exp, mu=0 and 1
        """
        # print ( f"@@0 datasetOrder {self.data.globalInfo.datasetOrderForModel}" )
        # print ( f"@@1 origDataSetOrder {self.data.origDataSetOrder}" )
        # print ( f"@@2 nsignals {self.nsignals} poi={poi_test} " )
        syields = []
        for ds in self.data.globalInfo.datasetOrderForModel.split(","):
            idx = self.data.origDataSetOrder.index ( ds )    
            tmp = float ( self.nsignals[idx]*poi_test )
            syields.append ( tmp )
        # syields = (np.array(self.nsignals)*poi_test).tolist()
        nzeroes = self.regressor_dim - len(syields)
        if nzeroes > 0:
            syields += [0]*nzeroes
        # print ( f"@@5 syields {syields}" )
        scaled_signal_yields = np.array( [syields], dtype=np.float32 )
        # print ( f"@@6 scaled_signal_yields {scaled_signal_yields}" )
        # print ( f"@@8 input {self.regressor.get_inputs()[0].shape}" )
        arr = self.regressor.run(None, {"input_1":scaled_signal_yields})
        arr = arr[0][0]
        # nLL_exp_mu0,nLL_exp_mu1,nLL_obs_mu0,nLL_obs_mu1
        ret = { "nll_exp_0": arr[0], "nll_exp_1": arr[1],
                "nll_obs_0": arr[2], "nll_obs_1": arr[3] }
        return ret

    def welcome(self):
        """
        greet the world
        """

        if nninfo["hasgreeted"]:
            return
        logger.info( f"NN interface, we are using xxx" )
        nninfo["hasgreeted"] = True

    def likelihood( self, mu=1.0, return_nll=False, expected=False):
        """
        Returns the value of the likelihood. \
        Inspired by the 'pyhf.infer.mle' module but for non-log likelihood

        :param return_nll: if true, return nll, not llhd
        :param expected: if False, compute expected values, if True, \
            compute a priori expected, if "posteriori" compute posteriori \
            expected
        """
        ret = self.negative_log_likelihood ( mu )
        lbl = "nll_obs_1"
        if expected:
            lbl = "nll_exp_1"
        nll = ret[lbl]
        logger.debug( f"Calling likelihood")
        return self.exponentiateNLL ( nll, not return_nll )

    def exponentiateNLL(self, nll, doIt):
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
              allowNegativeSignals=False):
        """
        Returns the negative log max likelihood

        :param return_nll: if true, return nll, not llhd
        :param expected: if False, compute expected values, if True,
        compute a priori expected, if "posteriori" compute posteriori
        expected
        :param allowNegativeSignals: if False, then negative nsigs are 
        replaced with 0.
        """
        # logger.error("expected flag needs to be heeded!!!")
        lmax, muhat, sigma_mu = float("nan"),float("nan"),float("nan")
        logger.error("Calling lmax")
        # print ( f"@@5 before minimize!" )
        def getNLL ( mu ):
            d = self.negative_log_likelihood ( mu )
            lbl = "nll_obs_1"
            if expected:
                lbl = "nll_exp_1"
            ret = d[lbl]
            # print ( f"@@X getNLL {mu}={ret}" )
            return ret
        for mu0 in [ 1.0, 0.0, 3.0, -1.0, 10.0, 0.1 ]:
            o = optimize.minimize ( getNLL, mu0, tol=1e-9 )
            # print ( f"@@6 o={o}" )
            if o.success == True:
                muhat = o.x
                lmax = self.exponentiateNLL ( o.fun, not return_nll )
                if not allowNegativeSignals and muhat < 0.:
                    muhat = 0.
                    lmax = self.likelihood ( 0., return_nll = return_nll,
                        expected = expected )
                sigma_mu = np.sqrt ( o.hess_inv[0][0] )
                ret = { "lmax": lmax, "muhat": muhat, 
                        "sigma_mu": sigma_mu }
                # print ( f"@@7 ret={ret}" )
                #sys.exit()
                return ret
        ret = { "lmax": lmax, "muhat": muhat, "sigma_mu": sigma_mu }
        #print ( f"@@9 ret {ret}" )
        return ret

    def getUpperLimitOnSigmaTimesEff(self, expected=False ):
        """
        Compute the upper limit on the fiducial cross section sigma times efficiency:

        :param expected:  - if set to 'True': uses expected SM backgrounds as signals
                          - else: uses 'self.nsignals'
        :return: the upper limit on sigma times eff at 'self.cl' level (0.95 by default)
        """
        if self.data.totalYield == 0.:
            return None
        else:
            ul = self.getUpperLimitOnMu( expected=expected )
            if ul == None:
                return ul
            if self.lumi is None:
                logger.error(f"asked for upper limit on fiducial xsec, but no lumi given with the data")
                return ul
            xsec = self.data.totalYield / self.lumi
            return ul * xsec

    def getUpperLimitOnMu(self, expected=False, 
            allowNegativeSignals = False ):
        """
        Compute the upper limit on the signal strength modifier with:
            - by default, the combination of the workspaces contained into self.workspaces

        :param expected:  - if set to 'True': uses expected SM backgrounds as signals
                          - else: uses 'self.nsignals'
        :return: the upper limit at 'self.cl' level (0.95 by default)
        """
        # print ( f"@@5 getUpperLimitOnMu" )
        mu_hat, sigma_mu, clsRoot = self.getCLsRootFunc(expected=expected,
                              allowNegativeSignals=allowNegativeSignals)
        a, b = determineBrentBracket(mu_hat, sigma_mu, clsRoot,
                allowNegative = allowNegativeSignals )
        mu_lim = optimize.brentq(clsRoot, a, b, rtol=1e-03, xtol=1e-06)
        return mu_lim

    def getCLsRootFunc(self, expected: bool = False, allowNegativeSignals : bool = False) -> Tuple[float, float, Callable]:
        """
        Obtain the function "CLs-alpha[0.05]" whose root defines the upper limit,
        plus mu_hat and sigma_mu

        :param expected: if True, compute expected likelihood, else observed
        """
        fmh = self.lmax(expected=expected, allowNegativeSignals=allowNegativeSignals)
        mu_hat, sigma_mu, _ = fmh["muhat"], fmh["sigma_mu"], fmh["lmax"]
        mu_hat = mu_hat if mu_hat is not None else 0.0
        nll0 = self.likelihood(mu_hat, expected=expected, return_nll=True)
        # a posteriori expected is needed here
        # mu_hat is mu_hat for signal_rel
        fmh = self.lmax(expected="posteriori", allowNegativeSignals=allowNegativeSignals,
                             return_nll=True)
        _, _, nll0A = fmh["muhat"], fmh["sigma_mu"], fmh["lmax"]

        # logger.error ( f"COMB nll0A {nll0A:.3f} mu_hatA {mu_hatA:.3f}" )
        # return 1.

        def clsRoot(mu: float, return_type: Text = "CLs-alpha") -> float:
            # at - infinity this should be .95,
            # at + infinity it should -.05
            # Make sure to always compute the correct llhd value (from theoryPrediction)
            # and not used the cached value (which is constant for mu~=1 an mu~=0)
            nll = self.likelihood(mu, return_nll=True, expected=expected)
            nllA = self.likelihood(mu, expected="posteriori", return_nll=True)
            return CLsfromNLL(nllA, nll0A, nll, nll0, return_type=return_type) if nll is not None else None

        return mu_hat, sigma_mu, clsRoot


    def transform ( self, expected : Union [ Text, bool ] ):
        """ replace the actual observations with backgrounds,
            if expected is True or "posteriori" """
        # always start from scratch
        # self.model = copy.deepcopy ( self.origModel )
        if expected == False:
            return
        self.data.observed = self.model.backgrounds
