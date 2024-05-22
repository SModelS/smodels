#!/usr/bin/env python3

"""
.. module:: nnInterface
   :synopsis: Code that delegates the computation of limits and likelihoods to machine learned models

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from typing import Union, Text
import copy
import onnxruntime
from smodels.base.smodelsLogging import logger

nninfo = {
    "hasgreeted": False,
}

class NNData:
    """
    Holds data for use in the machine learned models
    :ivar nsignals: signal predictions list divided into sublists, one for each json file
    :ivar modelFile: path to onnx model file
    """

    def __init__(self, nsignals, globalInfo ):
        self.nsignals = nsignals  # fb
        self.getTotalYield()
        self.globalInfo = globalInfo
        #filename = globalInfo.modelFile
        #filename = os.path.join ( os.path.dirname ( globalInfo.path ), filename )
        #self.modelFile = modelFile
        # print ( f"@@0 modelFile {hasattr(globalInfo,'onnx')}" )
        self.cached_likelihoods = {}  ## cache of likelihoods (actually twice_nlls)
        self.cached_lmaxes = {}  # cache of lmaxes (actually twice_nlls)
        self.cachedULs = {False: {}, True: {}, "posteriori": {}}

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

    def welcome(self):
        """
        greet the world
        """

        if nninfo["hasgreeted"]:
            return
        logger.info( f"NN interface, we are using xxx" )
        nninfo["hasgreeted"] = True

    def likelihood( self, mu=1.0, workspace_index=None, return_nll=False,
                    expected=False):
        """
        Returns the value of the likelihood. \
        Inspired by the 'pyhf.infer.mle' module but for non-log likelihood

        :param workspace_index: supply index of workspace to use. If None, \
                                choose index of best combo
        :param return_nll: if true, return nll, not llhd
        :param expected: if False, compute expected values, if True, \
            compute a priori expected, if "posteriori" compute posteriori \
            expected
        """

        logger.debug("Calling likelihood")
        return 1.

    def exponentiateNLL(self, twice_nll, doIt):
        """if doIt, then compute likelihood from nll,
        else return nll"""
        if twice_nll == None:
            return None
            #if doIt:
            #    return 0.0
            #return 9000.0
        if doIt:
            return np.exp(-twice_nll / 2.0)
        return twice_nll / 2.0

    def lmax( self, workspace_index=None, return_nll=False, expected=False,
              allowNegativeSignals=False):
        """
        Returns the negative log max likelihood

        :param return_nll: if true, return nll, not llhd
        :param workspace_index: supply index of workspace to use. If None, \
            choose index of best combo
        :param expected: if False, compute expected values, if True, \
            compute a priori expected, if "posteriori" compute posteriori \
            expected
        :param allowNegativeSignals: if False, then negative nsigs are replaced \
            with 0.
        """
        # logger.error("expected flag needs to be heeded!!!")
        logger.debug("Calling lmax")
        lmax, muhat, sigma_mu = 1.0, 1.0, 1.0
        ret = { "lmax": lmax, "muhat": muhat, "sigma_mu": sigma_mu }
        return ret

    def getUpperLimitOnSigmaTimesEff(self, expected=False, workspace_index=None):
        """
        Compute the upper limit on the fiducial cross section sigma times efficiency:
            - by default, the combination of the workspaces contained into self.workspaces
            - if workspace_index is specified, self.workspace[workspace_index]
              (useful for computation of the best upper limit)

        :param expected:  - if set to 'True': uses expected SM backgrounds as signals
                          - else: uses 'self.nsignals'
        :param workspace_index: - if different from 'None': index of the workspace to use
                                  for upper limit
                                - else: choose best combo
        :return: the upper limit on sigma times eff at 'self.cl' level (0.95 by default)
        """
        if self.data.totalYield == 0.:
            return None
        else:
            ul = self.getUpperLimitOnMu( expected=expected, workspace_index=workspace_index)
            if ul == None:
                return ul
            if self.lumi is None:
                logger.error(f"asked for upper limit on fiducial xsec, but no lumi given with the data")
                return ul
            xsec = self.data.totalYield / self.lumi
            return ul * xsec

    # Trying a new method for upper limit computation :
    # re-scaling the signal predictions so that mu falls in [0, 10] instead of
    # looking for mu bounds
    # Usage of the index allows for rescaling
    def getUpperLimitOnMu(self, expected=False, workspace_index=None):
        """
        Compute the upper limit on the signal strength modifier with:
            - by default, the combination of the workspaces contained into self.workspaces
            - if workspace_index is specified, self.workspace[workspace_index]
              (useful for computation of the best upper limit)

        :param expected:  - if set to 'True': uses expected SM backgrounds as signals
                          - else: uses 'self.nsignals'
        :param workspace_index: - if different from 'None': index of the workspace to use
                                  for upper limit
                                - else: choose best combo
        :return: the upper limit at 'self.cl' level (0.95 by default)
        """
        return 1.0

    def transform ( self, expected : Union [ Text, bool ] ):
        """ replace the actual observations with backgrounds,
            if expected is True or "posteriori" """
        # always start from scratch
        # self.model = copy.deepcopy ( self.origModel )
        if expected == False:
            return
        self.data.observed = self.model.backgrounds
