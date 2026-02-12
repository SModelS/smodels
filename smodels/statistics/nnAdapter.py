#!/usr/bin/env python3

"""
.. module:: nnAdapter
   :synopsis: An Adapter class that wraps around the neural networks (e.g. onnx
   files), handle all the pre and post processing. This adapter is
   meant to be published with the ML paper, and the SModels database, not with
   SModelS per se.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os, onnx, json, math, onnxruntime
import numpy as np
from typing import Union

class NNAdapter:
    """
    Adapter that wraps around a neural network
    """

    def __init__( self, onnxfile : os.PathLike,
                  allowsSyntheticData : bool = False ):
        """
        :param onnxfile: path to onnx file
        to the ML model's SRs
        :param allowsSyntheticData: if true, then also synthetic
        data can be supplied
        """
        self.onnxfile = onnxfile
        self.allowsSyntheticData = allowsSyntheticData
        self.mlModel = onnx.load ( onnxfile )
        self.parseMetaData ()
        self.getSROrder()
        self.instantiateRegressor()

    def instantiateRegressor ( self ):
        """ create the actual inference session object """
        sess = onnxruntime.InferenceSession ( self.onnxfile )
        self.regressor={ "session": sess,
                         "dim": sess.get_inputs()[0].shape[1] }

    def fillValues ( self, container, value ):
        """ given <value> fill in <container>, if value is sensible
        :param container: the container to fill
        :param value: the container, value to copy from
        """
        tmp = json.loads(value)
        if type(tmp) in [ list, tuple ] and len(tmp)==2:
            if math.isfinite(tmp[1]):
                container = tmp
            if container[0]==None:
                container[0]=tmp[0]


    def getSROrder ( self ):
        """ get the order of the signal regions as specified in the
        onnx meta information """
        self.srOrder = []
        for srname in self.onnxMeta["smYields"]:
            self.srOrder.append ( srname )

    def parseMetaData ( self ):
        """ parse the model's meta data """
        data = { "inputMeans": [], "inputErrors": [],
            "nLL_exp_mu0": [ None ]*2, "nLL_obs_mu0": [ None ]*2,
            "nLLA_exp_mu0": [ None ]*2, "nLLA_obs_mu0": [ None ]*2 }
        data [ "smYields" ] = {}
        data [ "obsYields" ] = {}
        data["nLL_obs_max"]= [ None ] * 2
        data["nLL_exp_max"]= [ None ] * 2
        data["nLLA_obs_max"]= [ None ] * 2
        data["nLLA_exp_max"]= [ None ] * 2
        remove_channels=[]
        import json, math
        for em in self.mlModel.metadata_props:
            if em.key == "remove_channels":
                # remove these channels at the end, so that order does not matter
                remove_channels = eval(em.value)
                data["remove_channels"] = remove_channels
            if em.key == "obs_yields":
                st = eval(em.value)
                for l in st: ## the sm yields are tuple of (name,value)
                    data["obsYields"][ l[0] ]= int ( l[1] )
            if em.key == "bkg_yields":
                st = eval(em.value)
                for l in st: ## the sm yields are tuple of (name,value)
                    data["smYields"][ l[0] ] = l[1]
            if em.key == "standardization_mean":
                data["inputMeans"] = eval(em.value)
            elif em.key == "standardization_std":
                data["inputErrors"] = eval(em.value)
            elif em.key  in [ 'nLL_exp_mu0', 'nLL_obs_mu0', 'nLLA_exp_mu0', \
                              'nLLA_obs_mu0' ]:
                data[em.key] = json.loads(em.value)
            elif em.key in [ 'nLL_exp_max', 'nLL_obs_max', 'nLLA_exp_max', \
                             'nLLA_obs_max', 'nLL_exp_mu0', ]:
                self.fillValues ( em.key, em.value )
            elif em.key == 'y_min':
                values = json.loads(em.value)
                if len(values)<7:
                    logger.error ( f"'y_min' in {onnxFile} has only {len(values)} entries, need 7." )
                    import sys; sys.exit(-1)
                indices = { "nLLA_obs_max": -1, "nLLA_exp_max": -3,
                            "nLL_obs_max" : -5, "nLL_exp_max": -7 }
                for name,index in indices.items():
                    if data[name] != None:
                        data[name] = [None,values[index]]
        if len(remove_channels)>0:
            data["smYields"]=removeSignalRegions ( remove_channels, data["smYields"] )
            data["obsYields"]=removeSignalRegions ( remove_channels, data["obsYields"] )
        # ['smYields', 'obsYields', 'inputMeans', 'inputErrors' ]
        self.onnxMeta={}
        for key,value in data.items():
            self.onnxMeta[key]=value

    def scaleYields ( self, yields : list ) -> np.array:
        """ scale the (total) yields

        :returns: the scaled total yields
        """
        scaled_yields = np.array( [yields], dtype=np.float32 )
        for i,x in enumerate(scaled_yields[0]):
            t = 0. # x
            err = self.onnxMeta["inputErrors"][i]
            if err > 1e-20:
                t = (x - self.onnxMeta["inputMeans"][i])/err
            #else:
            #    t = # - self.data.globalInfo.inputMeans[i]
            scaled_yields[0][i]=t
        return scaled_yields

    def predictFromScaledYields ( self, scaled_yields : np.array ):
        """ get the prediction from the NN

        :param scaled_yields: the input of the neural network
        :returns: arr, the unscaled unshifted output of the neural network
        """
        #if poi_test == 0.:
        #    print ( f"@@NNX we evaluate at {scaled_yields}" )
        if len(scaled_yields[0])!=self.regressor["dim"]:
            dim_nn = self.regressor["dim"]
            dim_input = len(scaled_yields[0])
            line=f"the network wants {dim_nn} input dimensions, but we supply {dim_input}. fix it!"
            logger.error ( f"[nnInterface] {line}" )
            print ( f"[nnInterface] {line}" )
            sys.exit()
        arr = self.regressor["session"].run(None,
                {"input_1":scaled_yields})
        # print ( f"@@NNA arr {arr}" )
        arr = arr[0][0]
        return arr

    def nllsFromPrediction( self, arr ) -> dict:
        """ given the networks predictions, compute the NLLs

        :param arr: the neural network output
        :returns: { "nll_exp_0": ..., "nll_exp_1": ...,
                "nll_obs_0": ..., "nll_obs_1": ...,
                "nllA_exp_0": ..., "nllA_exp_1": ...,
                "nllA_obs_0": ..., "nllA_obs_1": ... }
        """
        nll0obs =  self.onnxMeta["nLL_obs_mu0"]
        nll0exp =  self.onnxMeta["nLL_exp_mu0"]
        nllA0obs =  self.onnxMeta["nLLA_obs_mu0"]
        nllA0exp =  self.onnxMeta["nLLA_exp_mu0"]
        i_exp, i_obs, i_expA, i_obsA = -4, -3, -2, -1 # the indices
        expDelta = self.onnxMeta["inputMeans"][i_exp]
        obsDelta = self.onnxMeta["inputMeans"][i_obs]
        expDeltaA = self.onnxMeta["inputMeans"][i_expA]
        obsDeltaA = self.onnxMeta["inputMeans"][i_obsA]
        expErr = self.onnxMeta["inputErrors"][i_exp]
        obsErr = self.onnxMeta["inputErrors"][i_obs]
        expErrA = self.onnxMeta["inputErrors"][i_expA]
        obsErrA = self.onnxMeta["inputErrors"][i_obsA]
        nll1exp = nll0exp + arr[i_exp]*expErr + expDelta
        nll1obs = nll0obs + arr[i_obs]*obsErr + obsDelta
        nllA1exp = nllA0exp + arr[i_expA]*expErrA + expDeltaA
        nllA1obs = nllA0obs + arr[i_obsA]*obsErrA + obsDeltaA

        ret = { "nll_exp_0": nll0exp, "nll_exp_1": float(nll1exp),
                "nll_obs_0": nll0obs, "nll_obs_1": float(nll1obs),
                "nllA_exp_0": nllA0exp, "nllA_exp_1": float(nllA1exp),
                "nllA_obs_0": nllA0obs, "nllA_obs_1": float(nllA1obs) }
        return ret

    def predict ( self, yields : dict ) -> dict:
        """ predict for yields
        :param yields: e.g. { "SR1": 3, "SR2": 5 }

        :returns: xxx
        """
        inp_list = self.inputDictToList ( yields )
        scaled_yields = self.scaleYields ( inp_list )
        out = self.predictFromScaledYields ( scaled_yields )
        return self.nllsFromPrediction ( out )

    def inputDictToList ( self, in_dict : dict ) -> list:
        """ translate a dictionary of input yields to a list
        of said yields in the canonical order specified in the onnx
        :raises: exception if input SR is missing or one too many
        :returns: list of yields
        """
        ret = []
        #if len(in_dict) != len ( self.srOrder ):
        #    raise Exception ( f"length of dict ({len(in_dict)} does not match srOrder ({len(self.srOrder)})" )
        for sr in self.srOrder:
            if not sr in in_dict:
                ret.append ( 0. )
                # raise Exception( f"signal region {sr} not in input_dict" )
            else:
                ret.append ( in_dict[sr] )
        return ret

    def interact ( self ):
        import IPython; IPython.embed( colors = "neutral" )

def getRegionsForExample():
    # the regions dictionary of the example
    ret = [
        {'smodels': 'SRhigh_0Jb', 'pyhf': 'SRhigh_0Jb_cuts'},
        {'smodels': 'SRhigh_0Jc', 'pyhf': 'SRhigh_0Jc_cuts'},
        {'smodels': 'SRhigh_0Jd', 'pyhf': 'SRhigh_0Jd_cuts'},
        {'smodels': 'SRhigh_0Je', 'pyhf': 'SRhigh_0Je_cuts'},
        {'smodels': 'SRhigh_0Jf1', 'pyhf': 'SRhigh_0Jf1_cuts'},
        {'smodels': 'SRhigh_0Jf2', 'pyhf': 'SRhigh_0Jf2_cuts'},
        {'smodels': 'SRhigh_0Jg1', 'pyhf': 'SRhigh_0Jg1_cuts'},
        {'smodels': 'SRhigh_0Jg2', 'pyhf': 'SRhigh_0Jg2_cuts'},
        {'smodels': 'SRhigh_nJa', 'pyhf': 'SRhigh_nJa_cuts'},
        {'smodels': 'SRhigh_nJb', 'pyhf': 'SRhigh_nJb_cuts'},
        {'smodels': 'SRhigh_nJc', 'pyhf': 'SRhigh_nJc_cuts'},
        {'smodels': 'SRhigh_nJd', 'pyhf': 'SRhigh_nJd_cuts'},
        {'smodels': 'SRhigh_nJe', 'pyhf': 'SRhigh_nJe_cuts'},
        {'smodels': 'SRhigh_nJf', 'pyhf': 'SRhigh_nJf_cuts'},
        {'smodels': 'SRhigh_nJg', 'pyhf': 'SRhigh_nJg_cuts'},
        {'smodels': 'SRlow_0Jb', 'pyhf': 'SRlow_0Jb_cuts'},
        {'smodels': 'SRlow_0Jc', 'pyhf': 'SRlow_0Jc_cuts'},
        {'smodels': 'SRlow_0Jd', 'pyhf': 'SRlow_0Jd_cuts'},
        {'smodels': 'SRlow_0Je', 'pyhf': 'SRlow_0Je_cuts'},
        {'smodels': 'SRlow_0Jf1', 'pyhf': 'SRlow_0Jf1_cuts'},
        {'smodels': 'SRlow_0Jf2', 'pyhf': 'SRlow_0Jf2_cuts'},
        {'smodels': 'SRlow_0Jg1', 'pyhf': 'SRlow_0Jg1_cuts'},
        {'smodels': 'SRlow_0Jg2', 'pyhf': 'SRlow_0Jg2_cuts'},
        {'smodels': 'SRlow_nJb', 'pyhf': 'SRlow_nJb_cuts'},
        {'smodels': 'SRlow_nJc', 'pyhf': 'SRlow_nJc_cuts'},
        {'smodels': 'SRlow_nJd', 'pyhf': 'SRlow_nJd_cuts'},
        {'smodels': 'SRlow_nJe', 'pyhf': 'SRlow_nJe_cuts'},
        {'smodels': 'SRlow_nJf1', 'pyhf': 'SRlow_nJf1_cuts'},
        {'smodels': 'SRlow_nJf2', 'pyhf': 'SRlow_nJf2_cuts'},
        {'smodels': 'SRlow_nJg1', 'pyhf': 'SRlow_nJg1_cuts'},
        {'smodels': 'SRlow_nJg2', 'pyhf': 'SRlow_nJg2_cuts'}]
    return ret

if __name__ == "__main__":
    onnxfile = "../../unittests/testFiles/test.onnx"
    regions = getRegionsForExample()
    adapter = NNAdapter ( onnxfile, False )
    yields = {}
    for region in regions:
        yields[ region["pyhf"] ] = 0.
    ret = adapter.predict ( yields )
    print ( ret )
    adapter.interact()
