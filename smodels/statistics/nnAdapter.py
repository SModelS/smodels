#!/usr/bin/env python3

"""
.. module:: nnAdapter
   :synopsis: An Adapter class that wraps around the neural networks (e.g. onnx
   files), handle all the pre and post processing. This adapter is
   meant to be published with the ML paper, not with e.g. SModelS.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

__all__ = [ "NNAdapter" ]

import os, onnx, json, math, onnxruntime, sys
import numpy as np
from typing import Union

class NNAdapter:
    """
    Adapter that wraps around a neural network
    """
    __slots__ = [ "allowsSyntheticData", "mlModel", "modelType",
                  "onnxMeta", "srOrder", "regressor" ]

    def __init__( self, mlModel : Union[bytes,str,onnx.ModelProto,os.PathLike],
                  allowsSyntheticData : bool = False ):
        """
        :param mlModel: the model, as a ModelProto, as a bytes stream,
        or as a path to an onnx file (needing to end with .onnx)
        :param allowsSyntheticData: if true, then also synthetic
        data can be supplied, not used yet
        """
        self.mlModel = mlModel
        if type(mlModel) == str and mlModel.endswith ( "onnx") and \
                os.path.exists ( mlModel ):
            self.mlModel = onnx.load ( mlModel )
        elif type(mlModel) in [ bytes, str ]:
            self.mlModel = onnx.load_model_from_string ( mlModel )
        self.allowsSyntheticData = allowsSyntheticData
        self.modelType = "rafal"
        self._parseMetaData ()
        self._getSROrder()
        self._instantiateRegressor()

    def _instantiateRegressor ( self ):
        """ create the actual inference session object """
        sess = onnxruntime.InferenceSession ( self.mlModel.SerializeToString() )
        self.regressor={ "session": sess,
                         "dim": sess.get_inputs()[0].shape[1] }

    def _fillValues ( self, container : str, values : str ) -> list:
        """ given <values> fill in <container>, if values are sensible
        :param container: the container to fill, e.g. nLL_obs_max
        :param values: the container values to copy from, e.g. [2,0]
        :returns: the values, but cast to proper type
        """
        tmp = json.loads(values)
        if type(tmp) in [ list, tuple ] and len(tmp)==2:
            if math.isfinite(tmp[1]):
                container = tmp
            if container[0]==None:
                container[0]=tmp[0]

        return container

    def _getSROrder ( self ):
        """ get the order of the signal regions as specified in the
        onnx meta information. We rely on smYields in the meta information
        to define the canonical order.
        """
        self.srOrder = []
        for srname in self.onnxMeta["smYields"]:
            self.srOrder.append ( srname )

    def _removeSignalRegions ( self, channels : list, dictionary : dict ) -> dict:
        """ remove a list of signal regions called "channels" 
        from the dictionary of values.
        :returns: pruned dictionary
        """
        newDict = {}
        for SRname,value in dictionary.items():
            if SRname in channels:
                continue
            p1 = SRname.rfind("-")
            if p1 > 0 and SRname[:p1] in channels:
                continue
            newDict[SRname]=value
        return newDict

    def _parseMetaData ( self ):
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
        data["featureMeans"]= []
        data["featureErrors"]= []
        data["nllMeans"]= []
        data["nllErrors"]= []
        remove_channels=[]
        import json, math
        isJoaquins = False
        for em in self.mlModel.metadata_props:
            if "rafal::" in str(em.key):
                isJoaquins = True
        if isJoaquins:
            self.modelType = "joaquin"
        for em in self.mlModel.metadata_props:
            emkey = em.key.replace ( "rafal::", "" )
            # print ( f"@@XXX emkey {emkey} value {em.value}" )
            if emkey == "remove_channels":
                # remove these channels at the end, so that order does not matter
                remove_channels = eval(em.value)
                data["remove_channels"] = remove_channels
            elif emkey == "obs_yields":
                st = eval(em.value)
                for l in st: ## the sm yields are tuple of (name,value)
                    data["obsYields"][ l[0] ]= int ( l[1] )
            elif emkey == "bkg_yields":
                st = eval(em.value)
                for l in st: ## the sm yields are tuple of (name,value)
                    data["smYields"][ l[0] ] = l[1]
            elif emkey == "standardization_mean":
                data["inputMeans"] = eval(em.value)
            elif emkey == "standardization_std":
                data["inputErrors"] = eval(em.value)
            elif emkey  in [ 'nLL_exp_mu0', 'nLL_obs_mu0', 'nLLA_exp_mu0', \
                              'nLLA_obs_mu0' ]:
                data[emkey] = json.loads(em.value)
            elif emkey in [ 'nLL_exp_max', 'nLL_obs_max', 'nLLA_exp_max', \
                             'nLLA_obs_max', 'nLL_exp_mu0', ]:
                data[emkey] = self._fillValues ( emkey, em.value )
            elif emkey == 'y_min':
                values = json.loads(em.value)
                if len(values)<7:
                    print ( f"[nnAdapter] 'y_min' in {onnxFile} has only {len(values)} entries, need 7." )
                    import sys; sys.exit(-1)
                indices = { "nLLA_obs_max": -1, "nLLA_exp_max": -3,
                            "nLL_obs_max" : -5, "nLL_exp_max": -7 }
                for name,index in indices.items():
                    if data[name] != None:
                        data[name] = [None,values[index]]
            elif emkey == "standardization":
                # joaquins values
                values = eval(em.value)
                data["featureMeans"] = values["features_mean"][0]*3
                data["featureErrors"] = values["features_std"][0]*3
                data["nllMeans"] = values["nLLs_mean"][0]
                data["nllErrors"] = values["nLLs_std"][0]


                #values = json.loads(em.value)
        if len(remove_channels)>0:
            data["smYields"]=self._removeSignalRegions ( remove_channels,
                                                   data["smYields"] )
            data["obsYields"]=self._removeSignalRegions ( remove_channels,
                                                    data["obsYields"] )
        self.onnxMeta={}
        for key,value in data.items():
            self.onnxMeta[key]=value

    def _scaleYieldsRafal ( self, yields : list ) -> np.array:
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

    def _scaleYieldsJoaquin ( self, yields : list ) -> np.array:
        """ scale the (total) yields Joaquin version

        :returns: the scaled total yields
        """
        scaled_yields = np.array( [yields], dtype=np.float32 )
        for i,x in enumerate(scaled_yields[0]):
            t = 0. # x
            err = self.onnxMeta["featureErrors"][i]
            if err > 1e-20:
                t = (x - self.onnxMeta["featureMeans"][i])/err
            #else:
            #    t = # - self.data.globalInfo.inputMeans[i]
            scaled_yields[0][i]=t
        return scaled_yields

    def _predictFromScaledYields ( self, scaled_yields : np.array ) -> np.array:
        """ get the prediction from the NN

        :param scaled_yields: the input of the neural network
        :returns: arr, the unscaled unshifted output of the neural network
        """
        if len(scaled_yields[0])!=self.regressor["dim"]:
            dim_nn = self.regressor["dim"]
            dim_input = len(scaled_yields[0])
            line=f"the network wants {dim_nn} input dimensions, but we supply {dim_input}. fix it!"
            print ( f"[nnAdapter] {line}" )
            sys.exit()
        if False:
            print ( f"@@NNX we evaluate at {scaled_yields}" )
            print ( f"@@NNX we evaluate for {self.regressor['dim']} dims" )
        dct = { "input_1": scaled_yields }
        if self.modelType == "joaquin":
            dct = { "features": scaled_yields }
        arr = self.regressor["session"].run(None, dct )
        arr = arr[0][0]
        return arr

    def _log_with_negatives ( self, x : np.array ) -> np.array:
        """ a pre-processing step that joaquin is performing 
        but rafal is not
        """
        return np.sign(x) * np.log1p ( np.abs(x) )

    def _undo_log_with_negatives ( self, x : np.array ) -> np.array:
        """ a post-processing step that joaquin is performing 
        but rafal is not
        """
        return np.sign(x) * np.expm1 ( np.abs(x) )

    def postprocess ( self, arr ) -> dict:
        if self.modelType == "joaquin":
            return self.postprocessJoaquin ( arr )
        return self.postprocessRafal ( arr )

    def postprocessRafal( self, arr ) -> dict:
        """ given the networks predictions, compute the NLLs

        :param arr: the neural network output
        :returns: { "nll_exp_0": ..., "nll_exp_1": ...,
                "nll_obs_0": ..., "nll_obs_1": ...,
                "nllA_exp_0": ..., "nllA_exp_1": ...,
                "nllA_obs_0": ..., "nllA_obs_1": ... }
        """
        #print ( f"@@postprocessRafal beginning arr {arr}" )
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
        # print ( f"@@postprocessRafal ret {ret}" )
        return ret

    def postprocessJoaquin( self, arr ) -> dict:
        """ given the networks predictions, compute the NLLs

        :param arr: the neural network output
        :returns: { "nll_exp_0": ..., "nll_exp_1": ...,
                "nll_obs_0": ..., "nll_obs_1": ...,
                "nllA_exp_0": ..., "nllA_exp_1": ...,
                "nllA_obs_0": ..., "nllA_obs_1": ... }
        """
        # print ( f"@@postprocessJoaquin: we start with {arr}" )
        # print ( f"@@postprocessJoaquin onnxMeta {self.onnxMeta.keys()}" )
        nll0obs =  self.onnxMeta["nLL_obs_mu0"]
        nll0exp =  self.onnxMeta["nLL_exp_mu0"]
        nllA0obs =  self.onnxMeta["nLLA_obs_mu0"]
        nllA0exp =  self.onnxMeta["nLLA_exp_mu0"]
        i_exp, i_obs, i_expA, i_obsA = -4, -3, -2, -1 # the indices
        # print ( "@@postprocessJoaquin nllMeans", self.onnxMeta["nllMeans"] )
        # print ( "@@postprocessJoaquin nllErrors", self.onnxMeta["nllErrors"] )
        d_nll1obs = self.onnxMeta["nllMeans"][i_obs] + self.onnxMeta["nllErrors"][i_obs]*arr[i_obs-4]
        d_nll1exp = self.onnxMeta["nllMeans"][i_exp] + self.onnxMeta["nllErrors"][i_exp]*arr[i_exp-4]
        d_nllA1obs = self.onnxMeta["nllMeans"][i_obsA] + self.onnxMeta["nllErrors"][i_obsA]*arr[i_obsA-4]
        d_nllA1exp = self.onnxMeta["nllMeans"][i_expA] + self.onnxMeta["nllErrors"][i_expA]*arr[i_expA-4]
        nll1obs = self._undo_log_with_negatives ( d_nll1obs ) + nll0obs
        nll1exp = self._undo_log_with_negatives ( d_nll1exp ) + nll0exp
        nllA1obs = self._undo_log_with_negatives ( d_nllA1obs ) + nllA0obs
        nllA1exp = self._undo_log_with_negatives ( d_nllA1exp ) + nllA0exp

        ret = { "nll_exp_0": nll0exp, "nll_exp_1": float(nll1exp),
                "nll_obs_0": nll0obs, "nll_obs_1": float(nll1obs),
                "nllA_exp_0": nllA0exp, "nllA_exp_1": float(nllA1exp),
                "nllA_obs_0": nllA0obs, "nllA_obs_1": float(nllA1obs) }
        ret["nll_obs_max"] = self.onnxMeta["nLL_obs_max"][1]
        # print ( f"@@postprocessJoaquin we return {ret}" )
        return ret

    def predict ( self, yields : Union[dict,list] ) -> dict:
        """ predict for yields, the main method
        :param yields: e.g. { "SR1": 3, "SR2": 5 }, or [3,5]
        (in which case the order must match the one in the json)

        :returns: { 'nll_exp_0': ..., 'nll_exp_1': ..., 'nll_obs_0': ...,
                    'nll_obs_1': ..., 'nllA_exp_0': ..., 'nllA_exp_1': ...,
                    'nllA_obs_0': ..., 'nllA_obs_1': ... }
        """
        # print ( f"@@predict for yields {yields}" )
        scaled_yields = self.preprocess ( yields )
        out = self._predictFromScaledYields ( scaled_yields )
        ret = self.postprocess ( out )
        #print ( f"@@predict ret {ret}" )
        return ret

    def preprocess ( self, yields : Union[dict,list] ) -> dict:
        """ perform all the preprocessing steps
        :param yields: e.g. { "SR1": 3, "SR2": 5 }, or [3,5]
        (in which case the order must match the one in the json)
        """
        if self.modelType == "joaquin":
            return self.preprocessJoaquin ( yields )
        return self.preprocessRafal ( yields )

    def preprocessJoaquinOld ( self, yields : Union[dict,list] ) -> dict:
        ## this was the previous version of preprocessing
        inp_list = yields
        if type(inp_list)==dict:
            inp_list = self._inputDictToList ( yields )
        # print ( f"@@preprocessJoaquin predict {yields} inp_list {inp_list}" )
        # print ( "@@postprocessJoaquin featureMeans", self.onnxMeta["featureMeans"] )
        inp_list = self._log_with_negatives ( inp_list )
        scaled_yields = self._scaleYieldsJoaquin ( inp_list )
        return scaled_yields

    def preprocessJoaquin ( self, yields : Union[dict,list] ) -> dict:
        print ( f"@@0 we need to process {yields}" )
        inp_list = yields
        if type(inp_list)==dict:
            inp_list = self._inputDictToList ( yields )
        # print ( f"@@preprocessJoaquin predict {yields} inp_list {inp_list}" )
        # print ( "@@postprocessJoaquin featureMeans", self.onnxMeta["featureMeans"] )
        inp_list = self._log_with_negatives ( inp_list )
        scaled_yields = self._scaleYieldsJoaquin ( inp_list )
        print ( f"@@1 the old code gave {scaled_yields}" )
        from smodels.statistics.joaquinsPreprocessing import preprocess_features, \
            log_with_negatives
        return_dict = True
        incl_fvs = True
        type_tokens = None
        trafos = None
        trafos = { "fvs_standardized": [
                "log_w_negatives", "standardization" ],
            "nLL_trafos": [ "log_w_negatives", "standardization"
            ]
        }
        re = preprocess_features ( yields, incl_fvs = incl_fvs,
            type_tokens = type_tokens, trafos = trafos,
            mean = self.onnxMeta["inputMeans"],
            std = self.onnxMeta["inputErrors"], return_dict = return_dict )
        print ( f"@@2 joaquins code gives {re}" )

        sys.exit()

    def preprocessRafal ( self, yields : Union[dict,list] ) -> dict:
        inp_list = yields
        if type(inp_list)==dict:
            inp_list = self._inputDictToList ( yields )
        # print ( f"@@preprocessRafal predict {yields} inp_list {inp_list}" )
        scaled_yields = self._scaleYieldsRafal ( inp_list )
        return scaled_yields

    def _inputDictToList ( self, in_dict : dict ) -> list:
        """ translate a dictionary of input yields to a list
        of said yields in the canonical order specified in the onnx
        :raises: exception if input SR is missing or one too many
        :returns: list of yields
        """
        ret = []
        if len(in_dict) != len ( self.srOrder ):
            raise Exception ( f"length of dict ({len(in_dict)} does not match srOrder ({len(self.srOrder)})" )
        for sr in self.srOrder:
            dsr = sr
            if dsr.endswith ( "-0" ):
                dsr = sr[:-2]
            if sr in in_dict:
                ret.append ( in_dict[sr] )
                continue
            if dsr in in_dict:
                ret.append ( in_dict[dsr] )
                continue
            print( f"signal region {sr} not in input_dict" )
            ret.append ( 0. )
        return ret

if __name__ == "__main__":
    regions = [ 'SRhigh_0Jb_cuts', 'SRhigh_0Jc_cuts', 'SRhigh_0Jd_cuts',
        'SRhigh_0Je_cuts', 'SRhigh_0Jf1_cuts', 'SRhigh_0Jf2_cuts',
        'SRhigh_0Jg1_cuts', 'SRhigh_0Jg2_cuts', 'SRhigh_nJa_cuts',
        'SRhigh_nJb_cuts', 'SRhigh_nJc_cuts', 'SRhigh_nJd_cuts',
        'SRhigh_nJe_cuts', 'SRhigh_nJf_cuts', 'SRhigh_nJg_cuts',
        'SRlow_0Jb_cuts', 'SRlow_0Jc_cuts', 'SRlow_0Jd_cuts',
        'SRlow_0Je_cuts', 'SRlow_0Jf1_cuts', 'SRlow_0Jf2_cuts',
        'SRlow_0Jg1_cuts', 'SRlow_0Jg2_cuts', 'SRlow_nJb_cuts',
        'SRlow_nJc_cuts', 'SRlow_nJd_cuts', 'SRlow_nJe_cuts',
        'SRlow_nJf1_cuts', 'SRlow_nJf2_cuts', 'SRlow_nJg1_cuts',
        'SRlow_nJg2_cuts', 'CR_0J_WZ_cuts', 'CR_nJ_WZ_cuts' ]
    # onnxFile = "../../unittests/testFiles/test.onnx"
    onnxFile = "test.onnx"

    adapter = NNAdapter ( onnxFile, False )

    yields = {}
    for region in regions: # predict for no yields
        yields[ region ] = 0.
    ret = adapter.predict ( yields )
    print ("\n".join( f"{key:10s}: {value:.1f}" for key,value in ret.items()))
    import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()
