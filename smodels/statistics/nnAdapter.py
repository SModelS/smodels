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
                  onnxfilename : str,
                  allowsSyntheticData : bool = False ):
        """
        :param mlModel: the model, as a ModelProto, as a bytes stream,
        or as a path to an onnx file (needing to end with .onnx)
        :param onnxfilename: filename of onnxfile, for debugging only
        :param allowsSyntheticData: if true, then also synthetic
        data can be supplied, not used yet
        """
        if type(mlModel) == str and mlModel.endswith ( "onnx") and \
                os.path.exists ( mlModel ):
            self.mlModel = onnx.load ( mlModel )
        elif type(mlModel) in [ bytes, str ]:
            try:
                self.mlModel = onnx.load_model_from_string ( mlModel )
            except Exception as e:
                print( f"[nnAdapter] could not load model {onnxfilename}" )
                sys.exit(-1)
        self.allowsSyntheticData = allowsSyntheticData
        self.modelType = "joaquin"
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
        for em in self.mlModel.metadata_props:
            if "rafal::" in str(em.key):
                self.modelType = "joaquin"
        for em in self.mlModel.metadata_props:
            emkey = em.key.replace ( "rafal::", "" )
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
            elif emkey == "run_config":
                import yaml
                content = yaml.safe_load ( em.value )
                data["run_config"] = content

                #values = json.loads(em.value)
        if len(remove_channels)>0:
            data["smYields"]=self._removeSignalRegions ( remove_channels,
                                                   data["smYields"] )
            data["obsYields"]=self._removeSignalRegions ( remove_channels,
                                                    data["obsYields"] )
        #self.onnxMeta={}
        #for key,value in data.items():
        #    self.onnxMeta[key]=value
        self.onnxMeta = data

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
            print ( f"[nnAdapter] srOrder: {self.srOrder}" )
            sys.exit()
        dct = { "features": scaled_yields }
        arr = self.regressor["session"].run(None, dct )
        arr = arr[0][0]
        return arr

    def postprocess( self, arr : np.ndarray,
           add_errors : bool = True ) -> dict:
        """ given the networks predictions, compute the NLLs

        :param arr: the neural network output
        :param add_errors: if true, then add errors also FIXME describe
        more
        :returns: { "nll_exp_0": ..., "nll_exp_1": ...,
                "nll_obs_0": ..., "nll_obs_1": ...,
                "nllA_exp_0": ..., "nllA_exp_1": ...,
                "nllA_obs_0": ..., "nllA_obs_1": ... }
        """
        from smodels.statistics.joaquinsPreprocessing import undo_preprocess_nLLs
        deltas_prepd = np.array(arr, dtype=np.float64)
        trafos = self.onnxMeta["run_config"]["data"]["nLL_trafos"]
        nll_means = self.onnxMeta["nllMeans"]
        nll_errors = self.onnxMeta["nllErrors"]
        deltas = undo_preprocess_nLLs ( deltas_prepd[:4],
                mean = nll_means, std = nll_errors, trafos = trafos )
        deltas = list ( map ( float, deltas ) )
        nll0exp  = self.onnxMeta["nLL_exp_mu0"]
        nll0obs  = self.onnxMeta["nLL_obs_mu0"]
        nllA0exp = self.onnxMeta["nLLA_exp_mu0"]
        nllA0obs = self.onnxMeta["nLLA_obs_mu0"]

        nll1exp  = nll0exp  + deltas[0]
        nll1obs  = nll0obs  + deltas[1]
        nllA1exp = nllA0exp + deltas[2]
        nllA1obs = nllA0obs + deltas[3]
        ## error propagation, fixme for now we just do it by hand:
        ## s_y = abs ( y * s_x )
        ## FIXME this needs to be changed for something generic
        #s_nll1exp  = abs ( deltas[4] * deltas[0] )
        #s_nll1obs  = abs ( deltas[5] * deltas[1] )
        #s_nllA1exp = abs ( deltas[6] * deltas[2] )
        #s_nllA1obs = abs ( deltas[7] * deltas[3] )

        ret = { "nll_exp_0": nll0exp,  "nll_exp_1": nll1exp,
                "nll_obs_0": nll0obs,  "nll_obs_1": nll1obs,
                "nllA_exp_0": nllA0exp, "nllA_exp_1": nllA1exp,
                "nllA_obs_0": nllA0obs, "nllA_obs_1": nllA1obs }
        if self.onnxMeta["nLL_obs_max"][1] is not None:
            ret["nll_obs_max"] = self.onnxMeta["nLL_obs_max"][1]
        if add_errors:
            from smodels.statistics.joaquinsPreprocessing import \
                undo_preprocess_nLLs_errors
            errs = undo_preprocess_nLLs_errors ( deltas_prepd[4:],
                    deltas_prepd[:4],
                    mean = nll_means,
                    std = nll_errors, trafos = trafos,
                    eps = 1e-5 )
            ret["sigma_exp"] = errs[0]
            ret["sigma_obs"] = errs[1]
            ret["sigma_expA"] = errs[2]
            ret["sigma_obsA"] = errs[3]
        return ret

    def predict ( self, yields : Union[dict,list] ) -> dict:
        """ predict for yields, the main method
        :param yields: e.g. { "SR1": 3, "SR2": 5 }, or [3,5]
        (in which case the order must match the one in the json)

        :returns: { 'nll_exp_0': ..., 'nll_exp_1': ..., 'nll_obs_0': ...,
                    'nll_obs_1': ..., 'nllA_exp_0': ..., 'nllA_exp_1': ...,
                    'nllA_obs_0': ..., 'nllA_obs_1': ... }
        """
        scaled_yields = self.preprocess ( yields )
        out = self._predictFromScaledYields ( scaled_yields )
        ret = self.postprocess ( out )
        return ret

    def preprocess ( self, yields : Union[dict,list] ) -> dict:
        inp_list = yields
        if type(inp_list)==dict:
            inp_list = self._inputDictToList ( yields )
        from smodels.statistics.joaquinsPreprocessing import preprocess_features
        trafos = self.onnxMeta["run_config"]["data"]["trafos"]
        nYields = len(yields)
        re = preprocess_features ( inp_list,
            trafos = trafos,
            mean = self.onnxMeta["featureMeans"][:nYields],
            std = self.onnxMeta["featureErrors"][:nYields] )
        ret = [ re[0] ]
        return ret

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
