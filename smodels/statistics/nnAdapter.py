#!/usr/bin/env python3

"""
.. module:: nnAdapter
   :synopsis: An Adapter class that wraps around the neural networks
   published in arXiv:XXXX, handling all the pre- and post-processing.

.. moduleauthor:: OLLL Collaboration

"""

__all__ = [ "NNAdapter" ]

import os
import onnx
import json
import math
import onnxruntime
import sys
import numpy as np
from typing import Union

class NNAdapter:
    """
    Adapter that wraps around a neural network
    """
    __slots__ = [ "mlModel", "modelType", "onnxMeta", "srOrder", "regressor",
                  "session_options", "onnxfilename", "crRegions" ]

    def __init__( self, mlModel : Union[bytes,str,onnx.ModelProto,os.PathLike],
                  onnxfilename : None|str = None, session_options : dict = {},
                  validate_metadata : bool = True ):
        """
        :param mlModel: the model, as a ModelProto, as a bytes stream,
        or as a path to an onnx file (needing to end with .onnx)
        :param onnxfilename: filename of onnxfile, for debugging only
        if None, then assume it is the same as mlModel
        :param session_options: options for the onnxruntime inference session,
        e.g. { "inter_op_num_threads": 1 }
        :param validate_metadata: if true, then validate metadata in onnx file,
        before usage
        """
        if onnxfilename == None:
            onnxfilename = str(mlModel)
        assert type(mlModel) in [ bytes, str, onnx.ModelProto,os.PathLike],\
            "mlModel needs to be one of: bytes, str, onnx.ModelProto, PathType"
        if type(mlModel) == str and mlModel.endswith ( "onnx") and \
                os.path.exists ( mlModel ):
            self.mlModel = onnx.load ( mlModel )
        elif type(mlModel) in [ bytes, str ]:
            try:
                self.mlModel = onnx.load_model_from_string ( mlModel )
            except Exception:
                print( f"[nnAdapter] could not load model {onnxfilename}" )
                sys.exit(-1)
        self.onnxfilename = onnxfilename
        self.session_options = session_options
        if validate_metadata:
            from smodels.statistics.metadataValidator import validateMetaData
            validateMetaData ( self.mlModel.metadata_props )
        self._parseMetaData ()
        self._getSROrder()
        self._cleanCRs()
        self._instantiateRegressor()

    def predict ( self, yields : Union[dict,list],
           yields_are_signal_yields : bool = True,
           obs_as_bg : list|None|str = [] ) -> dict:
        """ proposal for a slightly different API

        :param yields: e.g. { "SR1": 3, "SR2": 5 }, or [3,5]
        (in which case the order must match the one in the json)

        :param yields_are_signal_yields: if True, then yields are
        interpreted as signal yields, and the backgrounds get added.
        if False, yields are assumed to be total yields

        :param obs_as_bg: a list of signal regions for which we use 
        observations as background_yields ("postfit"), given 
        yields_are_signal_yields is True. If None or "default", then 
        use self.onnxMeta["crRegions"] as defined in the onnx file

        :returns: the negative log likelihoods (nlls) as a dictionary:
        { 'nll_exp_0': ..., 'nll_exp_1': ..., 'nll_obs_0': ...,
        'nll_obs_1': ..., 'nllA_exp_0': ..., 'nllA_exp_1': ...,
        'nllA_obs_0': ..., 'nllA_obs_1': ... }
        where 0, 1 means mu=0, 1, respectively. exp refers to a priori
        expectation, obs are the observed values. nllA means the 
        nll is evaluated for the Asimov dataset with mu' = 0.
        """
        if obs_as_bg in [ None, "default", "postfit" ]:
            obs_as_bg = self.onnxMeta["crRegions"]
        if yields_are_signal_yields:
            yields = self._totalYieldsFromSignals ( yields, obs_as_bg )
        scaled_yields = self._preprocess ( yields )
        out = self._predictFromScaledYields ( scaled_yields )
        ret = self._postprocess ( out )
        return ret

    def _getCRs( self, channels : list ) -> list:
        """ get a list of every signal region marked as a control region
        """
        crRegions = []
        for ch in channels:
            for regionName, regionType in ch.items():
                if regionType == "CR":
                    crRegions.append ( regionName )
        return crRegions

    def _cleanCRs ( self ):
        """ the meta information has all regions of all models,
        so we clean the list of control regions here, 
        possibly also adding the "-o" postfix to regio names
        """
        newCRs = []
        for r in self.onnxMeta["crRegions"]:
            if r in self.srOrder:
                newCRs.append ( r )
            if r+"-0" in self.srOrder:
                newCRs.append ( f"{r}-0" )
        self.onnxMeta["crRegions"] = newCRs

    def _instantiateRegressor ( self ):
        """ create the actual inference session object """
        so = onnxruntime.SessionOptions()
        for k,v in self.session_options.items():
            setattr ( so, k, v )
        sess = onnxruntime.InferenceSession ( self.mlModel.SerializeToString(),
               so )
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
        onnx meta information. We rely on bkg_yields in the meta information
        to define the canonical order.
        """
        self.srOrder = []
        for srname in self.onnxMeta["bkg_yields"]:
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
        data [ "bkg_yields" ] = {}
        data [ "obs_yields" ] = {}
        data["nLL_obs_max"]= [ None ] * 2
        data["nLL_exp_max"]= [ None ] * 2
        data["nLLA_obs_max"]= [ None ] * 2
        data["nLLA_exp_max"]= [ None ] * 2
        data["featureMeans"]= []
        data["featureErrors"]= []
        data["nllMeans"]= []
        data["nllErrors"]= []
        remove_channels=[]
        import json
        for em in self.mlModel.metadata_props:
            if em.key == "channels":
                data["crRegions"] = self._getCRs ( eval ( em.value ) )
            if em.key == "remove_channels":
                # remove these channels at the end, so that order does not matter
                remove_channels = eval(em.value)
                data["remove_channels"] = remove_channels
            elif em.key == "obs_yields":
                st = eval(em.value)
                for l in st: ## the sm yields are tuple of (name,value)
                    data["obs_yields"][ l[0] ]= int ( l[1] )
            elif em.key == "bkg_yields":
                st = eval(em.value)
                for l in st: ## the sm yields are tuple of (name,value)
                    data["bkg_yields"][ l[0] ] = l[1]
            elif em.key == "standardization_mean":
                data["inputMeans"] = eval(em.value)
            elif em.key == "standardization_std":
                data["inputErrors"] = eval(em.value)
            elif em.key  in [ 'nLL_exp_mu0', 'nLL_obs_mu0', 'nLLA_exp_mu0', \
                              'nLLA_obs_mu0' ]:
                data[em.key] = json.loads(em.value)
            elif em.key in [ 'nLL_exp_max', 'nLL_obs_max', 'nLLA_exp_max', \
                             'nLLA_obs_max', 'nLL_exp_mu0', ]:
                data[em.key] = self._fillValues ( em.key, em.value )
            elif em.key == "standardization":
                values = eval(em.value)
                data["featureMeans"] = values["features_mean"][0]*3
                data["featureErrors"] = values["features_std"][0]*3
                data["nllMeans"] = values["nLLs_mean"][0]
                data["nllErrors"] = values["nLLs_std"][0]
            elif em.key == "run_config":
                import yaml
                content = yaml.safe_load ( em.value )
                data["run_config"] = content

        if len(remove_channels)>0:
            data["bkg_yields"]=self._removeSignalRegions ( remove_channels,
                data["bkg_yields"] )
            data["obs_yields"]=self._removeSignalRegions ( remove_channels,
                data["obs_yields"] )
        self.onnxMeta = data

    def _predictFromScaledYields ( self, scaled_yields : np.array ) -> np.array:
        """ get the prediction from the NN

        :param scaled_yields: the input of the neural network
        :returns: arr, the unscaled unshifted output of the neural network
        """
        if len(scaled_yields[0])!=self.regressor["dim"]:
            dim_nn = self.regressor["dim"]
            dim_input = len(scaled_yields[0])
            line=f"the network of {self.onnxfilename} wants {dim_nn} input dimensions, but we supply {dim_input}. fix it!"
            print ( f"[nnAdapter] {line}" )
            print ( f"[nnAdapter] srOrder: {self.srOrder}" )
            sys.exit()
        dct = { "features": scaled_yields }
        arr = self.regressor["session"].run(None, dct )
        arr = arr[0][0]
        return arr

    def _postprocess( self, arr : np.ndarray,
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
        from smodels.statistics.nnPreprocessing import postprocess_nLLs
        deltas_prepd = np.array(arr, dtype=np.float64)
        trafos = self.onnxMeta["run_config"]["data"]["nLL_trafos"]
        nll_means = self.onnxMeta["nllMeans"]
        nll_errors = self.onnxMeta["nllErrors"]
        deltas = postprocess_nLLs ( deltas_prepd[:4],
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

        ret = { "nll_exp_0": nll0exp,  "nll_exp_1": nll1exp,
                "nll_obs_0": nll0obs,  "nll_obs_1": nll1obs,
                "nllA_exp_0": nllA0exp, "nllA_exp_1": nllA1exp,
                "nllA_obs_0": nllA0obs, "nllA_obs_1": nllA1obs }
        if self.onnxMeta["nLL_obs_max"][1] is not None:
            ret["nll_obs_max"] = self.onnxMeta["nLL_obs_max"][1]
        if add_errors:
            from smodels.statistics.nnPreprocessing import postprocess_nLLs_errors
            errs = postprocess_nLLs_errors ( deltas_prepd[4:],
                    deltas_prepd[:4],
                    mean = nll_means,
                    std = nll_errors, trafos = trafos,
                    eps = 1e-5 )
            errs = list ( map ( float, errs ) )
            ret["sigma_exp"] = errs[0]
            ret["sigma_obs"] = errs[1]
            ret["sigma_expA"] = errs[2]
            ret["sigma_obsA"] = errs[3]
        return ret

    def _totalYieldsFromSignals ( self, signal_yields : dict,
           obs_as_bg : list = [] ) -> dict:
        """ given the signal yields, return the total
        yields, signal + background

        :param signal_yields: the signal yields, as a (srname, yield) dictionary
        :param obs_abs_bg: a list of signal regions for which we use 
        observations as background_yields ("postfit")

        :returns: the total yields, as a dictionary
        """
        new_yields = {}
        account_for_crs = obs_as_bg[:]

        for srname,smyield in self.onnxMeta["bkg_yields"].items():
            assert srname in signal_yields, \
                f"nnInterface: cannot find sr name {srname} in '{ signal_yields }'"
            signal = signal_yields[srname]
            if srname in obs_as_bg:
                account_for_crs.remove ( srname )
                smyield = self.onnxMeta["obs_yields"][srname]
                signal = 0.
            tot = smyield + signal
            new_yields[srname] = tot
        if len(account_for_crs)>0:
            raise Exception ( f"signal region(s) {account_for_crs} unknown" )

        return new_yields


    def _preprocess ( self, yields : Union[dict,list] ) -> dict:
        if type(yields)==dict:
            yields = self._inputDictToList ( yields )
        inp_list = np.array ( yields )
        from smodels.statistics.nnPreprocessing import preprocess_features
        trafos = self.onnxMeta["run_config"]["data"]["trafos"]
        nYields = len(yields)
        re = preprocess_features ( inp_list,
            trafos = trafos,
            mean = np.array ( self.onnxMeta["featureMeans"][:nYields] ),
            std = np.array ( self.onnxMeta["featureErrors"][:nYields] ) )
        ret = [ re[0] ]
        return ret

    def _inputDictToList ( self, in_dict : dict ) -> list:
        """ translate a dictionary of input yields to a list
        of said yields in the canonical order specified in the onnx
        :raises: exception if input SR is missing or one too many
        :returns: list of yields
        """
        ret = []
        account_for_srs = list(in_dict.keys())
        if len(in_dict) != len ( self.srOrder ):
            raise Exception ( f"length of dict ({len(in_dict)} does not match srOrder ({len(self.srOrder)})" )
        for sr in self.srOrder:
            dsr = sr
            #if dsr.endswith ( "-0" ):
            #    dsr = sr[:-2]
            if sr in in_dict:
                ret.append ( in_dict[sr] )
                account_for_srs.remove ( sr )
                continue
            #if dsr in in_dict:
            #    ret.append ( in_dict[dsr] )
            #    continue
            raise Exception ( f"signal region {sr} not in input_dict" )
        if len(account_for_srs)>0:
            raise Exception ( f"signal region(s) {account_for_srs} unknown" )
        return ret
