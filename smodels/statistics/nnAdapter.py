#!/usr/bin/env python3

"""
.. module:: nnAdapter
   :synopsis: An Adapter class that wraps around the neural networks (e.g. onnx
   files), handle all the pre and post processing. This adapter is
   meant to be published with the ML paper, and the SModels database, not with
   SModelS per se.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os, onnx

class NNAdapter:
    """
    Adapter that wraps around a neural network
    """

    def __init__( self, onnxfile : os.PathLike, regions : list[dict],
                  allowsSyntheticData : bool = False ):
        """
        :param onnxfile: path to onnx file
        :param regions: a list of dictionaries of how SModelS SRs translate
        to the ML model's SRs
        :param allowsSyntheticData: if true, then also synthetic
        data can be supplied
        """
        self.onnxfile = onnxfile 
        self.regions = regions 
        self.mlModel = onnx.load ( onnxfile )
        self.allowsSyntheticData = allowsSyntheticData
        self.parseMetaData ()

    def parseMetaData ( self ):
        """ parse the model's meta data """
        data = { "inputMeans": [], "inputErrors": [], 
            "nLL_exp_mu0": [ None ]*2, "nLL_obs_mu0": [ None ]*2, 
            "nLLA_exp_mu0": [ None ]*2, "nLLA_obs_mu0": [ None ]*2 }
        # data [ "smYields" ] = {}
        # data [ "obsYields" ] = {}
        #data["nLL_obs_max"]= [ None ] * 2
        #data["nLL_exp_max"]= [ None ] * 2
        #data["nLLA_obs_max"]= [ None ] * 2
        #data["nLLA_exp_max"]= [ None ] * 2
        remove_channels=[]
        import json, math
        for em in self.mlModel.metadata_props:
            if em.key == "remove_channels":
                # remove these channels at the end, so that order does not matter
                remove_channels = eval(em.value)
                # data["remove_channels"] = remove_channels
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
                fillValues ( em.key, em.value )
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
        #if len(remove_channels)>0:
        #    data["smYields"]=removeSignalRegions ( remove_channels, data["smYields"] )
        #    data["obsYields"]=removeSignalRegions ( remove_channels, data["obsYields"] )
        # ['smYields', 'obsYields', 'inputMeans', 'inputErrors' ]
        self.onnxMeta={}
        for key,value in data.items():
            self.onnxMeta[key]=value

    def predict ( self, yields : dict ) -> dict:
        """ predict for yields 
        :param yields: e.g. { "SR1": 3, "SR2": 5 }

        :returns: xxx
        """
        ret = {}
        return ret

    def interact ( self ):
        import IPython; IPython.embed( colors = "neutral" )
        
def getRegionsForExample():
    # the regions dictionary of the example
    ret = [
        {'smodels': 'SRWZ_1', 'pyhf': 'SR1_WZ_cuts'},
        {'smodels': 'SRWZ_2', 'pyhf': 'SR2_WZ_cuts'},
        {'smodels': 'SRWZ_3', 'pyhf': 'SR3_WZ_cuts'},
        {'smodels': 'SRWZ_4', 'pyhf': 'SR4_WZ_cuts'},
        {'smodels': 'SRWZ_5', 'pyhf': 'SR5_WZ_cuts'},
        {'smodels': 'SRWZ_6', 'pyhf': 'SR6_WZ_cuts'},
        {'smodels': 'SRWZ_7', 'pyhf': 'SR7_WZ_cuts'},
        {'smodels': 'SRWZ_8', 'pyhf': 'SR8_WZ_cuts'},
        {'smodels': 'SRWZ_9', 'pyhf': 'SR9_WZ_cuts'},
        {'smodels': 'SRWZ_10', 'pyhf': 'SR10_WZ_cuts'},
        {'smodels': 'SRWZ_11', 'pyhf': 'SR11_WZ_cuts'},
        {'smodels': 'SRWZ_12', 'pyhf': 'SR12_WZ_cuts'},
        {'smodels': 'SRWZ_13', 'pyhf': 'SR13_WZ_cuts'},
        {'smodels': 'SRWZ_14', 'pyhf': 'SR14_WZ_cuts'},
        {'smodels': 'SRWZ_15', 'pyhf': 'SR15_WZ_cuts'},
        {'smodels': 'SRWZ_16', 'pyhf': 'SR16_WZ_cuts'},
        {'smodels': 'SRWZ_17', 'pyhf': 'SR17_WZ_cuts'},
        {'smodels': 'SRWZ_18', 'pyhf': 'SR18_WZ_cuts'},
        {'smodels': 'SRWZ_19', 'pyhf': 'SR19_WZ_cuts'},
        {'smodels': 'SRWZ_20', 'pyhf': 'SR20_WZ_cuts'}]
    return ret

if __name__ == "__main__":
    onnxfile = "../../unittests/test.onnx"
    regions = getRegionsForExample()
    adapter = NNAdapter ( onnxfile, regions, False )
    yields = {}
    for region in regions:
        yields[region["smodels"]] = 0.
    ret = adapter.predict ( yields )
    print ( ret )
    adapter.interact()
