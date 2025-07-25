#!/usr/bin/env python3

"""
.. module:: pyhfInterface
   :synopsis: Code that delegates the computation of limits and likelihoods to
              pyhf.

.. moduleauthor:: Gael Alguero <gaelalguero@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import jsonpatch
import warnings
import jsonschema
import copy
import numpy as np
from smodels.base.smodelsLogging import logger
from smodels.statistics.basicStats import findRoot
from smodels.tools.caching import roundCache, lru_cache
from smodels.matching.theoryPrediction import mu_digits
from smodels.statistics.basicStats import observed, apriori, aposteriori, NllEvalType
import logging
logging.getLogger("pyhf").setLevel(logging.CRITICAL)
# warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", r"invalid value encountered in log")
from typing import Dict, List, Union, Text

jsonver = ""
try:
    import importlib.metadata

    jsonver = int(importlib.metadata.version("jsonschema")[0])
except Exception as e:
    try:
        from jsonschema import __version__ as jsonver
    except Exception as e:
        pass
if jsonver < 3:
    # if jsonschema.__version__[0] == "2": ## deprecated
    print( f"[SModelS:pyhfInterface] jsonschema is version {jsonschema.__version__}, we need > 3.x.x" )
    sys.exit()

import time, sys, os

try:
    import pyhf
except ModuleNotFoundError:
    print("[SModelS:pyhfInterface] pyhf import failed. Is the module installed?")
    sys.exit(-1)

ver = pyhf.__version__

pyhfinfo = {
    "backend": "numpy",
    "hasgreeted": False,
    "backendver": np.version.full_version,
    "ver": ver,
#    "required": "0.6.1", # the required pyhf version
}

def setBackend ( backend : str ) -> bool:
    """
    try to setup backend to <backend>

    :param backend: one of: numpy (default), pytorch, jax, tensorflow
    :returns: True, if worked, False if failed
    """
    try:
        pyhf.set_backend( backend )
        pyhfinfo["backend"] = backend
        pyhfinfo["backendver"] = "?"
        module_name = backend
        if backend == "pytorch":
            module_name = "torch"
        from importlib.metadata import version
        pyhfinfo["backendver"] = version(module_name)
        return True
    except (pyhf.exceptions.ImportBackendError,pyhf.exceptions.InvalidBackend) as e:
        print( f"[SModelS:pyhfInterface] WARNING could not set {backend} as the pyhf backend: {e}" )
        print( f"[SModelS:pyhfInterface] falling back to {pyhfinfo['backend']}." )
        # print("[SModelS:pyhfInterface] We however recommend that pytorch be installed.")
    return False

# setBackend ( "pytorch" )

countWarning = {"llhdszero": 0}
# Sets the maximum number of attempts for determining the brent bracketing interval for mu:
nattempts_max = 10

def guessPyhfType ( name : str ) -> str:
    """ given the pyhf analysis region name,
    guess the type. awkward, phase this out! """
    if name.startswith ( "CR" ):
        return "CR"
    if name.startswith ( "VR" ):
        return "VR"
    if name.startswith ( "SR" ):
        return "SR"
    if "CR" in name: ## arggh!!
        logger.debug ( f"old jsonFiles format used, and 'CR' in the middle of the region name: {name}. I will assume it is a control region but do switch to the new format ASAP" )
        return "CR"
    return "SR"

class PyhfData:
    """
    Holds data for use in pyhf
    :ivar nsignals: signal predictions dictionary of dictionaries,
    one for each json file, one entry per signal region
    :ivar inputJsons: list of json instances
    :ivar jsonFiles: optional list of json files
    :ivar nWS: number of workspaces = number of json files
    """

    def __init__( self, nsignals : Dict[str, Dict], inputJsons, jsonFiles=None,
                  includeCRs=False, signalUncertainty=None):
        self.nsignals = nsignals
        self.getTotalYield()
        self.inputJsons = inputJsons
        if jsonFiles is None:   # If no name has been provided for the json file(s) and the channels, use fake ones
            jsonFiles = {}
            for jFile,sregions in nsignals.items():
                regions = []
                srname = "SR1"
                if len(sregions)==1:
                    regions.append ( { "smodels": srname, "type": "SR", "pyhf": srname } )
                else:
                    for j in range( len ( sregions ) ):
                        srname = f"SR1[{j}]"
                        regions.append ( { "smodels": srname, "type": "SR",
                                           "pyhf": srname } )
                jsonFiles[ jFile ] = regions
        self.jsonFiles = jsonFiles
        self.includeCRs = includeCRs
        self.signalUncertainty = signalUncertainty
        self.combinations = None
        if jsonFiles != None:
            self.combinations = [os.path.splitext(os.path.basename(js))[0] for js in jsonFiles]

        self.nWS = len(inputJsons)
        self.errorFlag = False
        self.checkConsistency()
        self.getWSInfo()

    def getTotalYield ( self ):
        """ the total yield in all signal regions """
        S = 0
        for dct in self.nsignals.values():
            for signal in dct.values():
                if isinstance(signal, list):
                    for sig in signal:
                        S += sig
                else:
                    S += signal
        self.totalYield = S

    def createPatchForRegion ( self, region, i_ch, ch, jsName ):
        chname = ch['name']
        chname2 = f'{ch["name"]}[0]' ## binned SRs
        if not region["pyhf"] in [ chname, chname2 ]:
            return None, None
        if (region['type'] == 'SR') or (region['type'] == 'CR' and self.includeCRs and region['smodels'] is not None):
            if region['smodels'] not in self.nsignals[jsName]:
                logger.error(f"Region {region['smodels']} of {jsName} not in the signal dictionary!")
                self.errorFlag = True
                return None, None
            nBins = len(ch["data"])
            # Find all smodels names if many share the same pyhf name (for multi-bin regions)
            smodelsName = []
            ctr = 0
            if region["pyhf"] == chname: # no bins, so easy
                smodelsName.append(region['smodels'])
            else: ## bins
                hasAdded = True
                while hasAdded:
                    hasAdded = False
                    for rregion in self.jsonFiles[jsName]:
                        if rregion["pyhf"] == f'{ch["name"]}[{ctr}]':
                            smodelsName.append(rregion['smodels'])
                            ctr+=1
                            hasAdded = True
            if len(smodelsName) != nBins:
                logger.error(f"Json region {region['pyhf']} has {nBins} bins, but only {len(smodelsName)} are implemented!")
                self.errorFlag = True
                return None, None

            signal = []
            for name in smodelsName: # In case of multiple signals for one region, the dict ordering within globalInfo.jsonFiles matters
                signal.append(self.nsignals[jsName][name])

            smodelsName = ";".join( smodelsName ) # Name of the corresponding region(s). Join all names if multiple bins.

            ret = {
                    "path": f"/channels/{i_ch}/samples/0",
                    "size": nBins,
                    "smodelsName": smodelsName,
                    "signal": signal
                }, "signalRegions"
            return ret
        ret = { 'path': f"/channels/{i_ch}", 'name': chname,
                'type': region['type']}, "otherRegions"
        return ( ret )

    def updatePyhfNames ( self, jsName : str, observations : List ):
        """ if no pyhf names are given, get them from the ws,
        in the order of the ws """
        if "pyhf" in self.jsonFiles[jsName][0]:
            return

        ## we dont have the mapping smodels<->pyhf
        ctr = 0
        #ic ( "---" )
        nJsonFiles = len(self.jsonFiles[jsName])
        #nObs = len(observations)
        #ic ( self.includeCRs, nObs, nJsonFiles )
        #ic ( observations )
        #ic ( self.jsonFiles[jsName] )
        nSRs = 0
        for observation in observations:
            name = observation["name"]
            regionType = guessPyhfType ( name )
            if regionType == "SR":
                nSRs += 1
        # ic ( nSRs, nCRs )
        for observation in observations:
            name = observation["name"]
            regionType = guessPyhfType ( name )
            if regionType in [ "VR" ]:
                region = { "pyhf": observation["name"], "smodels": None,
                           "type": regionType }
                self.jsonFiles[jsName].append ( region )
                continue
            if not self.includeCRs and regionType in [ "CR" ]:
                region = { "pyhf": observation["name"], "smodels": None,
                           "type": regionType }
                self.jsonFiles[jsName].append ( region )
                continue
            if self.includeCRs and regionType in [ "CR" ]:
                if nSRs == nJsonFiles: # and nSRs+nCRs == nObs:
                    ## the signal regions alone do it
                    region = { "pyhf": observation["name"], "smodels": None,
                               "type": regionType }
                    self.jsonFiles[jsName].append ( region )
                    continue

            if len(observation["data"])==1:
                if ctr < len(self.jsonFiles[jsName]):
                   self.jsonFiles[jsName][ctr]["pyhf"]=f"{name}"
                   self.jsonFiles[jsName][ctr]["type"]=regionType
                ctr += 1
            else:
                for i in range(len(observation["data"])):
                    if ctr < len(self.jsonFiles[jsName]):
                        self.jsonFiles[jsName][ctr]["pyhf"]=f"{name}[{i}]"
                        self.jsonFiles[jsName][ctr]["type"]=regionType
                    ctr += 1


    def getWSInfo(self):
        """
        Getting informations from the json files

        :ivar channelsInfo: list of dictionaries (one dictionary for each json file) containing useful information about the json files
            - :key signalRegions: list of dictonaries with 'json path' and 'size' (number of bins) of the 'signal regions' channels in the json files
            - :key otherRegions: list of dictionnaries indicating the path and the name of the control and/or validation region channels
        """
        if self.errorFlag:
            return

        # Identifying the path to the channels in the main workspace files
        self.channelsInfo = []  # workspace specifications
        if not isinstance(self.inputJsons, list):
            logger.error("The 'inputJsons' parameter must be of type list")
            self.errorFlag = True
            return

        for ws, jsName in zip(self.inputJsons, self.jsonFiles):
            wsChannelsInfo = {}
            wsChannelsInfo["signalRegions"] = []
            wsChannelsInfo["otherRegions"] = []
            if not "channels" in ws.keys():
                logger.error(
                    "Json file number {} is corrupted (channels are missing)".format(
                        self.inputJsons.index(ws)
                    )
                )
                self.channelsInfo = None
                return

            sigInCRs = False
            signalNames = self.nsignals[jsName].keys()
            for signalName in signalNames:
                for region in self.jsonFiles[jsName]:
                    if signalName == region['smodels'] and region['type'] == 'CR':
                        sigInCRs = True
            if sigInCRs and not self.includeCRs:
                logger.warning("Signal in CRs but includeCRs = False. CRs will still be removed.")

            # smodelsRegions = self.nsignals[jsName].values() # CR and SR names implemented in the database
            smodelsRegions = self.jsonFiles[jsName]
            #import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()
            if "observations" in ws:
                self.updatePyhfNames ( jsName, ws["observations"] )
                patchedChannels = set()
                allChannels = set ( [ x["name"] for x in ws["observations"] ] )
                for i_r, region in enumerate ( self.jsonFiles[jsName] ):
                    for i_ch, ch in enumerate(ws["observations"]):
                        ## create a patch for the region, but only if channel matches
                        patch, patchType = self.createPatchForRegion ( region, i_ch, ch, jsName )
                        if patch != None:
                            wsChannelsInfo[patchType].append(patch)
                            patchedChannels.add ( ch['name'] )
                            break
                if allChannels != patchedChannels:
                    logger.debug ( f"could not patch {' '.join(allChannels-patchedChannels)} for {jsName}. Check the database!" )
                    # sys.exit()


            wsChannelsInfo["otherRegions"].sort(
                key=lambda path: int(path['path'].split("/")[-1]), reverse=True
            )  # Need to sort correctly the paths to the channels to be removed
            self.channelsInfo.append(wsChannelsInfo)
            # ic ( wsChannelsInfo )

    def checkConsistency(self):
        """
        Check various inconsistencies of the PyhfData attributes

        :param zeroSignalsFlag: boolean identifying if all SRs of a single json are empty
        """
        if not isinstance(self.nsignals, dict):
            logger.error("The 'nsignals' parameter must be of type list")
            self.errorFlag = True

        if self.nWS != len(self.nsignals):
            logger.error(
                "The number of subsignals provided is different from the number of json files"
            )
            self.errorFlag = True
        self.zeroSignalsFlag = list()

        for jsName, subSig in zip(self.jsonFiles, self.nsignals.values()):
            if not isinstance(subSig, dict):
                logger.error("The 'nsignals' parameter must be a dictionary of dictionary")
                self.errorFlag = True
                return
            nBinsJson = 0
            for region in self.jsonFiles[jsName]:
                if (region['type'] == 'SR') or (region['type'] == 'CR' and self.includeCRs and region['smodels'] is not None):
                    nBinsJson += 1
            if nBinsJson != len(subSig):
                logger.error(
                    f"The number of signals ({len(subSig)}) provided is different from the number of signal bins for json ({nBinsJson}) for {jsName}" )
                self.errorFlag = True
            allZero = all([s == 0 for s in subSig.values()])
            # Checking if all signals matching this json are zero
            self.zeroSignalsFlag.append(allZero)


class PyhfUpperLimitComputer:
    """
    Class that computes the upper limit using the jsons files and signal informations in the 'data' instance of 'PyhfData'
    """

    def __init__(self, data, cl=0.95, lumi=None ):
        """

        :param data: instance of 'PyhfData' holding the signals information
        :param cl: confdence level at which the upper limit is desired to be computed
        :ivar data: created from data
        :ivar nsignals: signal predictions list divided into sublists, one for each json file
        :ivar inputJsons: list of input json files as python json instances
        :ivar channelsInfo: list of channels information for the json files
        :ivar zeroSignalsFlag: list boolean flags in case all signals are zero for a specific json
        :ivar nWS: number of workspaces = number of json files
        :ivar patches: list of patches to be applied to the inputJsons as python dictionary instances
        :ivar workspaces: list of workspaces resulting from the patched inputJsons
        :ivar workspaces_expected: list of patched workspaces with observation yields replaced by the evaluationType ones
        :ivar cl: created from cl
        :ivar scale: scale that is applied to the signal predictions, dynamically changes throughout the upper limit calculation
        :ivar alreadyBeenThere: boolean flag that identifies when nsignals accidentally passes twice at two identical values
        """

        self.data = data
        self.lumi = lumi
        self.nsignals = copy.deepcopy ( self.data.nsignals )
        logger.debug(f"Signals : {self.nsignals}")
        self.inputJsons = self.data.inputJsons
        self.channelsInfo = None
        if hasattr ( self.data, "channelsInfo" ):
            self.channelsInfo = self.data.channelsInfo
        self.zeroSignalsFlag = self.data.zeroSignalsFlag
        self.nWS = self.data.nWS
        self.includeCRs = data.includeCRs
        self.patches = self.patchMaker()
        self.workspaces = self.wsMaker()
        self.workspaces_expected = self.wsMaker(apriori=True)
        self.cl = cl
        self.scale = 1.0
        self.sigma_mu = None
        self.alreadyBeenThere = (
            False  # boolean to detect wether self.signals has returned to an older value
        )
        # self.checkPyhfVersion()
        self.welcome()

    def welcome(self):
        """
        greet the world
        """

        if pyhfinfo["hasgreeted"]:
            return
        logger.info(
            f"Pyhf interface, we are using v{str(pyhfinfo['ver'])}, with {pyhfinfo['backend']} v{str(pyhfinfo['backendver'])} as backend."
        )
        pyhfinfo["hasgreeted"] = True

    def checkPyhfVersion(self):
        """
        check the pyhf version, currently we need 0.6.1+
        """

        if pyhfinfo["ver"] < pyhfinfo["required"]:
            logger.warning(
                f"pyhf version is {str(pyhfinfo['ver'])}. SModelS currently requires pyhf>={str(pyhfinfo['required'])}. You have been warned."
            )

    def rescale(self, factor):
        """
        Rescales the signal predictions (self.nsignals) and processes again the patches and workspaces

        :return: updated list of patches and workspaces (self.patches, self.workspaces and self.workspaces_expected)
        """

        for jsName in self.nsignals.keys():
            for regionName in self.nsignals[jsName].keys():
                self.nsignals[jsName][regionName] = self.nsignals[jsName][regionName]*factor
        try:
            self.alreadyBeenThere = self.nsignals == self.nsignals_2
        except AttributeError:
            pass
        self.scale *= factor
        logger.debug(f"new signal scale : {self.scale}")
        self.patches = self.patchMaker()
        self.workspaces = self.wsMaker()
        self.workspaces_expected = self.wsMaker(apriori=True)
        try:
            self.nsignals_2 = self.nsignals_1.copy()  # nsignals at previous-to-previous loop
        except AttributeError:
            pass
        self.nsignals_1 = self.nsignals.copy()  # nsignals at previous loop

    def changeChannelName ( self, srInfo ):
        """ changes the channel names in the json to match the SModelS name.
        FIXME this is only a hack for now, should be replaced by a
        self-respecting dataIdMap in the database itself,
        that maps the smodels names to the pyhf names explicitly.
        This method will then go away!
        """
        operators= []
        operator = {} # Operator for renaming the channels according to their region name from the database
        operator["op"] = "replace"
        operator["path"] = srInfo["path"].replace('samples/0','name')
        operator["value"] = srInfo["smodelsName"]
        operators.append(operator)

        operator = {} # Operator for renaming the observations according to their region name from the database
        operator["op"] = "replace"
        operator["path"] = srInfo["path"].replace('channels','observations').replace('samples/0','name')
        operator["value"] = srInfo["smodelsName"]
        operators.append(operator)
        return operators


    def patchMaker(self):
        """
        Method that creates the list of patches to be applied to the self.inputJsons workspaces, one for each region given the self.nsignals and the informations available in self.channelsInfo and the content of the self.inputJsons
        NB: It seems we need to include the change of the "modifiers" in the patches as well

        :return: the list of patches, one for each workspace
        """

        if self.channelsInfo == None:
            return None
        # Constructing the patches to be applied on the main workspace files
        patches = []
        for ws, info, (jsFileName,jsFile) in zip(self.inputJsons, self.channelsInfo, self.data.jsonFiles.items() ):
            patch = []
            for srInfo in info["signalRegions"]:
                operator = {} # Operator for patching the signal
                operator["op"] = "add"
                operator["path"] = srInfo["path"]
                value = {}
                #value["data"] = srInfo['signal']
                sr_order = srInfo["smodelsName"].split(";")
                nsignals = self.nsignals[jsFileName]
                # ic ( nsignals, sr_order, srInfo, jsFileName )
                value["data"] = [ nsignals[x] for x in sr_order ]
                #import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()
                # sys.exit()
                value["modifiers"] = []
                value["modifiers"].append({"data": None, "type": "normfactor", "name": "mu_SIG"})
                value["modifiers"].append({"data": None, "type": "lumi", "name": "lumi"})
                if self.data.signalUncertainty is not None and self.data.signalUncertainty != 0:
                    # Uncomment the line below to add a MC statistical uncertainty.
                    # value["modifiers"].append({"data": [sigBin*self.data.signalUncertainty for sigBin in value["data"]], "type": "staterror", "name": "MCError"})
                    value["modifiers"].append({"data": {"hi_data": [sigBin*(1+self.data.signalUncertainty) for sigBin in value["data"]],
                                                        "lo_data": [sigBin*(1-self.data.signalUncertainty) for sigBin in value["data"]]
                                                       },
                                               "type": "histosys",
                                               "name": "signalUncertainty"
                                              })
                value["name"] = "bsm"
                operator["value"] = value
                patch.append(operator)

                ## FIXME this if block is a hack, only to be used until
                ## smodels-database has proper dataIdMaps, mapping the smodels
                ## SR names to the pyhf ones. once these dataIdMaps are in place,
                ## they should be used instead of this hack that rewrites
                ## the pyhf channel names
                if False: # srInfo["smodelsName"]: # If the CRs/SRs have a name in the database (it is always True when running SModelS the usual way)
                    operators = self.changeChannelName ( srInfo )
                    for operator in operators:
                        patch.append(operator)

            for region in info["otherRegions"]:
                if region['type'] == 'CR' and self.includeCRs:
                    continue
                else:
                    patch.append({"op": "remove", "path": region['path']}) # operator for removing useless regions
            patches.append(patch)

        return patches

    def wsMaker(self, apriori=False):
        """
        Apply each region patch (self.patches) to his associated json (self.inputJsons) to obtain the complete workspaces

        :param apriori: - If set to 'True': Replace the observation data entries of each workspace by the corresponding sum of the evaluationType yields \
                        - Else: The observed yields put in the workspace are the ones written in the corresponfing json dictionary

        :returns: the list of patched workspaces
        """
        if self.patches == None:
            return None
        if self.nWS == 1:
            try:
                wsDict = jsonpatch.apply_patch(self.inputJsons[0], self.patches[0])
                if apriori == True:
                    # Replace the observation data entries by the corresponding sum of the evaluationType yields
                    for obs in wsDict["observations"]:
                        for ch in wsDict["channels"]:
                            # Finding matching observation and bkg channel
                            if obs["name"] == ch["name"]:
                                bkg = [0.0] * len(obs["data"])
                                for sp in ch["samples"]:
                                    if sp["name"] == "bsm":
                                        continue
                                    for iSR in range(len(obs["data"])):
                                        # Summing over all bkg samples for each bin/SR
                                        bkg[iSR] += sp["data"][iSR]
                                # logger.debug('bkgs for channel {} :\n{}'.format(obs['name'], bkg))
                                obs["data"] = bkg
                return [pyhf.Workspace(wsDict)]
            except (pyhf.exceptions.InvalidSpecification, KeyError) as e:
                logger.error(f"The json file is corrupted:\n{e}")
                return None
        else:
            workspaces = []
            for js, patch in zip(self.inputJsons, self.patches):
                wsDict = jsonpatch.apply_patch(js, patch)
                if apriori == True:
                    # Replace the observation data entries by the corresponding sum of the evaluationType yields
                    for obs in wsDict["observations"]:
                        for ch in wsDict["channels"]:
                            # Finding matching observation and bkg channel
                            if obs["name"] == ch["name"]:
                                bkg = [0.0] * len(obs["data"])
                                for sp in ch["samples"]:
                                    if sp["name"] == "bsm":
                                        continue
                                    for iSR in range(len(obs["data"])):
                                        # Summing over all bkg samples for each bin/SR
                                        bkg[iSR] += sp["data"][iSR]
                                # logger.debug('bkgs for channel {} :\n{}'.format(obs['name'], bkg))
                                obs["data"] = bkg
                try:
                    ws = pyhf.Workspace(wsDict)
                except (pyhf.exceptions.InvalidSpecification, KeyError) as e:
                    logger.error(
                        "Json file number {} is corrupted:\n{}".format(
                            self.inputJsons.index(json), e
                        )
                    )
                    return None
                workspaces.append(ws)

            return workspaces

    def backup(self):
        self.bu_signal = copy.deepcopy(self.data.nsignals)

    def restore(self):
        if not hasattr(self, "bu_signal"):
            return
        self.data.nsignals = copy.deepcopy(self.bu_signal)
        del self.bu_signal

    def get_position(self, name, model):
        """
        :param name: name of the parameter one wants to increase
        :param model: the pyhf model
        :return: the position of the parameter that has to be modified in order to turn positive the negative total yield
        """
        position = 0
        for par in model.config.par_order:
            if name == par:
                return position
            else:
                position += model.config.param_set(par).n_parameters

    def rescaleBgYields(self, init_pars, workspace, model):
        """
        Increase initial value of nuisance parameters until the starting value of the total yield (mu*signal + background) is positive

        :param init_pars: list of initial parameters values one wants to increase in order to turn positive the negative total yields
        :param workspace: the pyhf workspace
        :param model: the pyhf model
        :return: the list of initial parameters values that gives positive total yields
        """
        for pos,yld in enumerate(model.expected_actualdata(init_pars)):
            if yld < 0:
                sum_bins = 0
                positive = False
                # Find the SR and the bin (the position inside the SR) corresponding to the negative total yield
                for channel,nbins in model.config.channel_nbins.items():
                    sum_bins += nbins
                    if pos < sum_bins:
                        SR = channel
                        bin_num = pos-sum_bins+nbins # The position of the bin in the SR
                        break
                # Find the name of the parameter that modifies the background yield
                for channel in workspace['channels']:
                    if channel['name'] == SR:
                        for sample in range(len(channel['samples'])):
                            for modifier in channel['samples'][sample]['modifiers']:
                                # The multiplicative modifier that will be increased is of type 'staterror' (others may work as well)
                                # In case of a simplified json, let's take the only parameter called 'totalError'
                                if ('type' in modifier.keys() and modifier['type'] == 'staterror') or ('type' in modifier.keys() and modifier['type'] == 'normsys') or modifier['name'] == 'totalError':
                                    name = modifier['name']
                                    # Find the position of the parameter within the pyhf list of parameters
                                    position = self.get_position(name,model)+bin_num
                                    # Find the upper bound of the parameter that will be increased
                                    max_bound = model.config.suggested_bounds()[position][1]
                                    # If the parameter one wants to increase is not fixed, add 0.1 to its value (it is an arbitrary step) until
                                    # the total yield becomes positive or the parameter upper bound is reached
                                    if not model.config.suggested_fixed()[position]:
                                        while model.expected_actualdata(init_pars)[pos] < 0:
                                            if init_pars[position] + 0.1 < max_bound:
                                                init_pars[position] += 0.1
                                            else:
                                                init_pars[position] = max_bound
                                                break
                                        positive = model.expected_actualdata(init_pars)[pos] >= 0
                                if positive:
                                    break
                            if positive:
                                break
        return init_pars

    @roundCache(argname='mu',argpos=1,digits=mu_digits)
    def likelihood( self, mu : float = 1.0,
            workspace_index : Union[None,int] =None,
            return_nll : bool = False, evaluationType : NllEvalType=observed,
            asimov : Union[None,float] = None ):
        """
        Returns the value of the likelihood. \
        Inspired by the 'pyhf.infer.mle' module but for non-log likelihood

        :param workspace_index: supply index of workspace to use. If None, \
                                choose index of best combo
        :param return_nll: if true, return nll, not llhd
        :param expected: one of: observed, apriori, aposteriori
        """
        logger.debug("Calling likelihood")
        if type(workspace_index) == float:
            logger.error("workspace index is float")
        # logger.error("expected flag needs to be heeded!!!")
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                "Values in x were outside bounds during a minimize step, clipping to bounds",
            )
            # warnings.filterwarnings ( "ignore", "", module="pyhf.exceptions" )
            old_index = workspace_index
            if workspace_index == None:
                workspace_index = self.getBestCombinationIndex()
            if workspace_index == None:
                return None
            self.backup()
            try:
                if abs(mu - 1.0) > 1e-6:
                    for jsName in self.data.nsignals.keys():
                        for regionName in self.data.nsignals[jsName].keys():
                            self.data.nsignals[jsName][regionName] = self.data.nsignals[jsName][regionName]*mu
                self.__init__(self.data, self.cl, self.lumi)
                model,data,workspace = self.generateAsimovData ( asimov, workspace_index =
                       workspace_index, evaluationType = evaluationType )

                _, nllh = pyhf.infer.mle.fixed_poi_fit(
                    1.0, data, model, return_fitted_val=True, maxiter=200
                )
            except (pyhf.exceptions.FailedMinimization, ValueError) as e:
                logger.error(f"pyhf fixed_poi_fit failed for {list(self.data.jsonFiles)[workspace_index]} for mu={mu}: {e}")
                # lets try with different initialisation
                init, n_ = pyhf.infer.mle.fixed_poi_fit(
                    0.0, data, model, return_fitted_val=True, maxiter=200
                )
                initpars = init.tolist()
                initpars[model.config.poi_index] = 1.
                # Try to turn positive all the negative total yields (mu*signal + background) evaluated with the initial parameters
                initpars = self.rescaleBgYields(initpars, workspace, model)
                # If the a total yield is still negative with the increased initial parameters, print a message
                if not all([True if yld >= 0 else False for yld in model.expected_actualdata(initpars)]):
                    logger.debug(f'Negative total yield after increasing the initial parameters for mu={mu}')
                try:
                    bestFitParam, nllh = pyhf.infer.mle.fixed_poi_fit(
                        1.0,
                        data,
                        model,
                        return_fitted_val=True,
                        init_pars=initpars,
                        maxiter=2000,
                    )
                    # If a total yield is negative with the best profiled parameters, return None
                    if not all([True if yld >= 0 else False for yld in model.expected_actualdata(bestFitParam)]):
                        self.restore()
                        return self.exponentiateNLL(None, not return_nll)
                except (pyhf.exceptions.FailedMinimization, ValueError) as e:
                    logger.info(f"pyhf fixed_poi_fit failed twice for {list(self.data.jsonFiles)[workspace_index]} for mu={mu}: {e}")

                    self.restore()
                    return self.exponentiateNLL(None, not return_nll)

            ret = nllh.tolist()
            try:
                ret = float(ret)
            except:
                ret = float(ret[0])
            ret = self.exponentiateNLL(ret, not return_nll)
            # print ( "now leaving the fit mu=", mu, "llhd", ret, "nsig was", self.data.nsignals )
            self.restore()
            # print ( f"@@PI likelihood {mu} {workspace_index} {expected} {asimov} {ret}" )
            return ret

    @lru_cache
    def getBestCombinationIndex(self):
        """find the index of the best evaluationType combination"""
        if self.nWS == 1:
            return 0
        logger.debug( f"Finding best evaluationType combination among {self.nWS} workspace(s)" )
        ulMin = float("+inf")
        i_best = None
        for i_ws in range(self.nWS):
            if self.data.totalYield == 0.:
                continue
            if self.zeroSignalsFlag[i_ws] == True:
                logger.debug( f"Workspace number {i_ws} has zero signals" )
                continue
            else:
                ul = self.getUpperLimitOnMu(evaluationType=apriori, workspace_index=i_ws)
            if ul == None:
                continue
            if ul < ulMin:
                ulMin = ul
                i_best = i_ws
        return i_best

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

    def compute_invhess(self, x, data, model, index, epsilon=1e-05):
        '''
        if inv_hess is not given by the optimiser, calculate numerically by evaluating second order partial derivatives using 2 point central finite differences method
        :param x: parameter values given to pyhf.infer.mle.twice_nll taken from pyhf.infer.mle.fit - optimizer.x (best_fit parameter values)
        :param data: workspace.data(model) passed to pyhf.infer.mle.fit
        :param model: model passed to pyhf.infer.mle.fit
        :param index: index of the POI
        Note : If len(x) <=5, compute the entire hessian matrix and ind its inverse. Else, compute the hessian at the index of the POI and return its inverse (diagonal approximation)
        returns the inverse hessian at the index of the poi
        '''

        n = len(x)

        if n<=5:

            hessian = np.zeros((n,n))
            for i in range(n):
                for j in range(i+1):
                    eps_i = epsilon*np.eye(n)[:,i] #identity along ith column
                    eps_j = epsilon*np.eye(n)[:,j]

                    #twice_nll is the objective function, hence need to find its hessian
                    par_11 = pyhf.infer.mle.twice_nll(x + eps_i + eps_j, data, model)
                    par_12 = pyhf.infer.mle.twice_nll(x - eps_i + eps_j, data, model)
                    par_21 = pyhf.infer.mle.twice_nll(x + eps_i - eps_j, data, model)
                    par_22 = pyhf.infer.mle.twice_nll(x - eps_i - eps_j, data, model)

                    partial_xi_xj = (par_11 - par_12 - par_21 +par_22)/(4*epsilon**2)
                    hessian[i,j] = partial_xi_xj
                    if i!=j: hessian[j,i] = partial_xi_xj

            def is_positive_definite(matrix):
                eigenvalues = np.linalg.eigvals(matrix)
                return all(eig > 0 for eig in eigenvalues)

            if not is_positive_definite(hessian):
                #raise ValueError("Hessian Matrix is not positive definite")
                logger.warning("Hessian Matrix is not positive definite")

            inverse_hessian = np.linalg.inv(hessian)

            #return the inverse hessian at the poi
            return inverse_hessian[index][index]

        #calculate only the hessian at the poi and return its inverse (an approximation!)
        eps_i = epsilon*np.eye(n)[:,index]
        par_11 = pyhf.infer.mle.twice_nll(x + 2*eps_i, data, model)
        par_12 = pyhf.infer.mle.twice_nll(x, data, model)
        par_22 = pyhf.infer.mle.twice_nll(x - 2*eps_i, data, model)
        hessian = (par_11 - 2*par_12 + par_22)/(4*epsilon**2)

        #return the inverse hessian at the poi
        return 1.0/hessian

    @lru_cache
    def lmax( self, workspace_index=None, return_nll=False, evaluationType : NllEvalType=observed,
              allowNegativeSignals=False):
        """
        Returns the negative log max likelihood

        :param return_nll: if true, return nll, not llhd
        :param workspace_index: supply index of workspace to use. If None, \
            choose index of best combo
        :param expected: one of: observed, apriori, aposteriori
        :param allowNegativeSignals: if False, then negative nsigs are replaced \
            with 0.
        """
        # logger.error("expected flag needs to be heeded!!!")
        logger.debug("Calling lmax")
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                "Values in x were outside bounds during a minimize step, clipping to bounds",
            )

            self.__init__(self.data, self.cl, self.lumi)
            old_index = workspace_index
            if workspace_index == None:
                workspace_index = self.getBestCombinationIndex()
            if workspace_index != None:
                if self.zeroSignalsFlag[workspace_index] == True:
                    logger.warning( f"Workspace number {workspace_index} has zero signals" )
                    return None
                else:
                    workspace = self.updateWorkspace(workspace_index, evaluationType=evaluationType)
            else:
                return None
            # Same modifiers_settings as those used when running the 'pyhf cls' command line
            msettings = {"normsys": {"interpcode": "code4"}, "histosys": {"interpcode": "code4p"}}
            model = workspace.model(modifier_settings=msettings)
            muhat, maxNllh = model.config.suggested_init(), float("nan")
            sigma_mu = float("nan")
            # obs = workspace.data(model)
            try:
                bounds = model.config.suggested_bounds()
                if allowNegativeSignals:
                    bounds[model.config.poi_index] = (-5., 10. )

                o = None
                try:
                    muhat, maxNllh, o = pyhf.infer.mle.fit(workspace.data(model), model,
                            return_fitted_val=True, par_bounds = bounds, return_result_obj = True )
                    #removed jacobain way of computing sigma_mu

                except (pyhf.exceptions.FailedMinimization,ValueError) as e:
                    pass

                if hasattr ( o, "hess_inv" ): # maybe the backend gets changed
                    sigma_mu = float ( np.sqrt ( o.hess_inv[model.config.poi_index][model.config.poi_index] ) ) * self.scale
                else:
                    sigma_mu_temp = 1.
                    try:
                        warnings.filterwarnings( "ignore", category=RuntimeWarning )
                        _1, _2, o = pyhf.infer.mle.fit(workspace.data(model), model,
                                return_fitted_val=True, return_result_obj = True, init_pars = list(muhat), method="BFGS" )
                        sigma_mu_temp = float ( np.sqrt ( o.hess_inv[model.config.poi_index][model.config.poi_index] ) )
                    except (pyhf.exceptions.FailedMinimization,ValueError) as e:
                        pass
                    if abs ( sigma_mu_temp - 1.0 ) > 1e-5:
                        sigma_mu = sigma_mu_temp * self.scale
                    else:
                        warnings.filterwarnings( "ignore", category=RuntimeWarning )
                        _, _, o = pyhf.infer.mle.fit(workspace.data(model), model,
                            return_fitted_val=True, return_result_obj = True, method="L-BFGS-B" )
                        sigma_mu_temp = float ( np.sqrt ( o.hess_inv.todense()[model.config.poi_index][model.config.poi_index] ) )
                        if abs ( sigma_mu_temp - 1.0 ) < 1e-8: # Fischer information is nonsense
                            #Calculate inv_hess numerically
                            inv_hess = self.compute_invhess(o.x, workspace.data(model), model, model.config.poi_index)
                            sigma_mu_temp = float ( np.sqrt ( inv_hess))
                        if abs ( sigma_mu_temp - 1.0 ) > 1e-8:
                            sigma_mu = sigma_mu_temp * self.scale
                        else: # sigma_mu is still nonsense
                            logger.warning("sigma_mu computation failed, even with inv_hess numercially computed. sigma_mu will be set to 1.")
                            sigma_mu = 1.

            except (pyhf.exceptions.FailedMinimization, ValueError) as e:
                if pyhfinfo["backend"] == "pytorch":
                    logger.warning(f"pyhf mle.fit failed {e} with {pyhfinfo['backend']} v{pyhfinfo['backendver']}. Calculating inv_hess numerically.")
                if pyhfinfo["backend"] == "numpy":
                    logger.debug(f"pyhf mle.fit failed {e} with {pyhfinfo['backend']} v{pyhfinfo['backendver']}. Calculating inv_hess numerically.")

                #Calculate inv_hess numerically
                inv_hess = self.compute_invhess(o.x, workspace.data(model), model, model.config.poi_index)
                sigma_mu = float ( np.sqrt ( inv_hess)) * self.scale

            muhat = float ( muhat[model.config.poi_index]*self.scale )
            try:
                lmax = maxNllh.tolist()
            except:
                lmax= maxNllh
            try:
                lmax = float(lmax)
            except:
                lmax = float(lmax[0])
            lmax = self.exponentiateNLL(lmax, not return_nll)

            ret = { "lmax": lmax, "muhat": muhat, "sigma_mu": sigma_mu }
            return ret

    def updateWorkspace(self, workspace_index=None, evaluationType : NllEvalType=observed ):
        """
        Small method used to return the appropriate workspace

        :param workspace_index: the index of the workspace to retrieve from the corresponding list
        :param expected: if False, retuns the unmodified (but patched) workspace. Used for computing observed or aposteriori evaluationType limits.
                        if True, retuns the modified (and patched) workspace, where obs = sum(bkg). Used for computing apriori evaluationType limit.
        """
        if self.nWS == 1:
            if evaluationType == apriori:
                return self.workspaces_expected[0]
            else:
                return self.workspaces[0]
        else:
            if workspace_index == None:
                logger.error("No workspace index was provided.")
            if evaluationType == apriori:
                return self.workspaces_expected[workspace_index]
            else:
                return self.workspaces[workspace_index]

    def getUpperLimitOnSigmaTimesEff(self, evaluationType : NllEvalType=observed, workspace_index=None):
        """
        Compute the upper limit on the fiducial cross section sigma times efficiency:
            - by default, the combination of the workspaces contained into self.workspaces
            - if workspace_index is specified, self.workspace[workspace_index]
              (useful for computation of the best upper limit)

        :param expected:  - if set to apriori: uses evaluationType SM backgrounds as signals
                          - else: uses 'self.nsignals'
        :param workspace_index: - if different from 'None': index of the workspace to use
                                  for upper limit
                                - else: choose best combo
        :return: the upper limit on sigma times eff at 'self.cl' level (0.95 by default)
        """
        if self.data.totalYield == 0.:
            return None
        else:
            ul = self.getUpperLimitOnMu( evaluationType=evaluationType, workspace_index=workspace_index)
            if ul == None:
                return ul
            if self.lumi is None:
                logger.error(f"asked for upper limit on fiducial xsec, but no lumi given with the data")
                return ul
            xsec = self.data.totalYield / self.lumi
            return ul * xsec

    def generateAsimovData ( self, mu : Union[None,float] = 0.,
            workspace_index : Union[int,None] = None,
            evaluationType : NllEvalType=observed ) -> list:
        """ generate asimov data for the model for signal strength mu
        :param mu: if None, then actually return the original data
        :returns: tuple with workspace and list with asimov data
        """
        msettings = {
            "normsys": {"interpcode": "code4"},
            "histosys": {"interpcode": "code4p"},
        }
        if workspace_index == None:
            workspace_index = self.getBestCombinationIndex()
        if workspace_index == None:
            return None
        workspace = self.updateWorkspace(workspace_index, evaluationType=evaluationType)
        #with warnings.catch_warnings():
        #    warnings.simplefilter("ignore", category=(DeprecationWarning,UserWarning))
        model = workspace.model(modifier_settings=msettings)
        data = workspace.data(model)
        if mu == None:
            return (model,data,workspace)
        ad = pyhf.infer.calculators.generate_asimov_data(mu, data, model, None, None, None)
        return (model,ad,workspace)

    @roundCache(argname='mu',argpos=1,digits=mu_digits)
    def CLs( self, mu : float, evaluationType : NllEvalType,
             return_type: Text = "CLs",
             workspace_index : Union[int,None] = None ) -> float:
        """
        This is the CLs method that heeds self.scale, i.e. you give it the
        _absolute_ mu

        :param mu: compute for the parameter of interest mu
        :param expected: if observed, compute observed, apriori: compute a priori expected
        :param return_type: (Text) can be one of:
        "CLs-alpha", "1-CLs", "CLs" "alpha-CLs"
        CLs-alpha: returns CLs - 0.05
        alpha-CLs: return 0.05 - CLs
        1-CLs: returns 1-CLs value
        CLs: returns CLs value
        """
        return self._CLs ( mu / self.scale, evaluationType, return_type, workspace_index )

    def _CLs( self, mu_rel : float, evaluationType : NllEvalType,
             return_type: Text = "CLs",
             workspace_index : Union[int,None] = None ) -> float:
        """
        This is our internal method to compute CLs.

        :param mu_rel: compute for the parameter of interest mu_rel
        :param expected: one of: observed, apriori, aposteriori
        :param return_type: (Text) can be one of:
        "CLs-alpha", "1-CLs", "CLs" "alpha-CLs"
        CLs-alpha: returns CLs - 0.05
        alpha-CLs: return 0.05 - CLs
        1-CLs: returns 1-CLs value
        CLs: returns CLs value
        """
        assert return_type in ["CLs-alpha", "alpha-CLs", "1-CLs", "CLs"], f"Unknown return type: {return_type}."
        # If evaluationType == False, use unmodified (but patched) workspace
        # If evaluationType == True, use modified workspace where observations = sum(bkg) (and patched)
        # If evaluationType == posteriori, use unmodified (but patched) workspace
        if workspace_index == None:
            workspace_index = self.getBestCombinationIndex()
        if workspace_index == None:
            return None
        workspace = self.updateWorkspace(workspace_index, evaluationType=evaluationType)
        # Same modifiers_settings as those use when running the 'pyhf cls' command line
        msettings = {
            "normsys": {"interpcode": "code4"},
            "histosys": {"interpcode": "code4p"},
        }
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=(DeprecationWarning,UserWarning))
            model = workspace.model(modifier_settings=msettings)
            bounds = model.config.suggested_bounds()
            bounds[model.config.poi_index] = (0, 10)
            start = time.time()
            args = {}
            args["return_expected"] = evaluationType == aposteriori
            args["par_bounds"] = bounds
            # args["maxiter"]=100000
            pver = float(pyhf.__version__[:3])
            stat = "qtilde"
            if pver < 0.6:
                args["qtilde"] = True
            else:
                args["test_stat"] = stat
            with np.testing.suppress_warnings() as sup:
                if pyhfinfo["backend"] == "numpy":
                    sup.filter(RuntimeWarning, r"invalid value encountered in log")
                # print ("expected", evaluationType, "return_expected", args["return_expected"], "mu", mu, "\nworkspace.data(model) :", workspace.data(model, include_auxdata = False), "\nworkspace.observations :", workspace.observations, "\nobs[data] :", workspace['observations'])
                # ic ( workspace["channels"][0]["samples"][0]["data"] )
                # import sys, IPython; IPython.embed( colors = "neutral" ); sys.exit()
                try:
                    result = pyhf.infer.hypotest(mu_rel, workspace.data(model), model, **args)
                except Exception as e:
                    logger.error(f"when testing hypothesis {mu_rel}, caught exception: {e}")
                    result = float("nan")
                    if evaluationType == aposteriori:
                        result = [float("nan")] * 2
            end = time.time()
            logger.debug( f"Hypotest elapsed time : {end-start:1.4f} secs" )
            logger.debug(f"result for {mu_rel} {result}")
            if evaluationType == aposteriori:
                logger.debug("computing a-posteriori evaluationType limit")
                logger.debug(f"expected = {expected}, mu_rel = {mu_rel}, result = {result}")
                try:
                    CLs = float(result[1].tolist())
                except TypeError:
                    CLs = float(result[1][0])
            else:
                logger.debug(f"expected = {evaluationType}, mu_rel = {mu_rel}, result = {result}")
                CLs = float(result)
            # print ( f"@@PI0 CLs for {mu_rel}*{self.scale} {expected} {workspace_index} {self.nsignals} = {CLs}" )
            if return_type == "1-CLs":
                return 1.0 - CLs
            elif return_type == "CLs":
                return CLs
            elif return_type == "CLs-alpha":
                return CLs - 1 + self.cl
            # return_type == "alpha-CLs"
            return 1 - self.cl - CLs

    # Trying a new method for upper limit computation :
    # re-scaling the signal predictions so that mu falls in [0, 10] instead of
    # looking for mu bounds
    # Usage of the index allows for rescaling
    @lru_cache
    def getUpperLimitOnMu(self, evaluationType : NllEvalType=observed, workspace_index=None):
        """
        Compute the upper limit on the signal strength modifier with:
            - by default, the combination of the workspaces contained into self.workspaces
            - if workspace_index is specified, self.workspace[workspace_index]
              (useful for computation of the best upper limit)

        :param expected:  - if set to apriori: uses evaluationType SM backgrounds as signals
                          - else: uses 'self.nsignals'
        :param workspace_index: - if different from 'None': index of the workspace to use
                                  for upper limit
                                - else: choose best combo
        :return: the upper limit at 'self.cl' level (0.95 by default)
        """
        self.__init__(self.data, self.cl, self.lumi)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                "Values in x were outside bounds during a minimize step, clipping to bounds",
            )
            startUL = time.time()
            logger.debug("Calling getUpperLimitOnMu")
            if self.data.errorFlag or self.workspaces == None:
                # For now, this flag can only be turned on by PyhfData.checkConsistency
                return None

            if (
                all([self.zeroSignalsFlag[workspace_index] for workspace_index in range(self.nWS)])
                == True
            ):
                logger.debug(
                    f"There is (are) {self.nWS} workspace(s) and no signal(s) was (were) found" )
                return None
            if workspace_index == None:
                workspace_index = self.getBestCombinationIndex()
            if workspace_index == None:
                logger.debug("Best combination index not found")
                return None

            # Rescaling signals so that mu is in [0, 10]
            factor = 3.0
            wereBothLarge = False
            wereBothTiny = False
            nattempts = 0
            nNan = 0
            lo_mu, med_mu, hi_mu = 0.2, 1.0, 5.0
            # print ( "starting with expected", evaluationType )
            while "mu is not in [lo_mu,hi_mu]":
                nattempts += 1
                if nNan > 5:
                    # logger.warning("encountered NaN 5 times while trying to determine the bounds for brent bracketing. now trying with q instead of qtilde test statistic")
                    return None
                    # nattempts = 0
                if nattempts > nattempts_max:
                    logger.warning(
                        f"tried {nattempts_max} times to determine the bounds for brent bracketing. we abort now."
                    )
                    return None
                # Computing CL(1) - 0.95 and CL(10) - 0.95 once and for all
                rt1 = self._CLs(lo_mu, evaluationType, "alpha-CLs", workspace_index )
                # rt5 = CLs(med_mu)
                rt10 = self._CLs(hi_mu, evaluationType, "alpha-CLs", workspace_index )
                # print ( "we are at",lo_mu,med_mu,hi_mu,"values at", rt1, rt5, rt10, "scale at", self.scale,"factor at", factor )
                if rt1 < 0.0 and 0.0 < rt10:  # Here's the real while condition
                    break
                if self.alreadyBeenThere:
                    factor = 1 + .75 * (factor - 1)
                    logger.debug("Diminishing rescaling factor")
                if np.isnan(rt1):
                    rt5 = self._CLs(med_mu, evaluationType, "alpha-CLs", workspace_index )
                    if rt5 < 0.0 and rt10 > 0.0:
                        lo_mu = med_mu
                        med_mu = np.sqrt(lo_mu * hi_mu)
                        continue
                    if rt10 < 0.0:  ## also try to increase hi_mu
                        hi_mu = hi_mu + (10.0 - hi_mu) * 0.5
                        med_mu = np.sqrt(lo_mu * hi_mu)
                    nNan += 1
                    self.rescale(factor)
                    continue
                if np.isnan(rt10):
                    rt5 = self._CLs(med_mu, evaluationType, "alpha-CLs", workspace_index )
                    if rt5 > 0.0 and rt1 < 0.0:
                        hi_mu = med_mu
                        med_mu = np.sqrt(lo_mu * hi_mu)
                        continue
                    if rt1 > 0.0:  ## also try to decrease lo_mu
                        lo_mu = lo_mu * 0.5
                        med_mu = np.sqrt(lo_mu * hi_mu)
                    nNan += 1
                    self.rescale(1 / factor)
                    continue
                # Analyzing previous values of wereBoth***
                if rt10 < 0 and rt1 < 0 and wereBothLarge:
                    factor = 1 + (factor - 1) / 2
                    logger.debug("Diminishing rescaling factor")
                if rt10 > 0 and rt1 > 0 and wereBothTiny:
                    factor = 1 + (factor - 1) / 2
                    logger.debug("Diminishing rescaling factor")
                # Preparing next values of wereBoth***
                wereBothTiny = rt10 < 0 and rt1 < 0
                wereBothLarge = rt10 > 0 and rt1 > 0
                # Main rescaling code
                if rt10 < 0.0:
                    self.rescale(factor)
                    continue
                if rt1 > 0.0:
                    self.rescale(1 / factor)
                    continue
            # Finding the root (Brent bracketing part)
            logger.debug( f"Final scale : {self.scale}" )
            ul = findRoot ( self._CLs, lo_mu, hi_mu, args=(evaluationType, "alpha-CLs", workspace_index), rtol=1e-3, xtol=1e-3 )
            endUL = time.time()
            logger.debug( f"getUpperLimitOnMu elapsed time : {endUL-startUL:1.4f} secs" )
            ul = ul * self.scale
            return ul  # self.scale has been updated within self.rescale() method

if __name__ == "__main__":
    ### Needs to be updated
    print("Needs to be updated")
    # C = [
    #     18774.2,
    #     -2866.97,
    #     -5807.3,
    #     -4460.52,
    #     -2777.25,
    #     -1572.97,
    #     -846.653,
    #     -442.531,
    #     -2866.97,
    #     496.273,
    #     900.195,
    #     667.591,
    #     403.92,
    #     222.614,
    #     116.779,
    #     59.5958,
    #     -5807.3,
    #     900.195,
    #     1799.56,
    #     1376.77,
    #     854.448,
    #     482.435,
    #     258.92,
    #     134.975,
    #     -4460.52,
    #     667.591,
    #     1376.77,
    #     1063.03,
    #     664.527,
    #     377.714,
    #     203.967,
    #     106.926,
    #     -2777.25,
    #     403.92,
    #     854.448,
    #     664.527,
    #     417.837,
    #     238.76,
    #     129.55,
    #     68.2075,
    #     -1572.97,
    #     222.614,
    #     482.435,
    #     377.714,
    #     238.76,
    #     137.151,
    #     74.7665,
    #     39.5247,
    #     -846.653,
    #     116.779,
    #     258.92,
    #     203.967,
    #     129.55,
    #     74.7665,
    #     40.9423,
    #     21.7285,
    #     -442.531,
    #     59.5958,
    #     134.975,
    #     106.926,
    #     68.2075,
    #     39.5247,
    #     21.7285,
    #     11.5732,
    # ]
    # nsignal = [x / 100.0 for x in [47, 29.4, 21.1, 14.3, 9.4, 7.1, 4.7, 4.3]]
    # m = PyhfData(
    #     observed=[1964, 877, 354, 182, 82, 36, 15, 11],
    #     backgrounds=[2006.4, 836.4, 350.0, 147.1, 62.0, 26.2, 11.1, 4.7],
    #     covariance=C,
    #     #              third_moment = [ 0.1, 0.02, 0.1, 0.1, 0.003, 0.0001, 0.0002, 0.0005 ],
    #     third_moment=[0.0] * 8,
    #     nsignal=nsignal,
    #     name="ATLAS-SUSY-2018-31 model",
    # )
    # ulComp = PyhfUpperLimitComputer(cl=0.95)
    # # uls = ulComp.getUpperLimitOnMu ( Data ( 15,17.5,3.2,0.00454755 ) )
    # # print ( "uls=", uls )
    # ul_old = 131.828 * sum(
    #     nsignal
    # )  # With respect to the older refernece value one must normalize the xsec
    # print("old ul=", ul_old)
    # ul = ulComp.getUpperLimitOnMu(m)
    # print("ul", ul)
