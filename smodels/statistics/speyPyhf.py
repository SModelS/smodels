#!/usr/bin/env python3

"""
.. module:: speyPyhf
   :synopsis: Code that prepares the input for spey-pyhf, stuff that is
              too big and too pyhf specific for speyTools

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

__all__ = [ "SpeyPyhfData" ]

from smodels.base.smodelsLogging import logger
from smodels.experiment.datasetObj import CombinedDataSet
from typing import Optional

class SpeyPyhfData:
    """
    Holds data for use in pyhf
    :ivar nsignals: signal predictions list divided into sublists, one for each
                    json file
    :ivar inputJsons: list of json instances
    :ivar jsonFiles: optional list of json files
    :ivar nWS: number of workspaces = number of json files
    """
    __slots__ = [ "includeCRs", "nsignals", "inputJson", "jsonFile",
                  "errorFlag", "totalYield", "channelsInfo", "zeroSignalsFlag" ]

    def __init__( self, nsignals : list,
                  inputJson : dict,
                  jsonFile : Optional[str] = None,
                  includeCRs : bool = False ):
        # we dont want to be warned about deprecations within the pyhf code
        import warnings
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        self.includeCRs = includeCRs
        self.nsignals = nsignals  # fb
        self.getTotalYield()
        self.inputJson = inputJson
        self.jsonFile = jsonFile

        self.errorFlag = False
        self.getWSInfo()
        self.checkConsistency()

    def getTotalYield ( self ):
        """ the total yield in all signal regions """
        S = sum ( self.nsignals )
        self.totalYield = S

    def getWSInfo(self):
        """
        Getting information from the json files

        :ivar channelsInfo: list of dictionaries (one dictionary for each json
        file) containing useful information about the json files
            - :key signalRegions: list of dictonaries with 'json path' and
              'size' (number of bins) of the 'signal regions' channels in
              the json files
            - :key otherRegions: list of strings indicating the path to the
              control and validation region channels
        """
        # Identifying the path to the SR and VR channels in the main workspace files
        self.channelsInfo = None  # workspace specifications
        if not isinstance(self.inputJson, dict):
            logger.error("The `inputJson` parameter must be of type dict")
            self.errorFlag = True
            return
        ws = self.inputJson

        wsChannelsInfo = {}
        wsChannelsInfo["signalRegions"] = []
        wsChannelsInfo["otherRegions"] = []
        if "channels" not in ws.keys():
            idx = self.inputJsons.index(ws)
            logger.error (
                f"Json file number {idx} is corrupted (channels are missing)"
            )
            self.channelsInfo = None
            return
        for i_ch, ch in enumerate(ws["channels"]):
            if ch["name"][:2] == "SR":  # if channel name starts with 'SR'
                wsChannelsInfo["signalRegions"].append(
                    {
                        "path": "/channels/"
                        + str(i_ch)
                        + "/samples/0",
                        # Path of the new sample to add (signal prediction)
                        "size": len(ch["samples"][0]["data"]),
                    }
                )  # Number of bins
            else:
                wsChannelsInfo["otherRegions"].append("/channels/" + str(i_ch))
        wsChannelsInfo["otherRegions"].sort(
            key=lambda path: path.split("/")[-1], reverse=True
        )  # Need to sort correctly the paths to the channels to be removed
        self.channelsInfo = wsChannelsInfo

    @classmethod
    def createDataObject ( cls, dataset : CombinedDataSet, nsig : list,
           regionSetName : str ):
        """ an object creator method """

        globalInfo = dataset.globalInfo
        model_tuples = globalInfo.statModels[regionSetName]
        model_tuple = model_tuples[0]
        if "pyhf" not in model_tuple:
            return None # this is not a pyhf model we want here
        jsName = model_tuple[1]
        datasets = []
        regionSet = globalInfo.regionSets [ regionSetName ] # regionSetNames [ jsName ] ]
        # Constructing the list of signals with subsignals matching each json
        nsignals = []
        for region in regionSet:
            datasets.append ( globalInfo.regionMappings[region]["smodels"] )
            regionName = globalInfo.regionMappings[region]["smodels"]
            if regionName == None:
                continue
            if regionName not in nsig:
                logger.debug ( f"sr name {regionName} is not found in {nsig}" )
                continue
                #sys.exit(-1)
            sig = nsig[ regionName ]
            nsignals.append(sig)
        logger.error( f"list of datasets: {datasets}" )
        logger.error( f"jsonFile after filtering: {jsName}" )
        # Loading the jsonFiles, unless we already have them (because we pickled)
        json = globalInfo.cachedModels[jsName]
        return cls( nsignals , json, jsName)

    def checkConsistency(self):
        """
        Check various inconsistencies of the PyhfData attributes

        :param zeroSignalsFlag: boolean identifying if all SRs of a
            single json are empty
        """
        if not isinstance(self.nsignals, list):
            logger.error("The `nsignals` parameter must be of type list")
            self.errorFlag = True

        self.zeroSignalsFlag = list()
        if self.channelsInfo == None:
            return
        wsInfo = self.channelsInfo
        subSig = self.nsignals
        if not isinstance(subSig, list):
            logger.error("The `nsignals` parameter must be a two dimensional list")
            self.errorFlag = True
        nBinsJson = 0
        for sr in wsInfo["signalRegions"]:
            nBinsJson += sr["size"]
        if nBinsJson != len(subSig):
            logger.error(
                f"The number of signals provided is different from the number of bins for json number {self.channelsInfo.index(wsInfo)} and channel number {self.nsignals.index(subSig)}"
            )
            self.errorFlag = True
        allZero = all([s == 0 for s in subSig])
        # Checking if all signals matching this json are zero
        self.zeroSignalsFlag.append(allZero)

    def patchMaker(self):
        """
        Method that creates the list of patches to be applied to the
        `self.inputJsons` workspaces, one for each region given the
        `self.nsignals` and the information available in `self.channelsInfo`
        and the content of the `self.inputJsons` NB: It seems we need to
        include the change of the "modifiers" in the patches as well

        :return: the list of patches, one for each workspace
        """
        if self.channelsInfo == None:
            return None
        info = self.channelsInfo
        subSig = self.nsignals
        patch = []
        for srInfo in info["signalRegions"]:
            nBins = srInfo["size"]
            operator = {}
            operator["op"] = "add"
            operator["path"] = srInfo["path"]
            value = {}
            value["data"] = subSig[:nBins]
            subSig = subSig[nBins:]
            value["modifiers"] = []
            value["modifiers"].append({"data": None, "type": "normfactor",
                                       "name": "mu_SIG"})
            value["modifiers"].append({"data": None, "type": "lumi",
                                       "name": "lumi"})
            value["name"] = "bsm"
            operator["value"] = value
            patch.append(operator)
        if self.includeCRs:
            logger.debug("keeping the CRs")
        else:
            for path in info["otherRegions"]:
                patch.append({"op": "remove", "path": path})
        return patch
    
    def wsMaker(self, apriori=False):
        """
        Apply each region patch (self.patches) to his associated json
            (self.inputJsons) to obtain the complete workspaces
        :param apriori: - If set to `True`: Replace the observation data
            entries of each workspace by the corresponding sum of the
            expected yields \
            - Else: The observed yields put in the workspace are the ones
            written in the corresponfing json dictionary

        :returns: the list of patched workspaces
        """
        if self.patches == None:
            return None
        if self.nWS == 1:
            try:
                wsDict = jsonpatch.apply_patch(self.inputJsons[0], self.patches[0])
                if apriori == True:
                    # Replace the observation data entries by the
                    # corresponding sum of the expected yields
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
                logger.error("The json file is corrupted:\n{}".format(e))
                return None
        else:
            workspaces = []
            for js, patch in zip(self.inputJsons, self.patches):
                wsDict = jsonpatch.apply_patch(js, patch)
                if apriori == True:
                    # Replace the observation data entries by the
                    # corresponding sum of the expected yields
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
