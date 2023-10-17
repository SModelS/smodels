#!/usr/bin/env python3

"""
.. module:: speyPyhf
   :synopsis: Code that prepares the input for spey-pyhf, that is
              too big and too specific for speyTools

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from smodels.tools.smodelsLogging import logger
import os


class SpeyPyhfData:
    """
    Holds data for use in pyhf
    :ivar nsignals: signal predictions list divided into sublists, one for each
                    json file
    :ivar inputJsons: list of json instances
    :ivar jsonFiles: optional list of json files
    :ivar nWS: number of workspaces = number of json files
    """

    def __init__(self, nsignals : list, inputJsons, jsonFiles=None,
                 includeCRs : bool = False ):
        # we dont want to be warned about deprecations within the pyhf code
        import warnings
        warnings.filterwarnings("ignore", category=DeprecationWarning) 
        self.includeCRs = includeCRs
        self.nsignals = nsignals  # fb
        self.getTotalYield()
        self.inputJsons = inputJsons
        self.cached_likelihoods = {}  ## cache of likelihoods (actually twice_nlls)
        self.cached_lmaxes = {}  # cache of lmaxes (actually twice_nlls)
        self.cachedULs = {False: {}, True: {}, "posteriori": {}}
        self.jsonFiles = jsonFiles
        self.combinations = None
        if jsonFiles != None:
            self.combinations = [os.path.splitext(os.path.basename(js))[0] for js in jsonFiles]

        self.nWS = len(inputJsons)
        self.errorFlag = False
        self.getWSInfo()
        self.checkConsistency()

    def getTotalYield ( self ):
        """ the total yield in all signal regions """
        S = sum ( [ sum(x) for x in self.nsignals ] )
        self.totalYield = S

    def getWSInfo(self):
        """
        Getting information from the json files

        :ivar channelsInfo: list of dictionaries (one dictionary for each json file) containing useful information about the json files
            - :key signalRegions: list of dictonaries with 'json path' and 'size' (number of bins) of the 'signal regions' channels in the json files
            - :key otherRegions: list of strings indicating the path to the control and validation region channels
        """
        # Identifying the path to the SR and VR channels in the main workspace files
        self.channelsInfo = []  # workspace specifications
        if not isinstance(self.inputJsons, list):
            logger.error("The `inputJsons` parameter must be of type list")
            self.errorFlag = True
            return
        for ws in self.inputJsons:
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
            for i_ch, ch in enumerate(ws["channels"]):
                if ch["name"][:2] == "SR":  # if channel name starts with 'SR'
                    wsChannelsInfo["signalRegions"].append(
                        {
                            "path": "/channels/"
                            + str(i_ch)
                            + "/samples/0",  # Path of the new sample to add (signal prediction)
                            "size": len(ch["samples"][0]["data"]),
                        }
                    )  # Number of bins
                else:
                    wsChannelsInfo["otherRegions"].append("/channels/" + str(i_ch))
            wsChannelsInfo["otherRegions"].sort(
                key=lambda path: path.split("/")[-1], reverse=True
            )  # Need to sort correctly the paths to the channels to be removed
            self.channelsInfo.append(wsChannelsInfo)

    @classmethod
    def createDataObject ( cls, dataset, nsig : list ):
        """ an object creator method """
        jsonFiles = dataset.globalInfo.jsonFiles

        globalInfo = dataset.globalInfo
        jsonFiles = [js for js in globalInfo.jsonFiles]
        jsons = globalInfo.jsons.copy()
        # datasets = [ds.getID() for ds in dataset._datasets]
        datasets = [ds.getID() for ds in dataset.origdatasets]
        # Filtering the json files by looking at the available datasets
        for jsName in globalInfo.jsonFiles:
            if all([ds not in globalInfo.jsonFiles[jsName] for ds in datasets]):
                # No datasets found for this json combination
                jsIndex = jsonFiles.index(jsName)
                jsonFiles.pop(jsIndex)
                jsons.pop(jsIndex)
                continue
            if not all([ds in datasets for ds in globalInfo.jsonFiles[jsName]]):
                # Some SRs are missing for this json combination
                logger.error( "Wrong json definition in globalInfo.jsonFiles for json : %s" % jsName)
        logger.debug("list of datasets: {}".format(datasets))
        logger.debug("jsonFiles after filtering: {}".format(jsonFiles))
        # Constructing the list of signals with subsignals matching each json
        nsignals = list()
        for jsName in jsonFiles:
            subSig = list()
            for srName in globalInfo.jsonFiles[jsName]:
                try:
                    index = datasets.index(srName)
                except ValueError:
                    line = (
                        f"{srName} signal region provided in globalInfo is not in the list of datasets, {jsName}:{','.join(datasets)}"
                    )
                    raise ValueError(line)
                sig = nsig[index]
                subSig.append(sig)
            nsignals.append(subSig)
        # Loading the jsonFiles, unless we already have them (because we pickled)
        return cls(nsignals, jsons, jsonFiles)

    def checkConsistency(self):
        """
        Check various inconsistencies of the PyhfData attributes

        :param zeroSignalsFlag: boolean identifying if all SRs of a single json are empty
        """
        if not isinstance(self.nsignals, list):
            logger.error("The `nsignals` parameter must be of type list")
            self.errorFlag = True
        if self.nWS != len(self.nsignals):
            logger.error(
                "The number of subsignals provided is different from the number of json files"
            )
            self.errorFlag = True
        self.zeroSignalsFlag = list()
        if self.channelsInfo == None:
            return
        for wsInfo, subSig in zip(self.channelsInfo, self.nsignals):
            if not isinstance(subSig, list):
                logger.error("The `nsignals` parameter must be a two dimensional list")
                self.errorFlag = True
            nBinsJson = 0
            for sr in wsInfo["signalRegions"]:
                nBinsJson += sr["size"]
            if nBinsJson != len(subSig):
                logger.error(
                    "The number of signals provided is different from the number of bins for json number {} and channel number {}".format(
                        self.channelsInfo.index(wsInfo), self.nsignals.index(subSig)
                    )
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
        nsignals = self.nsignals
        # Constructing the patches to be applied on the main workspace files
        patches = []
        for ws, info, subSig in zip(self.inputJsons, self.channelsInfo, self.nsignals):
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
                value["modifiers"].append({"data": None, "type": "normfactor", "name": "mu_SIG"})
                value["modifiers"].append({"data": None, "type": "lumi", "name": "lumi"})
                value["name"] = "bsm"
                operator["value"] = value
                patch.append(operator)
            if self.includeCRs:
                logger.debug("keeping the CRs")
            else:
                for path in info["otherRegions"]:
                    patch.append({"op": "remove", "path": path})
            patches.append(patch)
        return patches

    def wsMaker(self, apriori=False):
        """
        Apply each region patch (self.patches) to his associated json (self.inputJsons) to obtain the complete workspaces
        :param apriori: - If set to `True`: Replace the observation data entries of each workspace by the corresponding sum of the expected yields \
                        - Else: The observed yields put in the workspace are the ones written in the corresponfing json dictionary

        :returns: the list of patched workspaces
        """
        if self.patches == None:
            return None
        if self.nWS == 1:
            try:
                wsDict = jsonpatch.apply_patch(self.inputJsons[0], self.patches[0])
                if apriori == True:
                    # Replace the observation data entries by the corresponding sum of the expected yields
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
                    # Replace the observation data entries by the corresponding sum of the expected yields
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

