#!/usr/bin/env python3

"""
.. module:: speyRoofit
   :synopsis: Code that prepares ROOT statistical models for spey-roofit.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

__all__ = [ "SpeyRoofitData" ]

import os
import re
import tempfile
from typing import Optional

from smodels.base.smodelsLogging import logger
from smodels.experiment.datasetObj import CombinedDataSet
from smodels.statistics.exceptions import SModelSStatisticsError as SModelSError


def _splitBinName(name: str) -> tuple:
    """Split a channel name into the channel and a possible bin index.

    ``"SR_meff[2]"`` becomes ``("SR_meff", 2)``; a name without a bin
    suffix becomes ``(name, None)``.
    """
    match = re.match(r"^(.*)\[(\d+)\]$", name)
    if match:
        return match.group(1), int(match.group(2))
    return name, None


def _channelName(region: dict) -> Optional[str]:
    """Return the channel name of a region for the roofit backend.

    The dedicated ``roofit`` key is used first. For backwards
    compatibility, the pyhf, HS3 and SModelS names are used as fallbacks.
    """
    for key in ("roofit", "pyhf", "hs3", "smodels"):
        value = region.get(key)
        if value is not None:
            return value
    return None


class SpeyRoofitData:
    """Holds the ROOT model and signal yields for spey-roofit.

    :ivar inputRoot: the ROOT file content as cached by :class:`Info`
    :ivar rootFile: path of the temporary ROOT file handed to spey-roofit
    :ivar rootFileName: name of the model in ``statModels``
    :ivar signal_yields: dictionary ``{channel: [yield per bin]}``
    :ivar nsignals: signal yields keyed by SModelS region name
    :ivar regions: the region mappings of the region set
    :ivar totalYield: total signal yield in all included regions
    """

    __slots__ = [ "inputRoot", "rootFile", "rootFileName", "signal_yields",
                  "nsignals", "regions", "channels", "totalYield",
                  "includeCRs", "errorFlag", "globalInfo" ]

    def __init__(self, nsignals: dict, inputRoot: bytes,
                 rootFileName: Optional[str] = None,
                 regions: Optional[list] = None,
                 includeCRs: bool = False,
                 globalInfo: Optional[object] = None):
        if not isinstance(inputRoot, (bytes, bytearray)):
            raise SModelSError(
                f"The ROOT statistical model {rootFileName} must be cached "
                f"as bytes, but is {type(inputRoot)}"
            )
        self.inputRoot = bytes(inputRoot)
        self.nsignals = nsignals
        self.regions = regions or []
        self.rootFileName = rootFileName
        self.rootFile = None
        self.includeCRs = includeCRs
        self.globalInfo = globalInfo
        self.errorFlag = False

        # spey-roofit takes a file name, so materialise the cached model.
        # The file is removed when this object is garbage collected; the
        # SpeyModelFacade keeps a reference for as long as the model lives.
        fd, path = tempfile.mkstemp(prefix="smodels_spey_roofit_", suffix=".root")
        os.close(fd)
        try:
            with open(path, "wb") as f:
                f.write(self.inputRoot)
        except Exception:
            try:
                os.unlink(path)
            except OSError:
                pass
            raise
        self.rootFile = path

        self.signal_yields = self.getSignalYields()

    def getSignalYields(self) -> dict:
        """Translate SModelS yields into spey-roofit channel yields."""
        signals = {}
        total = 0.0
        for region in self.regions:
            smodelsName = region.get("smodels")
            if smodelsName is None:
                continue
            if region.get("type") == "CR" and not self.includeCRs:
                continue
            channel = _channelName(region)
            if channel is None:
                logger.error(f"Region {region} has no roofit/pyhf/HS3/SModelS channel name")
                self.errorFlag = True
                continue
            channel, binIndex = _splitBinName(channel)
            label = region.get("_label", smodelsName)
            yield_ = float(self.nsignals.get(smodelsName,
                                             self.nsignals.get(label, 0.0)))
            values = signals.setdefault(channel, [])
            if binIndex is None:
                # Either one bin per channel, or one region per bin in the
                # order given by regionMappings.
                values.append(yield_)
            else:
                if binIndex >= len(values):
                    values.extend([0.0] * (binIndex + 1 - len(values)))
                values[binIndex] = yield_
            total += yield_
        self.channels = list(signals.keys())
        self.totalYield = total
        return signals

    def cleanup(self):
        """Remove the temporary ROOT file, if it still exists."""
        path = getattr(self, "rootFile", None)
        if path is None:
            return
        self.rootFile = None
        try:
            os.unlink(path)
        except OSError:
            pass

    def __del__(self):
        try:
            self.cleanup()
        except Exception:
            pass

    @classmethod
    def createDataObject(cls, dataset: CombinedDataSet, nsig: dict,
                         regionSetName: str) -> Optional["SpeyRoofitData"]:
        """Create the data object from a dataset and a region set.

        :param dataset: CombinedDataSet holding the cached ROOT model
        :param nsig: signal yields keyed by SModelS region name
        :param regionSetName: name of the region set to be combined
        :returns: the data object, or None if this is not a ROOT model
        """
        globalInfo = dataset.globalInfo
        if regionSetName not in globalInfo.statModels:
            raise SModelSError(
                f"regionSet {regionSetName} not found in statModels for "
                f"dataset {globalInfo.id}"
            )
        model_tuples = globalInfo.statModels[regionSetName]
        if not model_tuples:
            return None
        model_type, model_filename = model_tuples[0]
        if model_type not in ("roofit", "full_roofit"):
            return None
        if model_filename not in globalInfo.cachedModels:
            raise SModelSError(
                f"model {model_filename} has not been cached in dataset "
                f"{globalInfo.id}"
            )
        regions = []
        for label in globalInfo.regionSets[regionSetName]:
            if label not in globalInfo.regionMappings:
                raise SModelSError(
                    f"SR {label} defined in regionSet {regionSetName} not "
                    f"found in regionMappings for dataset {globalInfo.id}"
                )
            region = dict(globalInfo.regionMappings[label])
            region["_label"] = label
            regions.append(region)
        includeCRs = getattr(globalInfo, "includeCRs", False)
        return cls(nsig, globalInfo.cachedModels[model_filename],
                   rootFileName=model_filename, regions=regions,
                   includeCRs=includeCRs, globalInfo=globalInfo)
