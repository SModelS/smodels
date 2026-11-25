#!/usr/bin/env python3

"""
.. module:: testRoofit
   :synopsis: Test the ROOT/spey-roofit statistical model plumbing.
"""

import os
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, "../")

from smodels.base.physicsUnits import fb
from smodels.experiment.infoObj import Info
from smodels.statistics.speyTools import SpeyRetriever
from smodels.statistics.statsTools import CompRetriever


class _GlobalInfo:
    id = "test-roofit"
    lumi = 1. / fb
    regionMappings = {
        "SR1": {"smodels": "SR1", "roofit": "channelA[0]"},
        "SR2": {"smodels": "SR2", "roofit": "channelA[1]"},
        "CR1": {"smodels": None, "roofit": "channelCR"},
    }
    regionSets = {"all": ["SR1", "SR2", "CR1"]}
    statModels = {"all": [("roofit", "model.root")]}
    cachedModels = {"model.root": b"ROOTFILE"}


class _DataSet:
    globalInfo = _GlobalInfo


class RoofitTest(unittest.TestCase):

    def testRegionMappingDefaults(self):
        info = Info()
        regions = info.canonizeRegions({
            "SR1": {"smodels": "SR1", "pyhf": "pyhfChannel"},
            "SR2": {"smodels": "SR2"},
        })
        self.assertEqual(regions["SR1"]["roofit"], "pyhfChannel")
        self.assertEqual(regions["SR2"]["roofit"], "SR2")

    def testCacheRootModelAsBytes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with open(os.path.join(tmpdir, "model.root"), "wb") as f:
                f.write(b"ROOTFILE")
            info = Info()
            info.path = os.path.join(tmpdir, "globalInfo.txt")
            info.regionSets = {"all": ["SR1"]}
            info.statModels = {"all": [("full_roofit", "model.root")]}
            info.cacheStatsModels()
            self.assertEqual(info.cachedModels["model.root"], b"ROOTFILE")

    def testForRoofitPassesRootFileToSpey(self):
        kwargs = {}

        def fakeBackend(**kw):
            kwargs.update(kw)
            return "model"

        with patch("smodels.statistics.speyTools.get_backend",
                   return_value=fakeBackend) as getBackend:
            facade = SpeyRetriever.forRoofit("all", _DataSet(),
                                             {"SR1": 2., "SR2": 3.})
        getBackend.assert_called_once_with("root.roofit")
        self.assertEqual(kwargs["analysis"], "test-roofit")
        self.assertEqual(kwargs["signal_yields"], {"channelA": [2., 3.]})
        self.assertTrue(os.path.isfile(kwargs["workspace"]))
        with open(kwargs["workspace"], "rb") as f:
            self.assertEqual(f.read(), b"ROOTFILE")
        self.assertEqual(facade.dataType, "roofit")
        facade.roofitData.cleanup()

    def testNativeRetrieverDelegates(self):
        sentinel = object()
        with patch.object(SpeyRetriever, "forRoofit",
                          return_value=sentinel) as forRoofit:
            ret = CompRetriever.forRoofit("all", _DataSet(), {"SR1": 1.})
        self.assertIs(ret, sentinel)
        forRoofit.assert_called_once()


if __name__ == "__main__":
    unittest.main()
