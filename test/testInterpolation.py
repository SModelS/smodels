#!/usr/bin/env python3

"""
.. module:: testInterpolation
   :synopsis: Tests the retrieval of the upper limits, including the PCA
              and the triangulation.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys

sys.path.insert(0, "../")
import unittest
from smodels.experiment.txnameObj import TxNameData
from smodels.tools.physicsUnits import GeV, pb, fb
from databaseLoader import database


class InterpolationTest(unittest.TestCase):
    def testExpected(self):
        expRes = database.getExpResults(
            analysisIDs=["CMS-PAS-SUS-12-026"], datasetIDs=[None], txnames=["T1tttt"]
        )
        txname = expRes[0].datasets[0].txnameList[0]
        m = [[650.0 * GeV, 50.0 * GeV], [650.0 * GeV, 50.0 * GeV]]
        observed = txname.getULFor(m, expected=False)
        expected = txname.getULFor(m, expected=True)
        delta = abs(((observed - expected) / observed).asNumber())
        self.assertTrue(delta > 0.55 and delta < 0.60)

    def testExpectedFails(self):
        expRes = database.getExpResults(
            analysisIDs=["ATLAS-SUSY-2013-05"], datasetIDs=[None], txnames=["T2bb"]
        )
        txname = expRes[0].datasets[0].txnameList[0]
        m = [[650.0 * GeV, 50.0 * GeV], [650.0 * GeV, 50.0 * GeV]]
        expected = txname.getULFor(m, expected=True)
        self.assertTrue(expected is None)

    def testInterpolation(self):
        # print database
        expRes = database.getExpResults(
            analysisIDs=["ATLAS-SUSY-2013-05"], datasetIDs=[None], txnames=["T2bb"]
        )
        # expRes=listOfExpRes[0]   # ATLAS-SUSY-2013-05
        txname = expRes[0].datasets[0].txnameList[0]  # T2bb
        result = txname.txnameData.getValueFor(
            [[300.0 * GeV, 100.0 * GeV], [300.0 * GeV, 100.0 * GeV]]
        )
        self.assertAlmostEqual(result.asNumber(pb), 0.162457)
        result = txname.txnameData.getValueFor(
            [[300.0 * GeV, 125.0 * GeV], [300.0 * GeV, 125.0 * GeV]]
        )
        self.assertAlmostEqual(result.asNumber(pb), 0.237745)

    def test6D(self):
        # print database
        expRes = database.getExpResults(analysisIDs=["ATLAS-SUSY-2013-05"], txnames=["T6bbWW"])
        # expRes=listOfExpRes[0]   # ATLAS-SUSY-2013-05
        txname = expRes[0].datasets[0].txnameList[0]  # T6bbWW
        result = txname.txnameData.getValueFor(
            [[300.0 * GeV, 105.0 * GeV, 100.0 * GeV], [300.0 * GeV, 105.0 * GeV, 100.0 * GeV]]
        )
        self.assertAlmostEqual(result.asNumber(pb), 0.176266)
        result = txname.txnameData.getValueFor(
            [[300.0 * GeV, 270.0 * GeV, 200.0 * GeV], [300.0 * GeV, 270.0 * GeV, 200.0 * GeV]]
        )
        self.assertAlmostEqual(result.asNumber(pb), 87.0403)
        result = txname.txnameData.getValueFor(
            [[300.0 * GeV, 270.0 * GeV, 200.0 * GeV], [300.0 * GeV, 271.0 * GeV, 200.0 * GeV]]
        )
        self.assertAlmostEqual(result.asNumber(pb), 88.6505675)

    def testOutsidePlane(self):
        expRes = database.getExpResults(analysisIDs=["ATLAS-SUSY-2013-05"], txnames=["T2bb"])
        # expRes=listOfExpRes[0]   # ATLAS-SUSY-2013-05
        txname = expRes[0].datasets[0].txnameList[0]  # T6bbWW
        result = txname.txnameData.getValueFor(
            [[300.0 * GeV, 127.0 * GeV], [300.0 * GeV, 127.5 * GeV]]
        )
        self.assertAlmostEqual(result.asNumber(pb), 0.24452092)
        result = txname.txnameData.getValueFor(
            [[600.0 * GeV, 120.0 * GeV], [600.0 * GeV, 130.0 * GeV]]
        )
        self.assertAlmostEqual(result.asNumber(pb), 0.0197154)
        result = txname.txnameData.getValueFor(
            [[300.0 * GeV, 120.0 * GeV], [300.0 * GeV, 130.0 * GeV]]
        )
        self.assertTrue(result == None)

    def testWithDirectData(self):
        data = [
            [[[150.0 * GeV, 50.0 * GeV], [150.0 * GeV, 50.0 * GeV]], 3.0 * fb],
            [[[200.0 * GeV, 100.0 * GeV], [200.0 * GeV, 100.0 * GeV]], 5.0 * fb],
            [[[300.0 * GeV, 100.0 * GeV], [300.0 * GeV, 100.0 * GeV]], 10.0 * fb],
            [[[300.0 * GeV, 150.0 * GeV], [300.0 * GeV, 150.0 * GeV]], 13.0 * fb],
            [[[300.0 * GeV, 200.0 * GeV], [300.0 * GeV, 200.0 * GeV]], 15.0 * fb],
            [[[300.0 * GeV, 250.0 * GeV], [300.0 * GeV, 250.0 * GeV]], 20.0 * fb],
            [[[400.0 * GeV, 100.0 * GeV], [400.0 * GeV, 100.0 * GeV]], 8.0 * fb],
            [[[400.0 * GeV, 150.0 * GeV], [400.0 * GeV, 150.0 * GeV]], 10.0 * fb],
            [[[400.0 * GeV, 200.0 * GeV], [400.0 * GeV, 200.0 * GeV]], 12.0 * fb],
            [[[400.0 * GeV, 250.0 * GeV], [400.0 * GeV, 250.0 * GeV]], 15.0 * fb],
            [[[400.0 * GeV, 300.0 * GeV], [400.0 * GeV, 300.0 * GeV]], 17.0 * fb],
            [[[400.0 * GeV, 350.0 * GeV], [400.0 * GeV, 350.0 * GeV]], 19.0 * fb],
        ]
        txnameData = TxNameData(data, "upperLimit", sys._getframe().f_code.co_name)
        result = txnameData.getValueFor([[300.0 * GeV, 125.0 * GeV], [300.0 * GeV, 125.0 * GeV]])
        self.assertAlmostEqual(result.asNumber(pb), 0.0115)

    def testEfficiencyMaps(self):
        data = [
            [[[150.0 * GeV, 50.0 * GeV], [150.0 * GeV, 50.0 * GeV]], 0.03],
            [[[200.0 * GeV, 100.0 * GeV], [200.0 * GeV, 100.0 * GeV]], 0.05],
            [[[300.0 * GeV, 100.0 * GeV], [300.0 * GeV, 100.0 * GeV]], 0.10],
            [[[300.0 * GeV, 150.0 * GeV], [300.0 * GeV, 150.0 * GeV]], 0.13],
            [[[300.0 * GeV, 200.0 * GeV], [300.0 * GeV, 200.0 * GeV]], 0.15],
            [[[300.0 * GeV, 250.0 * GeV], [300.0 * GeV, 250.0 * GeV]], 0.20],
            [[[400.0 * GeV, 100.0 * GeV], [400.0 * GeV, 100.0 * GeV]], 0.08],
            [[[400.0 * GeV, 150.0 * GeV], [400.0 * GeV, 150.0 * GeV]], 0.10],
            [[[400.0 * GeV, 200.0 * GeV], [400.0 * GeV, 200.0 * GeV]], 0.12],
            [[[400.0 * GeV, 250.0 * GeV], [400.0 * GeV, 250.0 * GeV]], 0.15],
            [[[400.0 * GeV, 300.0 * GeV], [400.0 * GeV, 300.0 * GeV]], 0.17],
            [[[400.0 * GeV, 350.0 * GeV], [400.0 * GeV, 350.0 * GeV]], 0.19],
        ]
        txnameData = TxNameData(data, "efficiencyMap", sys._getframe().f_code.co_name)
        result = txnameData.getValueFor([[300.0 * GeV, 125.0 * GeV], [300.0 * GeV, 125.0 * GeV]])
        self.assertAlmostEqual(result, 0.115)

    def testOutsideConvexHull(self):
        data = [
            [[[150.0 * GeV, 50.0 * GeV], [150.0 * GeV, 50.0 * GeV]], 0.03],
            [[[200.0 * GeV, 100.0 * GeV], [200.0 * GeV, 100.0 * GeV]], 0.05],
            [[[300.0 * GeV, 100.0 * GeV], [300.0 * GeV, 100.0 * GeV]], 0.10],
            [[[300.0 * GeV, 150.0 * GeV], [300.0 * GeV, 150.0 * GeV]], 0.13],
            [[[300.0 * GeV, 200.0 * GeV], [300.0 * GeV, 200.0 * GeV]], 0.15],
            [[[300.0 * GeV, 250.0 * GeV], [300.0 * GeV, 250.0 * GeV]], 0.20],
            [[[400.0 * GeV, 100.0 * GeV], [400.0 * GeV, 100.0 * GeV]], 0.08],
            [[[400.0 * GeV, 150.0 * GeV], [400.0 * GeV, 150.0 * GeV]], 0.10],
            [[[400.0 * GeV, 200.0 * GeV], [400.0 * GeV, 200.0 * GeV]], 0.12],
            [[[400.0 * GeV, 250.0 * GeV], [400.0 * GeV, 250.0 * GeV]], 0.15],
            [[[400.0 * GeV, 300.0 * GeV], [400.0 * GeV, 300.0 * GeV]], 0.17],
            [[[400.0 * GeV, 350.0 * GeV], [400.0 * GeV, 350.0 * GeV]], 0.19],
        ]
        txnameData = TxNameData(data, "efficiencyMap", sys._getframe().f_code.co_name)
        result = txnameData.getValueFor([[300.0 * GeV, 125.0 * GeV], [300.0 * GeV, 123.0 * GeV]])
        self.assertAlmostEqual(result, 0.1144)

    def testOutsideConvexHull2(self):
        data = [
            [[[150.0 * GeV, 50.0 * GeV], [150.0 * GeV, 50.0 * GeV]], 0.03],
            [[[200.0 * GeV, 100.0 * GeV], [200.0 * GeV, 100.0 * GeV]], 0.05],
            [[[300.0 * GeV, 100.0 * GeV], [300.0 * GeV, 100.0 * GeV]], 0.10],
            [[[300.0 * GeV, 150.0 * GeV], [300.0 * GeV, 150.0 * GeV]], 0.13],
            [[[300.0 * GeV, 200.0 * GeV], [300.0 * GeV, 200.0 * GeV]], 0.15],
            [[[300.0 * GeV, 250.0 * GeV], [300.0 * GeV, 250.0 * GeV]], 0.20],
            [[[400.0 * GeV, 100.0 * GeV], [400.0 * GeV, 100.0 * GeV]], 0.08],
            [[[400.0 * GeV, 150.0 * GeV], [400.0 * GeV, 150.0 * GeV]], 0.10],
            [[[400.0 * GeV, 200.0 * GeV], [400.0 * GeV, 200.0 * GeV]], 0.12],
            [[[400.0 * GeV, 250.0 * GeV], [400.0 * GeV, 250.0 * GeV]], 0.15],
            [[[400.0 * GeV, 300.0 * GeV], [400.0 * GeV, 300.0 * GeV]], 0.17],
            [[[400.0 * GeV, 350.0 * GeV], [400.0 * GeV, 350.0 * GeV]], 0.19],
        ]
        txnameData = TxNameData(data, "efficiencyMap", sys._getframe().f_code.co_name)
        result = txnameData.getValueFor([[300.0 * GeV, 125.0 * GeV], [300.0 * GeV, 100.0 * GeV]])
        self.assertEqual(result, None)


if __name__ == "__main__":
    unittest.main()
