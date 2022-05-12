#!/usr/bin/env python3

"""
.. module:: testExtrapolation
   :synopsis: Tests the small extrapolation that we perform for
              sub-dimensional data.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
import sys

sys.path.insert(0, "../")
import unittest
import numpy
from smodels.experiment.txnameObj import TxNameData
from smodels.tools.physicsUnits import GeV, eV, fb, MeV, pb, keV
import unum


def pprint(energy):
    """return energy in pretty format"""
    unum.Unum.VALUE_FORMAT = "%.1f"
    unum.Unum.UNIT_FORMAT = "%s"
    units = [eV, keV, MeV, GeV]
    for unit in units:
        if energy.asNumber(unit) < 1000.0:
            return "%s" % energy.asUnit(unit)
    return energy


class ExtrapolationTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(ExtrapolationTest, self).__init__(*args, **kwargs)
        data = [
            [[[150.0 * GeV, 50.0 * GeV], [150.0 * GeV, 50.0 * GeV]], 3.0 * fb],
            [[[200.0 * GeV, 100.0 * GeV], [200.0 * GeV, 100.0 * GeV]], 5.0 * fb],
            [[[300.0 * GeV, 200.0 * GeV], [300.0 * GeV, 200.0 * GeV]], 15.0 * fb],
            [[[400.0 * GeV, 300.0 * GeV], [400.0 * GeV, 300.0 * GeV]], 17.0 * fb],
        ]
        self.txnameData = TxNameData(data, "upperLimit", sys._getframe().f_code.co_name, 0.05)

    def tryWith(self, masses):
        return self.txnameData.getValueFor(masses)

    def testWithDirectData(self):
        result = self.tryWith([[275.0 * GeV, 175.0 * GeV], [275.0 * GeV, 175.0 * GeV]])

        self.assertAlmostEqual(result.asNumber(pb), 0.0125)
        eps = 1 * keV
        result = self.tryWith([[275.0 * GeV, 175.0 * GeV + eps], [275.0 * GeV + eps, 175.0 * GeV]])

        self.assertAlmostEqual(result.asNumber(pb), 0.0125)

        result = self.tryWith([[275.0 * GeV, 185.0 * GeV], [275.0 * GeV, 165.0 * GeV]])
        self.assertTrue(result == None)

    def show(self):
        # txnameObj.nonZeroEps = 1e+4
        # txnameObj.simplexTolerance = 1e-2
        for eps in numpy.arange(2, 9, 0.3):
            e = (10**eps) * eV
            masses = [[275.0 * GeV + e, 175.0 * GeV], [275.0 * GeV, 175.0 * GeV - e]]
            print("%s: %s" % (pprint(e), self.tryWith(masses)))


if __name__ == "__main__":
    # ExtrapolationTest("testWithDirectData").show()
    unittest.main()
