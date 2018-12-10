#!/usr/bin/env python3

"""
.. module:: testExtrapolation
   :synopsis: Tests the small extrapolation that we perform for
              sub-dimensional data.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
import sys
sys.path.insert(0,"../")
import unittest
import numpy
from smodels.experiment.txnameObj import TxNameData
from smodels.tools.physicsUnits import GeV, eV, fb, MeV, pb, keV
import unum

def pprint ( energy ):
    """ return energy in pretty format """
    unum.Unum.VALUE_FORMAT = "%.1f"
    unum.Unum.UNIT_FORMAT = "%s"
    units = [ eV, keV, MeV, GeV ]
    for unit in units:
        if energy.asNumber ( unit ) < 1000.:
            return "%s" % energy.asUnit ( unit )
    return energy

class ExtrapolationTest(unittest.TestCase):
    def __init__ ( self, *args, **kwargs ):
        super(ExtrapolationTest, self).__init__(*args, **kwargs)
        data = [ [ [[ 150.*GeV, 50.*GeV], [ 150.*GeV, 50.*GeV] ], 3.*fb ],
             [ [[ 200.*GeV,100.*GeV], [ 200.*GeV,100.*GeV] ],  5.*fb ],
             [ [[ 300.*GeV,200.*GeV], [ 300.*GeV,200.*GeV] ], 15.*fb ],
             [ [[ 400.*GeV,300.*GeV], [ 400.*GeV,300.*GeV] ], 17.*fb ] ]
        self.txnameData = TxNameData( data, "upperLimit",
                sys._getframe().f_code.co_name, .05 )

    def tryWith ( self, masses ):
        return self.txnameData.getValueFor( masses )

    def testWithDirectData(self):
        result=self.tryWith ([[ 275.*GeV,175.*GeV], [ 275.*GeV,175.*GeV] ])

        self.assertAlmostEqual( result.asNumber(pb),0.0125 )
        eps = 1 * keV
        result=self.tryWith([[ 275.*GeV,175.*GeV + eps],
                             [ 275.*GeV + eps ,175.*GeV] ])
        
        self.assertAlmostEqual( result.asNumber(pb),0.0125 )

        result=self.tryWith([[ 275.*GeV,185.*GeV], [ 275.*GeV,165.*GeV] ])
        self.assertTrue ( result == None )

    def show ( self ):
        #txnameObj.nonZeroEps = 1e+4
        #txnameObj.simplexTolerance = 1e-2
        for eps in numpy.arange ( 2, 9, .3 ):
            e = (10**eps) * eV
            masses = [[ 275.*GeV + e, 175.*GeV ], [ 275.*GeV,175.*GeV - e ]]
            print ( "%s: %s" % ( pprint ( e ), self.tryWith ( masses ) ))

if __name__ == "__main__":
    # ExtrapolationTest("testWithDirectData").show()
    unittest.main()
