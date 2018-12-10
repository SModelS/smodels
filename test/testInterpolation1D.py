#!/usr/bin/env python3

"""
.. module:: testInterpolation
   :synopsis: Tests the retrieval of the upper limits, including the PCA
              and the triangulation.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.experiment.txnameObj import TxNameData
from smodels.tools.physicsUnits import GeV, pb
import numpy as np

class Interpolation1DTest(unittest.TestCase):

    def testWithDirectData(self):
        data = [[[[1.5640E+02*GeV],[1.5640E+02*GeV]],3.1127E-03*pb],
                [[[1.9728E+02*GeV],[1.9728E+02*GeV]],1.6769E-03*pb],
                [[[2.4632E+02*GeV],[2.4632E+02*GeV]],9.2590E-04*pb],
                [[[3.0899E+02*GeV],[3.0899E+02*GeV]],8.4880E-04*pb],
                [[[4.9155E+02*GeV],[4.9155E+02*GeV]],6.3010E-04*pb],
                [[[5.5967E+02*GeV],[5.5967E+02*GeV]],5.9950E-04*pb],
                [[[6.5504E+02*GeV],[6.5504E+02*GeV]],5.6320E-04*pb],
                [[[7.4496E+02*GeV],[7.4496E+02*GeV]],5.2260E-04*pb],
                [[[8.6757E+02*GeV],[8.6757E+02*GeV]],5.1580E-04*pb],
                [[[1.0283E+03*GeV],[1.0283E+03*GeV]],5.0900E-04*pb],
                [[[1.2191E+03*GeV],[1.2191E+03*GeV]],4.8380E-04*pb],
                [[[1.4098E+03*GeV],[1.4098E+03*GeV]],5.1410E-04*pb],
                [[[1.6005E+03*GeV],[1.6005E+03*GeV]],5.8110E-04*pb]]
        
        txnameData=TxNameData ( data, "upperLimit",
                sys._getframe().f_code.co_name )

        self.assertEqual(txnameData.dimensionality,1)

        #Tranformation to "PCA frame"        
        dataShift = (data[0][0][0][0]+data[-1][0][0][0])/2.
        dataShift = dataShift.asNumber(GeV)
        
        #Check inside the grid:
        result=txnameData.getValueFor([[ 1100.*GeV]]*2)
        self.assertAlmostEqual( result.asNumber(pb),0.00049953) 
        #Check outside the grid:
        result=txnameData.getValueFor([[ 2000.*GeV]]*2)
        self.assertEqual( result,None)
        
        #Check class methods:
        isimplex = txnameData.tri.find_simplex(np.array([350.-dataShift]))
        self.assertEqual(isimplex,3)
        simplex = txnameData.tri.simplices[isimplex]
        self.assertEqual(min(simplex),3)
        self.assertEqual(max(simplex),4)

if __name__ == "__main__":
    unittest.main()
