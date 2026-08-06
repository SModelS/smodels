#!/usr/bin/env python3

"""
.. module:: testStatistics
   :synopsis: Tests the statistics functionality.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Jory Sonneveld <jory@opmijnfiets.nl>

"""
import sys

sys.path.insert(0, "../")
import unittest

# from smodels.tools import statistics
from smodels.base import runtime

class ExperimentalFeaturesTest(unittest.TestCase):
    def testFeatures ( self ):
        """ a few simple tests for the experimentalFeatures feature """
        a = runtime.experimentalFeature ( "inexistentparameter" )
        self.assertTrue ( a == None )
        b = runtime.experimentalFeature ( "truncatedgaussians" )
        self.assertTrue ( b == False )

    def testParameterIni ( self ):
        from smodels.base.smodelsLogging import setLogLevel
        setLogLevel ( "error" )
        inifile = "testParameters_experimentalFeatures.ini"
        from smodels.matching.modelTester import getParameters, setExperimentalFeatures
        params = getParameters ( inifile )
        setExperimentalFeatures ( params )
        b = runtime.experimentalFeature ( "truncatedgaussians" )
        self.assertTrue ( b == True )
        a = runtime.experimentalFeature ( "inexistentparameter" )
        self.assertTrue ( a == None )


if __name__ == "__main__":
    unittest.main()
