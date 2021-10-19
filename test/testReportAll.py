#!/usr/bin/env python3

"""
.. module:: testReportAll
   :synopsis: Tests the reporting of all datasets

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os
import unittest
from unitTestHelpers import equalObjs, runMain, importModule

class RunSModelSTest(unittest.TestCase):
    definingRun = False ## meant only to adapt to changes in output format
    ## use with super great care!!

    def removeOutputs( self, f ):
        """ remove cruft outputfiles """
        for i in [ f, f.replace(".py",".pyc") ]:
            if os.path.exists( i ): os.remove( i )

    def testReportAll ( self ):
        filename = "./testFiles/slha/gluino_squarks.slha"
        inifile = "testReportAll.ini"
        outputfile = runMain(filename, inifile=inifile, suppressStdout = True)
        smodelsOutput = importModule ( outputfile )
        from default_reportAll import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'database version']
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.02,
                            ignore=ignoreFields, fname=outputfile,
	                          fname2 = "default_reportAll.py" )
        self.assertTrue(equals)
        self.removeOutputs ( outputfile )

if __name__ == "__main__":
    unittest.main()
