#!/usr/bin/env python3

"""
.. module:: testCombined
   :synopsis: Tests the combination code

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
import os
sys.path.insert(0, "../")
import unittest
from unitTestHelpers import equalObjs, runMain, importModule
from smodels.base.smodelsLogging import setLogLevel

class CombinedTest(unittest.TestCase):

    def testCombinedResult(self):
        filename = "./testFiles/slha/gluino_squarks.slha"
        from smodels.base.smodelsLogging import logger, setLogLevel
        setLogLevel ( "fatal" )
        outputfile = runMain(filename, inifile="testParameters_agg.ini", suppressStdout=True)
        smodelsOutput = importModule(outputfile)
        from gluino_squarks_default_agg import smodelsOutputDefault
        ignoreFields = ['input file', 'smodels version', 'ncpus', 'database version',
                        'model', 'promptwidth', 'stablewidth', 'checkinput',
                        'doinvisible', 'docompress', 'computestatistics']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                                                 key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput, smodelsOutputDefault, allowedRelDiff=0.02,
                           ignore=ignoreFields, fname=outputfile)
        if equals != True:
            logger.error( f"gluino_squarks_default_agg.py differs from {outputfile}!" )
        self.assertTrue(equals)
        for i in [outputfile, outputfile.replace(".py", ".pyc")]:
            if os.path.exists(i):
                os.remove(i)


if __name__ == "__main__":
    setLogLevel("debug")
    unittest.main()
