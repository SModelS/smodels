#!/usr/bin/env python3
 
"""
.. module:: testCombined
   :synopsis: Tests the combination code
 
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
 
"""
 
import sys,os,imp
sys.path.insert(0,"../")
import unittest
from os.path import join, basename
from smodels.installation import installDirectory as iDir
from databaseLoader import database ## to make sure the db exists
from unitTestHelpers import equalObjs
from smodels.tools.runSModelS import run
import redirector
import unum
 
from smodels.tools.smodelsLogging import logger, setLogLevel
 
class CombinedTest(unittest.TestCase):
    def runMain(self, filename, timeout = 0, suppressStdout=True, development=False,
                 inifile = "testParameters_agg.ini" ):
        to = None
        level = 'debug'
        if suppressStdout:
            level = 'error'
            to = os.devnull
        with redirector.stdout_redirected ( to = to ):
            out = join( iDir(), "test/unitTestOutput" )
            setLogLevel ( level )
            run(filename, parameterFile=join ( iDir(), "test/%s" % inifile ),
                 outputDir= out, db= database, timeout = timeout,
                 development = development)
            sfile = join(iDir(),"test/unitTestOutput/%s.py" % basename(filename))
            return sfile
 
    
    def testCombinedResult(self):
        filename = join ( iDir(), "inputFiles/slha/gluino_squarks.slha" )
        outputfile = self.runMain(filename)
        with open( outputfile, 'rb') as fp: ## imports file with dots in name
            output_module = imp.load_module("output",fp,outputfile, ('.py', 'rb', imp.PY_SOURCE) )
            smodelsOutput = output_module.smodelsOutput
            from gluino_squarks_default_agg import smodelsOutputDefault
            ignoreFields = ['input file','smodels version', 'ncpus', 'database version']
            smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                        key=lambda res: [res['theory prediction (fb)'],res['TxNames'],
                        res['AnalysisID'],res['DataSetID']])
            equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.02,
                               ignore=ignoreFields)
            self.assertTrue(equals)
        for i in [ outputfile, outputfile.replace(".py",".pyc") ]:
            if os.path.exists ( i ):
                os.remove ( i )
 
if __name__ == "__main__":
    unittest.main()
