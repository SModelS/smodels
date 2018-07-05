#!/usr/bin/env python

"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys,shutil,os,imp
sys.path.insert(0,"../")
import unittest
import glob
from os.path import join
from smodels.installation import installDirectory as iDir
from smodels.tools import crashReport
from smodels.tools.timeOut import NoTime
from unitTestHelpers import equalObjs, runMain
import time
from smodels.tools.smodelsLogging import logger, setLogLevel



class RunSModelSTest(unittest.TestCase):

    def testMultipleFiles( self ):
        out = join( iDir(), "test/unitTestOutput")
        for i in os.listdir( out ):
            if i[-8:]==".smodels":
                os.unlink( os.path.join ( out, i ))
        dirname = join( iDir(), "inputFiles/slha/" )
        runMain(dirname)
        nout = len([i for i in glob.iglob("unitTestOutput/*smodels") if not "~" in i])
        nin = len([i for i in glob.iglob("%s/*slha" % dirname) if not "~" in i])
        if nout != nin:
            logger.error("Number of output file (%d) differ from number of input files (%d)" % (nout, nin))
        self.assertTrue( nout == nin )
     
    def timeoutRun(self):
        filename = join ( iDir(), "inputFiles/slha/complicated.slha" )
        outputfile = runMain(filename, timeout=1, suppressStdout=True,
                             development=True, inifile = "timeout.ini" )
   
    def testTimeout(self):
        self.assertRaises(NoTime, self.timeoutRun)
 
    def testGoodFile(self):
        filename = join ( iDir(), "inputFiles/slha/gluino_squarks.slha" )
        outputfile = runMain(filename)
        with open( outputfile, 'rb') as fp: ## imports file with dots in name
            output_module = imp.load_module("output",fp,outputfile, ('.py', 'rb', imp.PY_SOURCE) )
            smodelsOutput = output_module.smodelsOutput

        from gluino_squarks_default import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: [res['theory prediction (fb)'],res['TxNames'],
                    res['AnalysisID'],res['DataSetID']])
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.02,
                           ignore=ignoreFields)
        for i in [ './output.py', './output.pyc' ]:
            if os.path.exists ( i ): os.remove ( i )
        self.assertTrue(equals)
    
    def testGoodFile13(self):
          
        filename = join ( iDir(), "inputFiles/slha/simplyGluino.slha" )
        outputfile = runMain(filename)
        with open( outputfile, 'rb') as fp: ## imports file with dots in name
            output_module = imp.load_module("output",fp,outputfile, ('.py', 'rb', imp.PY_SOURCE) )
            smodelsOutput = output_module.smodelsOutput
        
        from simplyGluino_default import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'Element']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: [res['theory prediction (fb)'],res['TxNames'],
                    res['AnalysisID'],res['DataSetID']])
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.07,
                           ignore=ignoreFields)
        if not equals:
            print "output13.py and simplyGluino_default.py differ!"
        for i in [ './output13.py', './output13.pyc' ]:
            if os.path.exists ( i ):
                continue
                os.remove ( i )
        self.assertTrue(equals)
        
    def testGoodFileHSCP(self):
        filename = join ( iDir(), "inputFiles/slha/longLived.slha" )
        outputfile = runMain(filename)
        with open( outputfile, 'rb') as fp: ## imports file with dots in name
            output_module = imp.load_module("output",fp,outputfile, ('.py', 'rb', imp.PY_SOURCE) )
            smodelsOutput = output_module.smodelsOutput

        from longLived_default import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: [res['theory prediction (fb)'],res['TxNames'],
                    res['AnalysisID'],res['DataSetID']])
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.02,
                           ignore=ignoreFields)
            
        for i in [ './outputHSC.py', './outputHSCP.pyc' ]:
            if os.path.exists ( i ): os.remove ( i )
        self.assertTrue(equals)               
  
    def testBadFile(self):
        filename = join (iDir(), "inputFiles/slha/I_dont_exist.slha" )
        outputfile = runMain(filename)
        with open( outputfile, 'rb') as fp: ## imports file with dots in name
            output_module = imp.load_module("output",fp,outputfile, ('.py', 'rb', imp.PY_SOURCE) )
            smodelsOutput = output_module.smodelsOutput        
        
        from bad_default import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus']
        equals = equalObjs( smodelsOutput,smodelsOutputDefault,allowedDiff=0.,
                            ignore=ignoreFields)
        for i in [ "./bad_output.py", "./bad_output.pyc" ]:
            if os.path.exists ( i ):
                os.remove( i )
        self.assertTrue( equals )
  
    def cleanUp ( self ):
        for f in os.listdir("."):
            if ".crash" in f: os.remove(f)
        for i in [ "crash_report_parameter", "crash_report_input" ]:
            if os.path.exists ( i ): os.remove ( i )
  
    def testCrash(self):
        filename = join ( iDir(), "inputFiles/slha/gluino_squarks.slha" )
        ctr=0
        crash_file = None
        self.cleanUp()
        outputfile = runMain(filename, timeout=1, suppressStdout=True,
                                   inifile= "timeout.ini" )
        try:
            ## trying to sync!!
            import ctypes
            libc = ctypes.CDLL("libc.so.6")
            libc.sync()
        except (OSError,AttributeError,ImportError) as e:
            pass
        time.sleep(.1)
        for f in os.listdir("."):
            if ".crash" in f:
                crash_file = f
                ctr+=1
        self.assertEqual ( ctr, 1 )
        inp, par = crashReport.readCrashReportFile(crash_file)
   
        with open(filename) as f:
            with open(inp) as g:
                self.assertEqual(f.readlines(), g.readlines())
        with open("timeout.ini") as f:
            with open(par) as g:
                self.assertEqual( f.readlines(), g.readlines())
        self.cleanUp()

if __name__ == "__main__":
    unittest.main()
