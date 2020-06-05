#!/usr/bin/env python3

"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys,os,imp
sys.path.insert(0,"../")
import unittest
import glob
from os.path import join
from smodels.installation import installDirectory as iDir
from smodels.tools import crashReport
from smodels.tools.timeOut import NoTime
from unitTestHelpers import equalObjs, runMain
import time

from smodels.tools.smodelsLogging import logger

class RunSModelSTest(unittest.TestCase):
    def testMultipleFiles( self ):
        out = join( iDir(), "test/unitTestOutput")
        for i in os.listdir( out ):
            if i[-8:]==".smodels":
                os.unlink( os.path.join ( out, i ))
        dirname = "./testFiles/slha/"
        runMain(dirname)
        nout = len([i for i in glob.iglob("unitTestOutput/*smodels") if not "~" in i])
        nin = len([i for i in glob.iglob("%s/*slha" % dirname) if not "~" in i])
        if nout != nin:
            logger.error("Number of output file (%d) differ from number of input files (%d)" % (nout, nin))
        self.assertTrue( nout == nin )

    def timeoutRun(self):
        filename = "./testFiles/slha/complicated.slha"
        outputfile = runMain(filename, timeout=1, suppressStdout=True,
                             development=True, inifile = "timeout.ini" )

    def testTimeout(self):
        self.assertRaises(NoTime, self.timeoutRun)

    def removeOutputs ( self, f ):
        """ remove cruft outputfiles """
        for i in [ f, f.replace(".py",".pyc") ]:
            if os.path.exists ( i ): os.remove ( i )

    def testGoodFile(self):
        filename = "./testFiles/slha/gluino_squarks.slha"
        outputfile = runMain(filename)
        with open( outputfile, 'rb') as fp: ## imports file with dots in name
            output_module = imp.load_module("output",fp,outputfile, ('.py', 'rb', imp.PY_SOURCE) )
            smodelsOutput = output_module.smodelsOutput
        from gluino_squarks_default import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'database version']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.02,
                           ignore=ignoreFields)
        self.assertTrue(equals)
        self.removeOutputs ( outputfile )

    def testPyhfCombination(self):
        filename = "./testFiles/slha/T6bbHH_pyhf.slha"
        inifile = "./testParameters_pyhf.ini"
        outputfile = runMain(filename, inifile=inifile)
        smodelsOutput = importModule(outputfile)
        from pyhf_default import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'database version']
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.02,
                           ignore=ignoreFields, fname=outputfile )
        self.assertTrue(equals)

    def testGoodFile13(self):

        filename = "./testFiles/slha/simplyGluino.slha"
        outputfile = runMain(filename,suppressStdout = True )
        with open( outputfile, 'rb') as fp: ## imports file with dots in name
            output_module = imp.load_module("output",fp,outputfile, ('.py', 'rb', imp.PY_SOURCE) )
            smodelsOutput = output_module.smodelsOutput
        from simplyGluino_default import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'Element', 'database version' ]
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.02,
                           ignore=ignoreFields)
        if not equals:
            e =  "output13.py and simplyGluino_default.py differ!"
            logger.error ( e )
            # raise AssertionError ( e )

        self.assertTrue(equals)

        ## test went through, so remove the output files
        self.removeOutputs ( outputfile )

    def testBadFile(self):
        # since 112 we skip non-existing slha files!
        filename = "./testFiles/slha/I_dont_exist.slha"
        of="unitTestOutput/I_dont_exist.slha.py"
        self.removeOutputs ( of )
        outputfile = runMain(filename  )
        self.assertTrue ( of in outputfile )
        self.assertTrue ( not os.path.exists ( outputfile ) )

    def cleanUp ( self ):
        for f in os.listdir("."):
            if ".crash" in f: os.remove(f)
        for i in [ "crash_report_parameter", "crash_report_input" ]:
            if os.path.exists ( i ): os.remove ( i )

    def testCrash(self):
        filename = "./testFiles/slha/gluino_squarks.slha"
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
