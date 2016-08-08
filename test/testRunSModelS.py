#!/usr/bin/env python

"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""

import sys
sys.path.insert(0,"../")
import unittest
import os
import glob
from os.path import join, basename
from smodels.installation import installDirectory as iDir
from smodels.tools import summaryReader, crashReport
from smodels.tools.timeOut import NoTime
from runSModelS import main

class RunSModelSTest(unittest.TestCase):
    def runMain(self, filename, timeout = 0, suppressStdout=True, development=False ):
        if suppressStdout:
            a=sys.stdout
            sys.stdout = open ( "stdout.log", "w" )
        out = join ( iDir(), "test/unitTestOutput" )
        main(filename, parameterFile=join ( iDir(), "test/testParameters.ini" ),
             outputDir= out, verbosity = 'error', db= None, timeout = timeout,
             development = development )
        outputfile = None
        bname = basename ( filename )
        sfile = join ( iDir(), "test/unitTestOutput/%s.smodels" % bname )
        if os.path.exists ( sfile ):
            outputfile = summaryReader.Summary( sfile )
        if suppressStdout:
            sys.stdout = a
        return outputfile

    def testMultipleFiles ( self ):
        out = join ( iDir(), "test/unitTestOutput" )
        for i in os.listdir ( out ):
            if i[-8:]==".smodels":
                os.unlink ( os.path.join ( out, i ) )
        filename = join ( iDir(), "inputFiles/slha/" )
        self.runMain ( filename )
        nout = len( list ( glob.iglob ( "unitTestOutput/*smodels" )) )
        nin = len( list ( glob.iglob ( "%s/*slha" % filename )) )
        self.assertTrue ( nout == nin )

    def timeoutRun(self):
        filename = join ( iDir(), "inputFiles/slha/gluino_squarks.slha" )
        outputfile = self.runMain(filename, timeout=1, suppressStdout=True,
                     development=True )

    def testTimeout(self):
        self.assertRaises ( NoTime, self.timeoutRun )

    def testGoodFile(self):
        filename = join ( iDir(), "inputFiles/slha/gluino_squarks.slha" )
        outputfile = self.runMain(filename )
        sfile = join ( iDir(), "test/gluino_squarks_default.txt" )
        sample = summaryReader.Summary( sfile )
        try:
            self.assertEquals(sample, outputfile )
        except AssertionError,e:
            msg = "%s != %s" % ( sample, outputfile )
            raise AssertionError ( msg )
        
    def testBadFile(self):

        filename = join ( iDir(), "inputFiles/slha/I_dont_exist.slha" )
        outputfile = self.runMain (filename ) 
        sfile = join ( iDir(), "test/summary_bad_default.txt" )
        sample = summaryReader.Summary( sfile )
        try:
            self.assertEquals(sample, outputfile )
        except AssertionError,e:
            msg = "%s != %s" % ( sample, outputfile )
            raise AssertionError ( msg )

    def testCrash(self):
        filename = join ( iDir(), "inputFiles/slha/broken.slha" )
        outputfile = self.runMain (filename )
        ts = 0
        cf = None
        for f in os.listdir("."):
            if not ".crash" in f: continue
            nts = f.replace("smodels-","").replace(".crash","")
            if int(nts) > ts: cf = f
        inp, par = crashReport.readCrashReportFile(cf)
        print inp, par
        self.assertEquals(open(filename).readlines(), open(inp).readlines())
        self.assertEquals(open("testParameters.ini").readlines(), open(par).readlines())

if __name__ == "__main__":
    unittest.main()
