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
from smodels.tools import summaryReader
from runSModelS import main

class RunSModelSTest(unittest.TestCase):
    def runMain(self, filename ):
        suppressStdout = True
        if suppressStdout:
            a=sys.stdout
            sys.stdout = open ( "stdout.log", "w" )
        out = join ( iDir(), "test/unitTestOutput" )
        main(filename, parameterFile=join ( iDir(), "test/testParameters.ini" ),
             outputDir= out, verbosity = 'error' )
        sfile = join ( iDir(), 
                "test/unitTestOutput/%s.smodels" % basename ( filename ) )
        outputfile = summaryReader.Summary( sfile )
        if suppressStdout:
            sys.stdout = a
        return outputfile

    def testMultipleFiles ( self ):
        filename = join ( iDir(), "inputFiles/slha/" )
        suppressStdout = True
        if suppressStdout:
            a=sys.stdout
            sys.stdout = open ( "stdout.log", "w" )
        out = join ( iDir(), "test/unitTestOutput" )
        for i in os.listdir ( out ):
            if i[-8:]==".smodels":
                os.unlink ( os.path.join ( out, i ) )
        mypid = os.getpid()
        try:
            main(filename, parameterFile=join ( iDir(), "test/testParameters.ini" ),
                 outputDir= out, verbosity = 'error' )
        except SystemExit,e:
            # the daughters!
            return True
            # pass
        if mypid != os.getpid():
            return True
        if suppressStdout:
            sys.stdout = a
        nout = len( list ( glob.iglob ( "unitTestOutput/*smodels" )) )
        nin = len( list ( glob.iglob ( "%s/*slha" % filename )) )
        self.assertTrue ( nout == nin )

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

if __name__ == "__main__":
    unittest.main()
