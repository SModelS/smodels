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
