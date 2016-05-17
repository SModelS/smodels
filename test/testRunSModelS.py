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
from smodels.installation import installDirectory
from smodels.tools import summaryReader
from runSModelS import main

class RunSModelSTest(unittest.TestCase):
    def runMain(self, filename ):
        suppressStdout = True
        if suppressStdout:
            a=sys.stdout
            sys.stdout = open ( "stdout.log", "w" )
        out = "%s/test/unitTestOutput" % installDirectory()
        #if os.path.exists ( out ): 
        #    os.unlink ( out )
        main(filename, 
             parameterFile="%s/test/testParameters.ini" %installDirectory(), 
             outputDir= out,
             verbosity = 'error' )
        outputfile = summaryReader.Summary(
                "%s/test/unitTestOutput/%s.smodels" % (installDirectory(),os.path.basename(filename)))
        if suppressStdout:
            sys.stdout = a
        return outputfile

    def testGoodFile(self):

        filename = "%s/inputFiles/slha/gluino_squarks.slha" % \
                    (installDirectory() )
        outputfile = self.runMain(filename )
        sample = summaryReader.Summary(
                "%s/test/gluino_squarks_default.txt" %installDirectory())
        #if not ( sample==outputfile ):
        #    print
        #    print "%s != %s" % ( os.path.basename(outputfile.filename), 
        #                         os.path.basename(sample.filename) )
        try:
            self.assertEquals(sample, outputfile )
        except AssertionError,e:
            msg = "%s != %s" % ( sample, outputfile )
            raise AssertionError ( msg )
        

    def testBadFile(self):

        filename = "%s/inputFiles/slha/I_dont_exist.slha" % \
                    (installDirectory() )
        outputfile = self.runMain (filename ) 
        sample = summaryReader.Summary(
                "%s/test/summary_bad_default.txt" %installDirectory())
        try:
            self.assertEquals(sample, outputfile )
        except AssertionError,e:
            msg = "%s != %s" % ( sample, outputfile )
            raise AssertionError ( msg )

if __name__ == "__main__":
    unittest.main()
