#!/usr/bin/env python

"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""
import unittest
import os
import sys
sys.path.append("../")
from smodels.installation import installDirectory
from smodels.tools import summaryReader
from runSModelS import main

class RunSModelSTest(unittest.TestCase):
    def testGoodFile(self):

        filename = "%s/inputFiles/slha/gluino_squarks.slha" % (installDirectory() )
#        print "filename=",filename
        main(filename, 
             parameterFile="%s/test/testParameters.ini" %installDirectory(), 
             outputFile="%s/test/unitTestOutput.txt" %installDirectory())
        sample = summaryReader.Summary(
                "%s/test/summary_default.txt" %installDirectory())
        outputfile = summaryReader.Summary(
                "%s/test/unitTestOutput.txt" %installDirectory())
        self.assertEquals(sample, outputfile)

    def testBadFile(self):

        filename = "%s/inputFiles/slha/I_dont_exist.slha" % (installDirectory() )
#        print "filename=",filename
        main(filename, 
             parameterFile="%s/test/testParameters.ini" %installDirectory(), 
             outputFile="%s/test/unitTestOutput.txt" %installDirectory())
        sample = summaryReader.Summary(
                "%s/test/summary_bad_default.txt" %installDirectory())
        outputfile = summaryReader.Summary(
                "%s/test/unitTestOutput.txt" %installDirectory())
        self.assertEquals(sample, outputfile)

if __name__ == "__main__":
    unittest.main()
