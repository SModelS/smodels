#!/usr/bin/env python

"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""
import unittest
import os
from smodels.installation import installDirectory
from smodels.tools import summaryReader
from runSModelS import main


class RunSModelSTest(unittest.TestCase):
    def testGoodFile(self):

        filename = "%sinputFiles/slha/T1_7TeV.slha" % (installDirectory() )
        main(filename, parameterfile="%s/test/testParameters.ini" %installDirectory(), outputfile="%s/test/unitTestOutput.txt" %installDirectory())
        sample = summaryReader.Summary("%s/test/summary_default.txt" %installDirectory())
        outputfile = summaryReader.Summary("%s/test/unitTestOutput.txt" %installDirectory())
        self.assertEquals(sample, outputfile)

if __name__ == "__main__":
    unittest.main()
