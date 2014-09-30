#!/usr/bin/env python

"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""
import unittest
import os
from smodels.installation import installDirectory
from runSModelS import main


class RunSModelSTest(unittest.TestCase):
    def testGoodFile(self):

        filename = "%sinputFiles/slha/andrePT4.slha" % (installDirectory() )
        main(filename, outputfile="%s/test/unitTestOutput.txt" %installDirectory())
        sample = open("%s/test/summary_default.txt" %installDirectory(),"r")
        outputfile = open("%s/test/unitTestOutput.txt" %installDirectory(), "r")
        self.assertEquals(outputfile.readlines(), sample.readlines())

if __name__ == "__main__":
    unittest.main()
