#!/usr/bin/env python3
 
"""
.. module:: testInteractivePlots
   :synopsis: Tests the interactive plot tool
 
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
 
"""
 
import sys,os,shutil
sys.path.insert(0,"../")
import unittest
from smodels.base.physicsUnits import GeV, fb, pb
from smodels.tools import databaseBrowser


class RunBrowserTest(unittest.TestCase):
    
    # this test corresponds to calling
    # ../smodelsTools.py interactive-plots -f ./testFiles/scanExample/smodels-output.tar.gz -s testFiles/scanExample/slhas.tar.gz -p iplots_parameters.py
    def testBrowser(self):

        browser = databaseBrowser.Browser("unittest")
        self.assertEqual(len(browser),16)

        # Check the upper limit for the HM200 signal region:
        ul = browser.getULForSR(expid='ATLAS-SUSY-2018-31', datasetID='SRA_H')
        self.assertAlmostEqual(ul.asNumber(fb),0.0368)

        # Check the upper limit for the LM300 signal region:
        ul = browser.getULForSR(expid='ATLAS-SUSY-2018-31', datasetID='SRC_26')
        self.assertAlmostEqual(ul.asNumber(fb),0.0455)

        # Check the upper limit for TChiWH simplified model:
        ul = browser.getULFor(expid='CMS-SUS-16-039', txname='TChiWH',
                              massarray = [[275*GeV,50*GeV]]*2)
        self.assertAlmostEqual(ul.asNumber(pb),1.12,2)

        # Check the upper limit for TChiWH simplified model:
        ul = browser.getULFor(expid='CMS-SUS-16-039', txname='TChiWH',
                              massarray = [[324*GeV,175*GeV]]*2)
        self.assertAlmostEqual(ul.asNumber(pb),1.75,1)


        # No results will be given if masses outside the data grid are used:
        ul = browser.getULFor(expid='CMS-SUS-16-039', txname='TChiWH',
                            massarray = [[400.*GeV,100.*GeV]]*2)
        self.assertTrue(ul is None)


if __name__ == "__main__":
    unittest.main()
