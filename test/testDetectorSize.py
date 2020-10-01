#!/usr/bin/env python3

"""
.. module:: testReweighting
   :synopsis: Tests the function of lifetime reweighting
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.share.models import mssm
from smodels.theory.branch import Branch
from smodels.theory.element import Element
from smodels.tools.reweighting import calculateProbabilities,defaultEffReweight
from smodels.tools.physicsUnits import GeV
from databaseLoader import database

class DetectorSizeTest(unittest.TestCase):

    def testProbabilities(self):

        Leff_inner = 0.000769
        Leff_outer = 7.0
        widths = [1e-20,1e-17,1e-16,1e-14,1e-11]
        fracs = [[0.9996,3.897e-08,0.0003546],
                [0.7014,3.897e-05,0.2986],
                [0.0288,0.0003896,0.9708],
                [0.0,0.03822,0.9618],
                [0,1,0]]
        for i,w in enumerate(widths):
            prob = calculateProbabilities(w,Leff_inner,Leff_outer)
            F_long, F_prompt, F_displaced = prob['F_long'],prob['F_prompt'],prob['F_displaced']
            f = [F_long, F_prompt, F_displaced]
            for j,ff in enumerate((f)):
                diff = abs(ff-fracs[i][j])
                if diff > 0.0:
                    diff = diff/(ff+fracs[i][j])
            self.assertLess(diff,1e-3)

        Leff_inner = 0.1
        Leff_outer = 15.0
        fracs = [[0.9992,5.068e-06,0.0007548],
                 [0.4676,0.005055,0.5274],
                 [0.0004997,0.04941,0.9501],
                 [0,0.9937,0.006297],
                 [0,1,0]]

        for i,w in enumerate(widths):
            prob = calculateProbabilities(w,Leff_inner,Leff_outer)
            F_long, F_prompt, F_displaced = prob['F_long'],prob['F_prompt'],prob['F_displaced']
            f = [F_long, F_prompt, F_displaced]
            for j,ff in enumerate((f)):
                diff = abs(ff-fracs[i][j])
                if diff > 0.0:
                    diff = diff/(ff+fracs[i][j])
            self.assertLess(diff,1e-3)

    def testEffReweight(self):

        n1 = mssm.n1.copy()
        n1.totalwidth = 1e-18*GeV
        st1 = mssm.st1.copy()
        st1.totalwidth = 5e-15*GeV
        gluino = mssm.gluino.copy()
        gluino.totalwidth = 1e-15*GeV
        branch1 = Branch()
        branch1.oddParticles = [st1,n1]
        branch2 = Branch()
        branch2.oddParticles = [gluino,n1]
        el1 = Element([branch1,branch2])
        r = defaultEffReweight(el1) #Use default values
        self.assertAlmostEqual(r*1e5,6.991,places=2)

        r = defaultEffReweight(el1,Leff_inner=0.1,Leff_outer=15.0)
        self.assertAlmostEqual(r,0.314,places=2)


    def testTxnameDataReweight(self):

        c1 = mssm.c1.copy()
        c1.totalwidth = 1e-17*GeV
        c1.mass = 100*GeV
        gluino = mssm.gluino.copy()
        gluino.totalwidth = 1e-15*GeV
        gluino.mass = 500*GeV
        branch1 = Branch()
        branch1.oddParticles = [gluino,c1]
        branch2 = Branch()
        branch2.oddParticles = [gluino,c1]
        el1 = Element([branch1,branch2])

        listofanalyses = database.getExpResults(
                analysisIDs= [ "CMS-EXO-13-006"], dataTypes = ['efficiencyMap'])

        exp = listofanalyses[0]
        ds = [d for d in exp.datasets if d.dataInfo.dataId == 'c000'][0]
        tx = [t for t in ds.txnameList if str(t) == 'THSCPM3'][0]
        effDefault = tx.getEfficiencyFor(el1)
        self.assertAlmostEqual(tx.txnameData.Leff_inner,0.007)
        self.assertAlmostEqual(tx.txnameData.Leff_outer,7.0)
        self.assertAlmostEqual(effDefault/1e-5,5.223348,places=3)

        ds = [d for d in exp.datasets if d.dataInfo.dataId == 'c000track'][0]
        tx = [t for t in ds.txnameList if str(t) == 'THSCPM3'][0]
        effNewSize = tx.getEfficiencyFor(el1)
        self.assertAlmostEqual(tx.txnameData.Leff_inner,0.007)
        self.assertAlmostEqual(tx.txnameData.Leff_outer,3.0)
        self.assertAlmostEqual(effNewSize/1e-5,7.83466,places=3)


        branch1 = Branch()
        branch1.oddParticles = [c1]
        branch2 = Branch()
        branch2.oddParticles = [c1]
        el2 = Element([branch1,branch2])

        ds = [d for d in exp.datasets if d.dataInfo.dataId == 'c000'][0]
        tx = [t for t in ds.txnameList if str(t) == 'THSCPM1'][0]
        effDefault = tx.getEfficiencyFor(el2)
        self.assertAlmostEqual(tx.txnameData.Leff_inner,0.007)
        self.assertAlmostEqual(tx.txnameData.Leff_outer,7.0)
        self.assertAlmostEqual(effDefault,0.073292924,places=3)

        ds = [d for d in exp.datasets if d.dataInfo.dataId == 'c000track'][0]
        tx = [t for t in ds.txnameList if str(t) == 'THSCPM1'][0]
        effNewSize = tx.getEfficiencyFor(el2)
        self.assertAlmostEqual(tx.txnameData.Leff_inner,0.1)
        self.assertAlmostEqual(tx.txnameData.Leff_outer,5.0)
        self.assertAlmostEqual(effNewSize,0.0897630,places=3)


if __name__ == "__main__":
    unittest.main()
