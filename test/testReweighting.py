#!/usr/bin/env python3

"""
.. module:: testReweighting
   :synopsis: Tests the function of lifetime reweighting
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.share.models import SMparticles, mssm
from smodels.theory.branch import Branch
from smodels.theory.element import Element
from smodels.tools.reweighting import calculateProbabilities,reweightFactorFor
from smodels.tools.physicsUnits import GeV

class ReweightingTest(unittest.TestCase):

    def testcalculateProbabilities(self):

        gluino = mssm.gluino.copy()
        gluino.totalwidth = 1.*10**(-30)*GeV
        prob = calculateProbabilities(gluino.totalwidth.asNumber(GeV),
                                        Leff_inner=0.000769,Leff_outer=7.0)
        F_long, F_prompt, F_displaced = prob['F_long'],prob['F_prompt'],prob['F_displaced']
        self.assertAlmostEqual(F_long, 1.)
        self.assertEqual(F_prompt, 0.)
        self.assertAlmostEqual(F_displaced, 0.)

    def testreweightFactorFor(self):

        n1 = mssm.n1.copy()
        n1.totalwidth = 0.*GeV
        st1 = mssm.st1.copy()
        st1.totalwidth = 1e-13*GeV
        F_prompt = 0.3228249017964917
        Fdisp = 0.6771750982035083
        gluino = mssm.gluino.copy()
        gluino.totalwidth = 1.*10**(-30)*GeV
        t = SMparticles.t

        branch1 = Branch()
        branch1.oddParticles = [n1]
        branch2 = Branch()
        branch2.oddParticles = [gluino]
        el1 = Element([branch1,branch2])
        f = reweightFactorFor(el1, 'prompt')
        self.assertAlmostEqual(f,1.,places=3)
        f = reweightFactorFor(el1, 'displaced')
        self.assertAlmostEqual(f,0.,places=3)

        branch3 = Branch()
        branch3.oddParticles = [st1,n1]
        branch3.evenParticles = [[t]]
        el2 = Element([branch1,branch3])
        f = reweightFactorFor(el2, resType='prompt')
        self.assertAlmostEqual(f,F_prompt,places=3)
        f = reweightFactorFor(el2, resType='displaced')
        self.assertAlmostEqual(f,Fdisp,places=3)


if __name__ == "__main__":
    unittest.main()
