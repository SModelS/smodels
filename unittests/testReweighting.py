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
from smodels.experiment.reweighting import calculateProbabilities,reweightFactorFor
from smodels.base.physicsUnits import GeV
from smodels.decomposition.theorySMS import TheorySMS
from smodels.experiment.expSMS import ExpSMS
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.particle import Particle
from smodels.base.model import Model
from collections import OrderedDict
from smodels.base.physicsUnits import fb, GeV, TeV
from smodels.base.crossSection import XSection,XSectionInfo,XSectionList
from smodels.share.models import SMparticles
from smodels.experiment.defaultFinalStates import finalStates
from unitTestHelpers import theorySMSFromString as fromString

slhafile = './testFiles/slha/lightEWinos.slha'

model = Model( BSMparticles=BSMList, SMparticles=SMList)
model.updateParticles(inputFile=slhafile,ignorePromptQNumbers=['spin'])
invisible = Particle(label='invisible',pdg=50000,mass=500*GeV,isSM=False)
model.BSMparticles.append(invisible)


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

        n1 = model.getParticle(label='N1')
        n1.totalwidth = 0.*GeV
        st1 = model.getParticle(label='st_1')
        st1.totalwidth = 1e-13*GeV
        gluino = model.getParticle(label='gluino')
        gluino.totalwidth = 1.*10**(-30)*GeV


        sms1 = fromString('(PV > N1,gluino)',model=model)
        f = reweightFactorFor(sms1, 'prompt')
        self.assertAlmostEqual(f,1.,places=3)
        f = reweightFactorFor(sms1, 'displaced')
        self.assertAlmostEqual(f,0.,places=3)

        sms2 = fromString('(PV > N1,st_1(1)), (st_1(1) > t+,N1)',model=model)
        F_prompt = 0.3228249017964917
        Fdisp = 0.6771750982035083

        f = reweightFactorFor(sms2, resType='prompt')
        self.assertAlmostEqual(f,F_prompt,places=3)
        f = reweightFactorFor(sms2, resType='displaced')
        self.assertAlmostEqual(f,Fdisp,places=3)


if __name__ == "__main__":
    unittest.main()
