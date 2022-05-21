#!/usr/bin/env python3

"""
.. module:: testQNumbers
   :synopsis: Tests QNUMBERS

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
sys.path.insert(0,"../")
import unittest
from smodels.tools.smodelsLogging import setLogLevel
from smodels.tools import runtime
from smodels.theory.model import Model
from smodels import particlesLoader
from smodels.share.models.SMparticles import SMList
from imp import reload
setLogLevel('error')


class ModelsTest(unittest.TestCase):
    #Check if models loaded through python module or through SLHA QNUMBERS give the same model
    def testRuntimeImport(self):
        filename = "./testFiles/slha/gluino_squarks.slha"
        runtime.modelFile = 'mssm'
        reload(particlesLoader)
        from smodels.particlesLoader import BSMList
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(filename)

        runtime.modelFile = 'mssmQNumbers.slha'
        reload(particlesLoader)
        from smodels.particlesLoader import BSMList
        modelB = Model(BSMparticles=BSMList, SMparticles=SMList)
        modelB.updateParticles(filename)

        for ptc in model.BSMparticles:
            ptcB = modelB.getParticlesWith(pdg = ptc.pdg)
            if not ptcB: #If particule is its own anti-particle, it should not appear in modelB
                ptcB = modelB.getParticlesWith(pdg = -ptc.pdg)
            self.assertEqual(len(ptcB),1)
            ptcB = ptcB[0]
            for attr,val in ptc.__dict__.items():
                if attr in ['_id','label','_comp','pdg','_isInvisible']:
                    continue
                valB = getattr(ptcB,attr)
                if attr == 'decays':
                    self.assertEqual(len(val),len(valB))
                    continue
                self.assertEqual(val,valB)

if __name__ == "__main__":
    unittest.main()
