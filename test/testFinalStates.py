#!/usr/bin/env python3

"""
.. module:: testFinalStates
   :synopsis: Checks if the finalStates (used in the database) are properly defined
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""
import sys
sys.path.insert(0,"../")
import unittest

from smodels.theory.model import Model
from smodels.share.models.SMparticles import SMList
from smodels.experiment.defaultFinalStates import finalStates


class FinalStateTest(unittest.TestCase):

    def testUniqueLabels(self):

        allLabels = finalStates.getValuesFor('label')
        #Check if for each label there is a unique particle object defined
        for label in allLabels:
            p = finalStates.getParticlesWith(label=label)
            self.assertTrue(len(p) == 1)

    def testUniquePDGs(self):

        model = Model(SMparticles=SMList,BSMparticles=[])
        allPDGs = model.getValuesFor('pdg')
        #Check if for each SM PDG there is a unique particle object defined
        for pdg in allPDGs:
            p = model.getParticlesWith(pdg=pdg)
            self.assertTrue(len(p) == 1)


if __name__ == "__main__":
    unittest.main()
