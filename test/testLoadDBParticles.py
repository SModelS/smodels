#!/usr/bin/env python3

"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
sys.path.insert(0,"../")
import unittest
from smodels.experiment.databaseObj import Database
from databaseLoader import database

class LoadDBParticlesTest(unittest.TestCase):

    def testLoadDBPaticles(self):
        modelA = database.databaseParticles
        self.assertEqual(modelA.label,'databaseParticles.py') #Simple test to see if databaseParticles.py is being used

    def testFallBack(self):
        # from smodels.installation import version
        # dbpath = "https://smodels.github.io/database/unittest%s" % version().replace(".","")
        # dbpath = "unittest"
        dbpath = "https://smodels.github.io/database/unittest200rc9"
        dbOld = Database( dbpath, discard_zeroes = False, force_load ='pcl')
        model = dbOld.databaseParticles
        self.assertEqual(model.label,'DB Final States (default)') #Simple fallback test


if __name__ == "__main__":
    unittest.main()
