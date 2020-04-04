#!/usr/bin/env python3

"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys,os
sys.path.insert(0,"../")
import unittest
from unitTestHelpers import equalObjs, runMain, importModule
from smodels.tools.smodelsLogging import setLogLevel
from smodels.tools import runtime
from smodels import particlesLoader
from imp import reload
import subprocess
from smodels.tools.smodelsLogging import logger
from smodels.experiment.databaseObj import Database
from databaseLoader import database

class LoadDBParticlesTest(unittest.TestCase):

    def testLoadDBPaticles(self):
        modelA = database.databaseParticles
        modelB = database.expResultList[0].globalInfo._databaseParticles
        self.assertEqual(modelA.label,modelB.label)
        self.assertEqual(modelA.label,'databaseParticles.py')

    def testFallBack(self):
        dbOld = Database( "https://smodels.github.io/database/unittest200rc9",
                            discard_zeroes = False, force_load ='pcl')
        model = dbOld.databaseParticles
        self.assertEqual(model.label,'DB Final States (default)')


if __name__ == "__main__":
    unittest.main()
