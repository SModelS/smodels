#!/usr/bin/env python3

"""
.. module:: testIntermediateBSM
   :synopsis: Test the implementation of specific quantum numbers for intermediate BSM states
              in the database.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
sys.path.insert(0,"../")
from smodels.theory.element import Element
from smodels.experiment.databaseParticles import MET,HSCP,anyOdd,gluino,chargino
from smodels.share.models import mssm
from smodels.experiment import infoObj
from smodels.experiment.txnameObj import TxName, rescaleWidth, unscaleWidth
from smodels.tools.physicsUnits import GeV
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.theory.auxiliaryFunctions import flattenArray
from databaseLoader import database
import numpy as np
import unittest


class IntermediateTest(unittest.TestCase):

    def testElement(self):
        """Test the creation of elements using intermediate states."""
        #Without intermediate or final states:
        el = Element(info="[[[q,q]],[[q,q]]]")
        self.assertEqual(el.oddParticles,[[anyOdd,MET],[anyOdd,MET]])

        #Only with final states:
        el = Element(info="[[[q,q]],[[q,q]]]", finalState = ['MET','HSCP'])
        self.assertEqual(el.oddParticles,[[anyOdd,MET],[anyOdd,HSCP]])

        #Only with intermediate states:
        el = Element(info="[[[q,q]],[[q,q]]]", intermediateState = [['gluino'],['C1+']])
        self.assertEqual(el.oddParticles,[[gluino,MET],[chargino,MET]])

        #With intermediate and final states:
        el = Element(info="[[[q,q]],[[q,q]]]", finalState = ['HSCP','HSCP'],
                        intermediateState = [['gluino'],['C1+']])
        self.assertEqual(el.oddParticles,[[gluino,HSCP],[chargino,HSCP]])

        #Slightly more complicated element:
        el = Element(info="[[[q],[q]],[[q,q]]]", finalState = ['HSCP','HSCP'],
                        intermediateState = [['gluino','C1+'],['C1+']])
        self.assertEqual(el.oddParticles,[[gluino,chargino,HSCP],[chargino,HSCP]])

    def testDatabase(self):
        """Test the database without intermediate states."""
        #Create test elements
        c1 = mssm.c1
        n1 = mssm.n1
        g = mssm.gluino
        n1.mass = 100*GeV
        g.mass = 1000*GeV
        c1.mass = 500*GeV
        el1 = Element(info="[[[q,q]],[[q,q]]]")
        el2 = Element(info="[[[q,q]],[[q,q]]]")
        el1.branches[0].oddParticles = [g,n1]
        el1.branches[1].oddParticles = [g,n1]
        el2.branches[0].oddParticles = [c1,n1]
        el2.branches[1].oddParticles = [c1,n1]

        #Test result without intermediate states:
        expRes = database.getExpResults(analysisIDs=["ATLAS-SUSY-2016-08"],
                    datasetIDs=[None], txnames=["T5Disp" ] )
        txname=expRes[0].datasets[0].txnameList[0]

        newEl = txname.hasElementAs(el1)
        self.assertEqual(newEl,el1)
        newEl = txname.hasElementAs(el2)
        self.assertEqual(newEl,el2)

    def testDatabaseInt(self):
        """Test the database with intermediate states."""
        #Create test elements
        c1 = mssm.c1
        n1 = mssm.n1
        g = mssm.gluino
        n1.mass = 100*GeV
        g.mass = 1000*GeV
        c1.mass = 500*GeV
        el1 = Element(info="[[[q,q]],[[q,q]]]")
        el2 = Element(info="[[[q,q]],[[q,q]]]")
        el1.branches[0].oddParticles = [g,n1]
        el1.branches[1].oddParticles = [g,n1]
        el2.branches[0].oddParticles = [c1,n1]
        el2.branches[1].oddParticles = [c1,n1]

        #Test result without intermediate states:
        expRes = database.getExpResults(analysisIDs=["ATLAS-SUSY-2016-081"],
                    datasetIDs=[None], txnames=["T5Disp" ] )
        txname=expRes[0].datasets[0].txnameList[0]

        #el1 should match, since txname requires gluinos
        newEl = txname.hasElementAs(el1)
        self.assertEqual(newEl,el1)
        self.assertNotEqual(newEl,el2)
        #el2 should not match, since txname requires gluinos
        newEl = txname.hasElementAs(el2)
        self.assertEqual(False,newEl)



if __name__ == "__main__":
    unittest.main()
