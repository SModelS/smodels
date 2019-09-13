#!/usr/bin/env python3

"""
.. module:: testAsciiGraph
   :synopsis: Tests the ascii grapher.
              Depends also on lheReader, lheDecomposer.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import unittest
import sys
sys.path.insert(0,"../")
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model
from smodels.tools import asciiGraph
from smodels.theory import decomposer

class AsciiTest(unittest.TestCase):
    def orig(self):
        return """ /------------\\
 |    q  q    |
 |    \ /     |
 | ----*----  |
 | ----*----  |
 |    / \     |
 |    q  q    |
 \------------/
"""

    def testGraph(self):
        """ draw ascii graph """
        
        filename = "./testFiles/lhe/simplyGluino.lhe"
        model = Model(BSMparticles = BSMList, SMparticles = SMList)
        model.updateParticles(filename)
        
        
        topList = decomposer.decompose(model, sigmacut=0)
        element = topList.getElements()[0]


        d1=self.orig().split("\n")
        d2=asciiGraph.asciidraw(element, border=True ).split("\n")
        self.assertEqual(d1,d2)



if __name__ == "__main__":
    unittest.main()
