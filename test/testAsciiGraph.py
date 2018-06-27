#!/usr/bin/env python

"""
.. module:: testAsciiGraph
   :synopsis: Tests the ascii grapher.
              Depends also on lheReader, lheDecomposer.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import unittest
import sys
sys.path.insert(0,"../")
from smodels.particleDefinitions import BSM
from smodels.theory.model import Model
from smodels.tools import asciiGraph
from smodels.installation import installDirectory
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
        
        filename = "%sinputFiles/lhe/simplyGluino.lhe" % (installDirectory() )        
        model = Model(filename, BSM)
        model.updateParticles()
        
        
        topList = decomposer.decompose(model)
        element = topList.getElements()[0]
        print(element)

        d1=self.orig().split("\n")
        d2=asciiGraph.asciidraw(element, border=True ).split("\n")
        print(d2)
#         self.assertEqual(d1,d2)



if __name__ == "__main__":
    unittest.main()
