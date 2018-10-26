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
from smodels.installation import installDirectory
from smodels.theory import decomposer

class AsciiTest(unittest.TestCase):
    def orig(self):
        return """ /------------\\
 |    d  d    |
 |    \ /     |
 | ----*----  |
 | ----*----  |
 |    / \     |
 |    d  d    |
 \------------/
"""

    def testGraph(self):
        """ draw ascii graph """
        
        filename = "%sinputFiles/lhe/simplyGluino.lhe" % (installDirectory() )        
        model = Model(inputFile=filename, BSMparticles = BSMList, SMparticles = SMList)
        model.updateParticles()
        
        
        topList = decomposer.decompose(model)
        element = topList.getElements()[0]

        d1=self.orig().split("\n")
        d2=asciiGraph.asciidraw(element, border=True ).split("\n")
        self.assertEqual(d1,d2)



if __name__ == "__main__":
    unittest.main()
