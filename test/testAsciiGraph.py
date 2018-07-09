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
        from smodels.tools import asciiGraph
        from smodels.theory import lheReader, lheDecomposer, crossSection
        from smodels.installation import installDirectory

        filename = "%sinputFiles/lhe/simplyGluino.lhe" % (installDirectory() )
        reader = lheReader.LheReader(filename)
        event = reader.next()
        element = lheDecomposer.elementFromEvent(event,
                                             crossSection.XSectionList())
        d1=self.orig().split("\n")
        d2=asciiGraph.asciidraw ( element, border=True ).split("\n")
        #for (idx,line) in enumerate(d1):
        #        print "%d >>%s<< >>%s<<" % (idx,line, d2[idx] )
        #print d1==d2
        reader.close()
        self.assertEqual(d1,d2)



if __name__ == "__main__":
    unittest.main()
