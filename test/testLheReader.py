#!/usr/bin/env python333

"""
.. module:: testLheReader
   :synopsis: Tests the lheReader
              Depends also on lheDecomposer.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
sys.path.insert(0,"../")
import unittest

class LheReaderTest(unittest.TestCase):
    def testReader(self):
        """ test the LheReader """
        from smodels.theory import lheReader, lheDecomposer, crossSection
        from smodels.installation import installDirectory
        from smodels.tools.physicsUnits import GeV

        filename = "%sinputFiles/lhe/simplyGluino.lhe" % (installDirectory() )
        reader = lheReader.LheReader(filename)
        event = reader.next()
        element = lheDecomposer.elementFromEvent(event,
                                             crossSection.XSectionList())
        s=str(element)
        assert ( s == "[[[q,q]],[[q,q]]]" )
        b0=element.branches[0]
        sb0=str(b0)
        assert ( sb0 == "[[q,q]]" )
        assert ( b0.masses[0]-675*GeV ) < .1*GeV
        assert ( b0.masses[1]-600*GeV ) < .1*GeV

if __name__ == "__main__":
    unittest.main()
