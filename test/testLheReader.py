#!/usr/bin/env python

"""
.. module:: testLheReader
   :synopsis: Tests the lheReader
              Depends also on lheDecomposer.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import unittest

class LheReaderTest(unittest.TestCase):
    def testReader(self):
        """ test the LheReader """
        from smodels.theory import lheReader, lheDecomposer, crossSection
        from smodels.installation import installDirectory
        from smodels.tools.physicsUnits import GeV

        filename = "%sinputFiles/lhe/T1_1.lhe" % (installDirectory() )
        reader = lheReader.LheReader(filename)
        event = reader.next()
        element = lheDecomposer.elementFromEvent(event,
                                             crossSection.XSectionList())
        s=str(element)
        assert ( s == "[[[jet,jet]],[[jet,jet]]]" )
        b0=element.branches[0]
        sb0=str(b0)
        assert ( sb0 == "[[jet,jet]]" )
        assert ( b0.masses[0]-475*GeV ) < .1*GeV
        assert ( b0.masses[1]-325*GeV ) < .1*GeV

if __name__ == "__main__":
    unittest.main()
