#!/usr/bin/env python

"""
.. module:: testClusterer
   :synopsis: Tests the mass clusterer
              Depends also on lheReader, lheDecomposer.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import unittest

class ClustererTest(unittest.TestCase):
    def testClusterer(self):
        """ test the mass clusterer """
        from smodels.theory import lheReader, lheDecomposer, crossSection
        from smodels.theory import clusterTools
        from smodels.installation import installDirectory
        from smodels.tools.physicsUnits import GeV, pb
        import copy

        filename = "%sinputFiles/lhe/T1_1.lhe" % (installDirectory() )
        reader = lheReader.LheReader(filename)
        event = reader.next()
        event_xsec=event.metainfo["totalxsec"]
        assert  ( event_xsec - 84.*pb ) < .1 *pb
        xsecs = crossSection.getXsecFromLHEFile(filename)
        element = lheDecomposer.elementFromEvent(event, xsecs )
                  #crossSection.XSectionList( { "8 TeV (NLL)": xsec } ))
        ## print "weight",element.weight
        e0=copy.deepcopy(element)
        ## make a second element with a slightly different gluino mass
        e1=copy.deepcopy(element)
        e1.branches[0].masses[0]=473*GeV
        e1.branches[1].masses[0]=473*GeV

        # lets now cluster the two different gluino masses.
        # yes, this is a very strange example :)
        newel=clusterTools.groupAll ( [e0,e1] )
        newmasses=newel.getAvgMass()
        self.assertAlmostEquals ( newmasses[0][0]/GeV, 474.008385295 ) 
        self.assertAlmostEquals ( newmasses[0][1]/GeV, 325. )

if __name__ == "__main__":
    unittest.main()
