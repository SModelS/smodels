#!/usr/bin/env python

"""
.. module:: testClusterer
   :synopsis: Tests the mass clusterer
              Depends also on lheReader, lheDecomposer.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import unittest
import sys
sys.path.append('../')

class ClustererTest(unittest.TestCase):
    def testClusterer(self):
        """ test the mass clusterer """
        from smodels.theory import lheReader, lheDecomposer, crossSection
        from smodels.theory import clusterTools
        from smodels.installation import installDirectory
        from smodels.tools.physicsUnits import GeV, pb
        import copy

        filename = "%sinputFiles/lhe/simplyGluino.lhe" % (installDirectory() )
        reader = lheReader.LheReader(filename)
        event = reader.next()
        event_xsec=event.metainfo["totalxsec"]
        self.assertTrue ( abs ( event_xsec - 0.262 * pb ) < .1 *pb )
        xsecs = crossSection.getXsecFromLHEFile(filename)
        element = lheDecomposer.elementFromEvent(event, xsecs )
                  #crossSection.XSectionList( { "8 TeV (NLL)": xsec } ))
        element.txname=None
        # print "w0=",element.branches[0].masses
        e0=copy.deepcopy(element)

        ## make a second element with a slightly different gluino mass
        e1=copy.deepcopy(element)
        e1.branches[0].masses[0]=725*GeV
        e1.branches[1].masses[0]=725*GeV

        # lets now cluster the two different gluino masses.
        # yes, this is a very strange example :)
        newel=clusterTools.groupAll ( [e0,e1] )
        newel.txname=True
        newmasses=newel.getAvgMass()
        print "newmasses=",newmasses
        self.assertAlmostEquals ( newmasses[0][0]/GeV, 700. ) 
        self.assertAlmostEquals ( newmasses[0][1]/GeV, 200. )

if __name__ == "__main__":
    unittest.main()
