#!/usr/bin/env python

"""
.. module:: testClusterer
   :synopsis: Tests the mass clusterer
              Depends also on lheReader, lheDecomposer.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest

class ClustererTest(unittest.TestCase):
    def testClusterer(self):
        """ test the mass clusterer """
        from smodels.theory import lheReader, lheDecomposer, crossSection
        from smodels.theory import clusterTools
        from smodels.experiment.txnameObj import TxName, TxNameData
        from smodels.installation import installDirectory
        from smodels.tools.physicsUnits import GeV, pb, fb
        import copy
        
        data = [ [ [[ 674.99*GeV, 199.999*GeV], [ 674.99*GeV, 199.999*GeV] ],  .03*fb ], 
               [ [[ 725.0001*GeV,200.*GeV], [ 725.0001*GeV,200.*GeV] ], .06*fb ] ,
               [ [[ 750.*GeV,250.*GeV], [ 750.*GeV,250.*GeV] ], .03*fb ] ]
        txnameData=TxNameData(data)
        txname=TxName("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/data/T2bb.txt","info")
        txname.txnameData = txnameData
        txname.txnameData.dataTag = 'efficiencyMap'

        filename = "%sinputFiles/lhe/simplyGluino.lhe" % (installDirectory() )
        reader = lheReader.LheReader(filename)
        event = reader.next()
        event_xsec=event.metainfo["totalxsec"]
        self.assertTrue ( abs ( event_xsec - 0.262 * pb ) < .1 *pb )
        xsecs = crossSection.getXsecFromLHEFile(filename)
        element = lheDecomposer.elementFromEvent(event, xsecs )
        element.txname=None
        e0=copy.deepcopy(element) ## has a gluino with mass of 675 GeV


        ## make a second element with a slightly different gluino mass
        e1=copy.deepcopy(element)
        e1.branches[0].masses[0]=725*GeV
        e1.branches[1].masses[0]=725*GeV
        e0.txname = txname
        e1.txname = txname

        # lets now cluster the two different gluino masses.
        newel=clusterTools.groupAll ( [e0,e1] )
        newmasses=newel.getAvgMass()
        self.assertTrue ( newmasses==None ) ## in the case of efficiency maps the avg mass is none
        ## since it makes no sense

        txname.txnameData.dataTag = 'upperLimits'
        newel=clusterTools.clusterElements ( [e0,e1], 5. )
        ## this example gives an avg cluster mass of 700 gev
        self.assertTrue ( newel[0].getAvgMass()[0][0] == 700. * GeV )
        
        newel=clusterTools.clusterElements ( [e0,e1], .5 )
        #in this example the distance is not in maxdist, so we dont cluster
        self.assertTrue ( len(newel)==2 )

if __name__ == "__main__":
    unittest.main()
