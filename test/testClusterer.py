#!/usr/bin/env python3

"""
.. module:: testClusterer
   :synopsis: Tests the mass clusterer
              Depends also on lheReader, lheDecomposer.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys,os
sys.path.insert(0,"../")
import unittest

class ClustererTest(unittest.TestCase):
    def testClusterer(self):
        """ test the mass clusterer """
        from smodels.theory.element import Element
        from smodels.theory.branch import Branch
        from smodels.theory import clusterTools
        from smodels.experiment.txnameObj import TxName, TxNameData
        from smodels.experiment.infoObj import Info
        from smodels.share.models import SMparticles, mssm
        from smodels.theory.crossSection import XSection,XSectionInfo,XSectionList        
        from smodels.tools.physicsUnits import GeV, TeV, fb
        import copy
        
        data = [[ [[ 674.99*GeV, 199.999*GeV], [ 674.99*GeV, 199.999*GeV] ],.03*fb ], 
               [ [[ 725.0001*GeV,200.*GeV], [ 725.0001*GeV,200.*GeV] ], .06*fb ] ,
               [ [[ 750.*GeV,250.*GeV], [ 750.*GeV,250.*GeV] ], .03*fb ] ]
        info = Info(os.path.join("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/data/","dataInfo.txt"))
        globalInfo = Info(os.path.join("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/","globalInfo.txt"))        
        txnameData=TxNameData(data, "efficiencyMap", Id=1)
        txname=TxName("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/data/T2bb.txt",globalInfo,info)
        txname.txnameData = txnameData
        
        u = SMparticles.u
        gluino = mssm.gluino
        gluino.__setattr__("mass", 675.*GeV)
        n1 = mssm.n1
        n1.__setattr__("mass", 200.*GeV)
        

        w1 = XSectionList()
        w1.xSections.append(XSection())
        w1.xSections[0].info = XSectionInfo()
        w1.xSections[0].info.sqrts = 8.*TeV
        w1.xSections[0].info.label = '8 TeV'
        w1.xSections[0].info.order = 0
        w1.xSections[0].value = 10.*fb

        b1 = Branch()
        b1.evenParticles = [[u,u]]
        b1.BSMparticles = [gluino, n1]
        b2 = b1.copy()
        el1 = Element()
        el1.branches=[b1,b2]  
        el1.weight = w1       
        el1.txname = txname           


        ## make a second element with a slightly different gluino mass
        el2=copy.deepcopy(el1)
        el2.branches[0].BSMparticles[0].__setattr__("mass", 725.*GeV) 
        el2.branches[1].BSMparticles[0].__setattr__("mass", 725.*GeV)  
        el2.txname = txname

        # lets now cluster the two different gluino masses.
        newel=clusterTools.groupAll ( [el1,el2] )
        newmasses=newel.getAvgMass()
        self.assertTrue ( newmasses==None ) ## in the case of efficiency maps the avg mass is none
        ## since it makes no sense

        txname.txnameData.dataTag = 'upperLimits'
        newel=clusterTools.clusterElements ( [el1,el2], 5. )
        ## this example gives an avg cluster mass of 700 gev
        self.assertTrue ( newel[0].getAvgMass()[0][0] == 700. * GeV )
        
        newel=clusterTools.clusterElements ( [el1,el2], .5 )
        #in this example the distance is not in maxdist, so we dont cluster
        self.assertTrue ( len(newel)==2 )

if __name__ == "__main__":
    unittest.main()
