#!/usr/bin/env python3

"""
.. module:: testClusterer
   :synopsis: Tests the mass clusterer
              Depends also on lheReader, lheDecomposer.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest

from smodels.theory.element import Element
from smodels.theory.branch import Branch
from smodels.theory import clusterTools
from smodels.experiment.txnameObj import TxName, TxNameData
from smodels.experiment.datasetObj import DataSet
from smodels.experiment.infoObj import Info
from smodels.share.models import SMparticles, mssm
from smodels.theory.crossSection import XSection,XSectionInfo,XSectionList        
from smodels.tools.physicsUnits import GeV, TeV, fb
import copy


class ClustererTest(unittest.TestCase):
    
    def testClusterer(self):
        """ test the mass clusterer """
        
        data = [[ [[ 674.99*GeV, 199.999*GeV], [ 674.99*GeV, 199.999*GeV] ],.03*fb ], 
               [ [[ 725.0001*GeV,200.*GeV], [ 725.0001*GeV,200.*GeV] ], .06*fb ] ,
               [ [[ 750.*GeV,250.*GeV], [ 750.*GeV,250.*GeV] ], .03*fb ] ]
        info = Info("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/data/dataInfo.txt")
        globalInfo = Info("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/globalInfo.txt")
        txnameData=TxNameData(data, "upperLimit", Id=1)
        txname=TxName("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/data/T2bb.txt",globalInfo,info)
        txname.txnameData = txnameData
        dataset = DataSet(info = globalInfo, createInfo = False)
        dataset.dataInfo = info
        dataset.txnameList = [txname]
        
        u = SMparticles.u
        gluino = mssm.gluino.copy()
        gluino.__setattr__("mass", 675.*GeV)
        gluino.__setattr__('totalwidth',float('inf')*GeV)
        n1 = mssm.n1.copy()
        n1.__setattr__("mass", 200.*GeV)
        n1.__setattr__('totalwidth',0.*GeV)
        

        w1 = XSectionList()
        w1.xSections.append(XSection())
        w1.xSections[0].info = XSectionInfo()
        w1.xSections[0].info.sqrts = 8.*TeV
        w1.xSections[0].info.label = '8 TeV'
        w1.xSections[0].info.order = 0
        w1.xSections[0].value = 10.*fb

        b1 = Branch()
        b1.evenParticles = [[u,u]]
        b1.oddParticles = [gluino, n1]
        b2 = b1.copy()
        el1 = Element()
        el1.branches=[b1,b2]  
        el1.weight = w1       
        el1.txname = txname           


        ## make a second element with a slightly different gluino mass
        el2=copy.deepcopy(el1)
        el2.txname = txname
        el2.branches[0].oddParticles[0].__setattr__("mass", 725.*GeV) 
        el2.branches[1].oddParticles[0].__setattr__("mass", 725.*GeV)

        #Cluster for upper limits (all elements close in upper limit should be clustered together)
        maxDist = 5. #Cluster all elements
        newel=clusterTools.clusterElements([el1,el2], maxDist, dataset)[0]
        newmasses=newel.getAvgMass()
        self.assertEqual(newmasses,[[700.*GeV,200.*GeV]]*2)
        
        maxDist = 0.5 #Elements differ and should not be clustered
        newel=clusterTools.clusterElements([el1,el2], maxDist, dataset)
        #in this example the distance is not in maxdist, so we dont cluster
        self.assertTrue(len(newel)==2)



        info = Info("./database/8TeV/CMS/CMS-SUS-13-012-eff/6NJet8_1000HT1250_200MHT300/dataInfo.txt")
        globalInfo = Info("./database/8TeV/CMS/CMS-SUS-13-012-eff/globalInfo.txt")
        txnameData=TxNameData(data, "efficiencyMap", Id=1)
        txname=TxName("./database/8TeV/CMS/CMS-SUS-13-012-eff/6NJet8_1000HT1250_200MHT300/T2.txt",globalInfo,info)
        txname.txnameData = txnameData
        dataset = DataSet(info = globalInfo, createInfo = False)
        dataset.dataInfo = info
        dataset.txnameList = [txname]

        #Cluster for efficiency maps (all elements should be clustered together independent of maxDist)
        maxDist = 0.001
        newel=clusterTools.clusterElements([el1,el2],maxDist,dataset)[0]
        newmasses=newel.getAvgMass()
        self.assertEqual(newmasses,[[700.*GeV,200.*GeV]]*2)


    def testClustererLifeTimes(self):
        """ test the clustering with distinct lifetimes"""
        
        
        data = [[ [[ 674.99*GeV, 199.999*GeV], [ 674.99*GeV, 199.999*GeV] ],.03*fb ], 
               [ [[ 725.0001*GeV,200.*GeV], [ 725.0001*GeV,200.*GeV] ], .06*fb ] ,
               [ [[ 750.*GeV,250.*GeV], [ 750.*GeV,250.*GeV] ], .03*fb ] ]
        info = Info("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/data/dataInfo.txt")
        globalInfo = Info("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/globalInfo.txt")
        txnameData=TxNameData(data, "upperLimit", Id=1)
        txname=TxName("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/data/T2bb.txt",globalInfo,info)
        txname.txnameData = txnameData
        dataset = DataSet(info = globalInfo, createInfo = False)
        dataset.dataInfo = info
        dataset.txnameList = [txname]
        
        u = SMparticles.u
        gluino = mssm.gluino.copy()
        gluino.__setattr__("mass", 675.*GeV)
        gluino.__setattr__('totalwidth',1e-15*GeV)
        n1 = mssm.n1.copy()
        n1.__setattr__("mass", 200.*GeV)
        n1.__setattr__('totalwidth',0.*GeV)
        

        w1 = XSectionList()
        w1.xSections.append(XSection())
        w1.xSections[0].info = XSectionInfo()
        w1.xSections[0].info.sqrts = 8.*TeV
        w1.xSections[0].info.label = '8 TeV'
        w1.xSections[0].info.order = 0
        w1.xSections[0].value = 10.*fb

        b1 = Branch()
        b1.evenParticles = [[u,u]]
        b1.oddParticles = [gluino, n1]
        b2 = b1.copy()
        el1 = Element()
        el1.branches=[b1,b2]  
        el1.weight = w1       
        el1.txname = txname           


        ## make a second element with a slightly different gluino width
        el2=copy.deepcopy(el1)
        el2.txname = txname

        el2.branches[0].oddParticles[0].__setattr__("mass", 675.*GeV) 
        el2.branches[1].oddParticles[0].__setattr__("mass", 675.*GeV)
        el2.branches[0].oddParticles[0].__setattr__("totalwidth", 0.9e-15*GeV) 
        el2.branches[1].oddParticles[0].__setattr__("totalwidth", 0.9e-15*GeV)
        
        newel=clusterTools.clusterElements([el1,el2], 5., dataset)
        ## this example gives an avg cluster mass of 700 gev
        self.assertEqual(newel[0].getAvgMass()[0][0],675.*GeV)
        self.assertAlmostEqual(newel[0].getAvgWidth()[0][0].asNumber(GeV)*1e15,0.95)
        
        newel=clusterTools.clusterElements([el1,el2], .5, dataset)
        #in this example the distance is in maxdist, so we cluster
        self.assertTrue(len(newel)==1)

        newel=clusterTools.clusterElements([el1,el2], .1, dataset)
        #in this example the distance is not in maxdist, so we dont cluster
        self.assertTrue(len(newel)==2)


if __name__ == "__main__":
    unittest.main()
