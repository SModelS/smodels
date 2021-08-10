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
from smodels.theory import decomposer
from smodels.particlesLoader import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model
from smodels.theory import clusterTools
from smodels.experiment.txnameObj import TxName, TxNameData
from smodels.experiment.datasetObj import DataSet
from smodels.experiment.infoObj import Info
from smodels.share.models import SMparticles, mssm
from smodels.theory.crossSection import XSection,XSectionInfo,XSectionList
from smodels.tools.physicsUnits import GeV, TeV, fb
from databaseLoader import database
from smodels.theory.clusterTools import clusterElements
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.theory.particle import ParticleList
from smodels.experiment.defaultFinalStates import finalStates
from smodels.theory.auxiliaryFunctions import roundValue


class ClustererTest(unittest.TestCase):

    def testSimpleCluster(self):
        """ test the mass clusterer """

        data = [[ [[ 674.99*GeV, 199.999*GeV], [ 674.99*GeV, 199.999*GeV] ],.03*fb ],
               [ [[ 725.0001*GeV,200.*GeV], [ 725.0001*GeV,200.*GeV] ], .06*fb ] ,
               [ [[ 750.*GeV,250.*GeV], [ 750.*GeV,250.*GeV] ], .03*fb ] ]
        info = Info("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/data/dataInfo.txt")
        globalInfo = Info("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/globalInfo.txt")
        txnameData=TxNameData(data, "upperLimit", Id=1)
        txname=TxName("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/data/T2bb.txt",
                        globalInfo,info,finalStates)
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
        el1.eff = 1. #(Used in clustering)


        ## make a second element with a slightly different gluino mass
        el2= el1.copy()
        el2.motherElements = [el2] #Enforce el2 and el1 not to be related
        el2.txname = txname
        el2.branches[0].oddParticles = [ptc.copy() for ptc in el1.branches[0].oddParticles]
        el2.branches[1].oddParticles = [ptc.copy() for ptc in el1.branches[1].oddParticles]
        el2.branches[0].oddParticles[0].__setattr__("mass", 725.*GeV)
        el2.branches[1].oddParticles[0].__setattr__("mass", 725.*GeV)
        el2.eff = 1. #(Used in clustering)

        #Cluster for upper limits (all elements close in upper limit should be clustered together)
        maxDist = 5. #Cluster all elements
        newel=clusterTools.clusterElements([el1,el2], maxDist, dataset)[0]
        newmasses=newel.averageElement().mass
        self.assertEqual(newmasses,[[700.*GeV,200.*GeV]]*2)

        maxDist = 0.5 #Elements differ and should not be clustered
        newel=clusterTools.clusterElements([el1,el2], maxDist, dataset)
        #in this example the distance is not in maxdist, so we dont cluster
        self.assertTrue(len(newel)==2)



        info = Info("./database/8TeV/CMS/CMS-SUS-13-012-eff/6NJet8_1000HT1250_200MHT300/dataInfo.txt")
        globalInfo = Info("./database/8TeV/CMS/CMS-SUS-13-012-eff/globalInfo.txt")
        txnameData=TxNameData(data, "efficiencyMap", Id=1)
        txname=TxName("./database/8TeV/CMS/CMS-SUS-13-012-eff/6NJet8_1000HT1250_200MHT300/T2.txt",
                        globalInfo,info,finalStates)
        txname.txnameData = txnameData
        dataset = DataSet(info = globalInfo, createInfo = False)
        dataset.dataInfo = info
        dataset.txnameList = [txname]

        #Cluster for efficiency maps (all elements should be clustered together independent of maxDist)
        maxDist = 0.001
        newel=clusterTools.clusterElements([el1,el2],maxDist,dataset)[0]
        newmasses=newel.averageElement().mass
        self.assertEqual(newmasses,[[700.*GeV,200.*GeV]]*2)

    def testClusteringEM(self):

        slhafile = 'testFiles/slha/lightEWinos.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(slhafile)
        sigmacut = 5.*fb
        mingap = 5.*GeV
        toplist = decomposer.decompose(model, sigmacut, doCompress=True, doInvisible=True, minmassgap=mingap)

        #Test clustering for EM results
        dataset = database.getExpResults(analysisIDs='CMS-SUS-13-012',
                                         datasetIDs='3NJet6_800HT1000_300MHT450')[0].getDataset('3NJet6_800HT1000_300MHT450')


        el1 = toplist[0].elementList[0].copy()
        el2 = toplist[0].elementList[1].copy()
        el3 = toplist[2].elementList[1].copy()
        el1.eff = 1. #(Used in clustering)
        el2.eff = 1. #(Used in clustering)
        el3.eff = 1. #(Used in clustering)
        #All elements have the same UL (for EM results)
        el1._upperLimit = el2._upperLimit = el3._upperLimit = 1.*fb
        #Clustering should not depend on the mass, width or txname:
        el1.txname = el2.txname = el3.txname = None
        el1.mass = [[1000.*GeV,10*GeV]]*2
        el2.mass = [[1500.*GeV,10*GeV]]*2
        el3.mass = [[200.*GeV,100*GeV,90.*GeV]]*2
        el3.totalwidth = [[1e-10*GeV,1e-15*GeV,0.*GeV]]*2
        clusters = clusterElements([el1,el2,el3],maxDist=0.2,dataset=dataset)
        self.assertEqual(len(clusters),1)
        self.assertEqual(sorted(clusters[0].elements),sorted([el1,el2,el3]))
        self.assertEqual(clusters[0].averageElement().mass,None)
        self.assertEqual(clusters[0].averageElement().totalwidth,None)

    def testClusteringUL(self):

        slhafile = 'testFiles/slha/lightEWinos.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(slhafile)
        sigmacut = 5.*fb
        mingap = 5.*GeV
        toplist = decomposer.decompose(model, sigmacut, doCompress=True, doInvisible=True, minmassgap=mingap)

        #Test clustering for UL results
        dataset = database.getExpResults(analysisIDs='ATLAS-SUSY-2013-02',datasetIDs=None)[0].getDataset(None)

        el1 = toplist[1].elementList[0]
        el2 = toplist[1].elementList[1]
        el3 = toplist[1].elementList[2]
        weights = [el.weight.getMaxXsec().asNumber(fb) for el in [el1,el2,el3]]

        #Overwrite masses and txname label
        el1.mass = [[1000.*GeV,100.*GeV]]*2
        el2.mass = [[1020.*GeV,100.*GeV]]*2
        el3.mass = [[500.*GeV,100.*GeV]]*2
        el1.txname = el2.txname = el3.txname = 'T1'
        el1.eff = 1. #(Used in clustering)
        el2.eff = 1. #(Used in clustering)
        el3.eff = 1. #(Used in clustering)

        #Check clustering with distinct elements
        clusters = clusterElements([el1,el2,el3],maxDist=0.2,dataset=dataset)
        self.assertEqual(len(clusters),2)
        averageMasses = [[[(1000.*GeV*weights[0]+1020.*GeV*weights[1])/(weights[0]+weights[1]),100.*GeV]]*2, el3.mass]
        for ibr,br in enumerate(averageMasses):
            for im,m in enumerate(br):
                averageMasses[ibr][im] = roundValue(m,nround=5)

        elClusters = [[el1,el2], [el3]]
        for ic,cluster in enumerate(clusters):
            avgEl = cluster.averageElement()
            self.assertEqual(sorted(cluster.elements),sorted(elClusters[ic]))
            for ibr,br in enumerate(avgEl.mass):
                for im,m in enumerate(br):
                    self.assertAlmostEqual(m.asNumber(GeV),averageMasses[ic][ibr][im].asNumber(GeV))
            self.assertEqual(avgEl.totalwidth,elClusters[ic][0].totalwidth)

        #Check clustering with distinct elements but no maxDist limit
        clusters = clusterElements([el1,el2,el3],maxDist=10.,dataset=dataset)
        self.assertEqual(len(clusters),1)
        cluster = clusters[0]
        avgEl = cluster.averageElement()
        averageMass = [[(1000.*GeV*weights[0]+1020.*GeV*weights[1] + 500.*GeV*weights[2])/sum(weights),100.*GeV]]*2
        for ibr,br in enumerate(averageMass):
            for im,m in enumerate(br):
                averageMass[ibr][im] = roundValue(m,nround=5)

        self.assertEqual(sorted(cluster.elements),sorted([el1,el2,el3]))
        for ibr,br in enumerate(avgEl.mass):
            for im,m in enumerate(br):
                self.assertAlmostEqual(m.asNumber(GeV),averageMass[ibr][im].asNumber(GeV))
        self.assertEqual(avgEl.totalwidth,el1.totalwidth)


        #Check clustering where elements have same upper limits, but not the average element:
        el1._upperLimit = 1.*fb
        el2._upperLimit = 1.*fb
        el3._upperLimit = 1.*fb
        clusters = clusterElements([el1,el2,el3],maxDist=0.1,dataset=dataset)
        self.assertEqual(len(clusters),2)

    def testClustererLifeTimes(self):
        """ test the clustering with distinct lifetimes"""


        data = [[ [[ 674.99*GeV, 199.999*GeV], [ 674.99*GeV, 199.999*GeV] ],.03*fb ],
               [ [[ 725.0001*GeV,200.*GeV], [ 725.0001*GeV,200.*GeV] ], .06*fb ] ,
               [ [[ 750.*GeV,250.*GeV], [ 750.*GeV,250.*GeV] ], .03*fb ] ]
        info = Info("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/data/dataInfo.txt")
        globalInfo = Info("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/globalInfo.txt")
        txnameData=TxNameData(data, "upperLimit", Id=1)
        txname=TxName("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/data/T2bb.txt",
                        globalInfo,info,finalStates)
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
        b1.evenParticles = [ParticleList([u,u])]
        b1.oddParticles = [gluino, n1]
        b2 = b1.copy()
        el1 = Element()
        el1.branches=[b1,b2]
        el1.weight = w1
        el1.txname = txname
        el1.eff = 1. #(Used in clustering)


        ## make a second element with a slightly different gluino width
        el2 = el1.copy()
        el2.motherElements = [el2]  #Enforce el2 and el1 not to be related
        el2.txname = txname
        el2.branches[0].oddParticles = [ptc.copy() for ptc in el1.branches[0].oddParticles]
        el2.branches[1].oddParticles = [ptc.copy() for ptc in el1.branches[1].oddParticles]
        el2.eff = 1. #(Used in clustering)

        el2.branches[0].oddParticles[0].__setattr__("mass", 675.*GeV)
        el2.branches[1].oddParticles[0].__setattr__("mass", 675.*GeV)
        el2.branches[0].oddParticles[0].__setattr__("totalwidth", 0.9e-15*GeV)
        el2.branches[1].oddParticles[0].__setattr__("totalwidth", 0.9e-15*GeV)

        newel=clusterTools.clusterElements([el1,el2], 5., dataset)
        ## this example gives an avg cluster mass of 700 gev
        self.assertEqual(newel[0].averageElement().mass[0][0],675.*GeV)
        self.assertAlmostEqual(newel[0].averageElement().totalwidth[0][0].asNumber(GeV)*1e15,0.95)

        newel=clusterTools.clusterElements([el1,el2], .5, dataset)
        #in this example the distance is in maxdist, so we cluster
        self.assertTrue(len(newel)==1)

        newel=clusterTools.clusterElements([el1,el2], .1, dataset)
        #in this example the distance is not in maxdist, so we dont cluster
        self.assertTrue(len(newel)==2)

    def testComplexCluster(self):
        """ test the mass clusterer """

        slhafile = 'testFiles/slha/416126634.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(slhafile)
        sigmacut = 0.03*fb
        mingap = 5.*GeV
        toplist = decomposer.decompose(model, sigmacut, doCompress=True, doInvisible=True, minmassgap=mingap)

        #Test clustering for UL results
        expResult = database.getExpResults(analysisIDs='CMS-SUS-16-039',dataTypes='upperLimit')[0]
        predictions = theoryPredictionsFor(expResult, toplist, combinedResults=False, marginalize=False)
        clusterSizes = sorted([len(p.elements) for p in predictions])
        self.assertEqual(clusterSizes, [1,16,24])

if __name__ == "__main__":
    unittest.main()
