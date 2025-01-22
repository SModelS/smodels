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

from smodels.decomposition import decomposer
from smodels.matching import clusterTools
from smodels.experiment.txnameObj import TxName
from smodels.experiment.txnameDataObj import TxNameData
from smodels.experiment.datasetObj import DataSet
from smodels.experiment.infoObj import Info
from smodels.base.crossSection import XSection,XSectionInfo,XSectionList
from smodels.base.physicsUnits import GeV, TeV, fb
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model
import numpy as np
from databaseLoader import database
from smodels.matching.theoryPrediction import theoryPredictionsFor
from smodels.experiment.defaultFinalStates import finalStates
from smodels.matching.matcherAuxiliaryFuncs import roundValue
from smodels.base.smodelsLogging import setLogLevel
setLogLevel("error")


from unitTestHelpers import theorySMSFromString as fromString
slhafile = './testFiles/slha/lightEWinos.slha'

model = Model( BSMparticles=BSMList, SMparticles=SMList)
model.updateParticles(inputFile=slhafile,ignorePromptQNumbers=['spin'])


class ClustererTest(unittest.TestCase):

    def testSimpleCluster(self):
        """ test the mass clusterer """

        data = [[ [[ 674.99, 199.999], [ 674.99, 199.999] ],.03 ],
               [ [[ 725.0001,200.], [ 725.0001,200.] ], .06 ] ,
               [ [[ 750.,250.], [ 750.,250.] ], .03 ] ]
        xvalues = [np.array(pt[0]).flatten() for pt in data]
        yvalues = [pt[1] for pt in data]
        info = Info("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/data/dataInfo.txt")
        globalInfo = Info("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/globalInfo.txt")
        txnameData=TxNameData(x=xvalues,y=yvalues, txdataId=1)
        txname=TxName("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/data/T2bb.txt",
                        globalInfo,info,finalStates)
        txname.txnameData = txnameData
        dataset = DataSet(info = globalInfo, createInfo = False)
        dataset.dataInfo = info
        dataset.txnameList = [txname]


        w1 = XSectionList()
        w1.xSections.append(XSection())
        w1.xSections[0].info = XSectionInfo()
        w1.xSections[0].info.sqrts = 8.*TeV
        w1.xSections[0].info.label = '8 TeV'
        w1.xSections[0].info.order = 0
        w1.xSections[0].value = 10.*fb

        sms = fromString("(PV > sb_1(1),sb_1(2)), (sb_1(1) > b,N1), (sb_1(2) > b,N1)",
                         model=model)

        sms =txname.hasSMSas(sms)
        sms.eff = 1.0
        sms.txname = txname
        sms.weight = w1.getMaxXsec()

        sb1 = model.getParticle(label='sb_1')
        sb1.mass = 675.*GeV
        sb1.totalwidth = float('inf')*GeV
        n1 = model.getParticle(label='N1')
        n1.mass = 200.*GeV
        n1.totalwidth = 0.*GeV


        # make a second SMS with a slightly different sbottom mass
        sms2 = fromString("(PV > sb_2(1),sb_2(2)), (sb_2(1) > b,N1), (sb_2(2) > b,N1)",
                         model=model)

        sms2 =txname.hasSMSas(sms2)
        sms2.eff = 1.0
        sms2.txname = txname
        sms2.weight = w1.getMaxXsec()

        sb2 = model.getParticle(label='sb_2')
        sb2.mass = 725*GeV
        sb2.totalwidth = float('inf')*GeV


        #Cluster for upper limits (all elements close in upper limit should be clustered together)
        maxDist = 5. #Cluster all elements
        cluster = clusterTools.clusterSMS([sms,sms2], maxDist, dataset)
        self.assertEqual(len(cluster),1)
        cluster = cluster[0]
        newmass = [node.mass for node in cluster.averageSMS.nodes if not node.isSM]
        self.assertEqual(newmass,[700.*GeV,700.*GeV,200.*GeV,200.*GeV])

        maxDist = 0.5 #Elements differ and should not be clustered
        cluster = clusterTools.clusterSMS([sms,sms2], maxDist, dataset)
        #in this example the distance is not in maxdist, so we dont cluster
        self.assertEqual(len(cluster),2)



        info = Info("./database/8TeV/CMS/CMS-SUS-13-012-eff/6NJet8_1000HT1250_200MHT300/dataInfo.txt")
        globalInfo = Info("./database/8TeV/CMS/CMS-SUS-13-012-eff/globalInfo.txt")
        txnameData=TxNameData(x=xvalues,y=yvalues, txdataId=1)
        txname = TxName("./database/8TeV/CMS/CMS-SUS-13-012-eff/6NJet8_1000HT1250_200MHT300/T2.txt",
                        globalInfo,info,finalStates)
        txname.txnameData = txnameData
        dataset = DataSet(info = globalInfo, createInfo = False)
        dataset.dataInfo = info
        dataset.txnameList = [txname]

        sms = fromString("(PV > sb_1(1),sb_1(2)), (sb_1(1) > q,N1), (sb_1(2) > q,N1)",
                         model=model)
        sms =txname.hasSMSas(sms)
        sms.eff = 1.0
        sms.txname = txname
        sms.weight = w1.getMaxXsec()
        # make a second SMS with a slightly different sbottom mass
        sms2 = fromString("(PV > sb_2(1),sb_2(2)), (sb_2(1) > q,N1), (sb_2(2) > q,N1)",
                         model=model)

        sms2 =txname.hasSMSas(sms2)
        sms2.eff = 1.0
        sms2.txname = txname
        sms2.weight = w1.getMaxXsec()

        sb2 = model.getParticle(label='sb_2')
        sb2.mass = 725*GeV
        sb2.totalwidth = float('inf')*GeV


        #Cluster for efficiency maps (all elements should be clustered together independent of maxDist)
        maxDist = 0.001
        cluster = clusterTools.clusterSMS([sms,sms2],maxDist,dataset)
        self.assertEqual(len(cluster),1)
        cluster = cluster[0]
        newmass = [node.mass for node in cluster.averageSMS.nodes if not node.isSM]
        self.assertEqual(newmass,[700.*GeV,700.*GeV,200.*GeV,200.*GeV])

    def testClusteringEM(self):

        slhafile = 'testFiles/slha/lightEWinos.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(slhafile)
        sigmacut = 5.*fb
        mingap = 5.*GeV
        topDict= decomposer.decompose(model, sigmacut, massCompress=True,
                                       invisibleCompress=True, minmassgap=mingap)

        #Test clustering for EM results
        database.selectExpResults(analysisIDs='CMS-SUS-13-012',
                                         datasetIDs='3NJet6_800HT1000_300MHT450')
        dataset = database.expResultList[0]
        dataset = dataset.getDataset('3NJet6_800HT1000_300MHT450')


        sms1 = topDict[110110101000][0]
        sms2 = topDict[110110101000][1]
        sms3 = topDict[11101010011011010000][1]

        sms1.eff = 1. #(Used in clustering)
        sms2.eff = 1. #(Used in clustering)
        sms3.eff = 1. #(Used in clustering)
        sms1.weight = 10.*fb #(Used in clustering)
        sms2.weight = 10.*fb #(Used in clustering)
        sms3.weight = 10.*fb #(Used in clustering)


        #All smsements have the same UL (for EM results)
        sms1._upperLimit = sms2._upperLimit = sms3._upperLimit = 1.*fb
        #Clustering should not depend on the mass, width or txname:
        sms1.txname = sms2.txname = sms3.txname = None
        clusters = clusterTools.clusterSMS([sms1,sms2,sms3],maxDist=0.2,dataset=dataset)
        self.assertEqual(len(clusters),1)
        self.assertEqual(sorted(clusters[0].smsList),sorted([sms1,sms2,sms3]))
        self.assertEqual(clusters[0].averageSMS,None)  # SMS have distinct topologies, no averageSMS can be built

    def testClusteringUL(self):

        slhafile = 'testFiles/slha/lightEWinos.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(slhafile)

        sms1 = fromString("(PV > N2(1),N2(2)), (N2(1) > u,u,N1), (N2(2) > u,u,N1)",
                         model=model)
        sms2 = fromString("(PV > N3(1),N3(2)), (N3(1) > d,d,N1), (N3(2) > d,d,N1)",
                         model=model)
        sms3 = fromString("(PV > N4(1),N4(2)), (N4(1) > c,c,N1), (N4(2) > c,c,N1)",
                         model=model)

        smsList = [sms1,sms2,sms3]
        n1 = model.getParticle(label='N1')
        n2 = model.getParticle(label='N2')
        n3 = model.getParticle(label='N3')
        n4 = model.getParticle(label='N4')
        n1.mass = 100*GeV
        n2.mass = 1000*GeV
        n3.mass = 1020*GeV
        n4.mass = 500*GeV

        #Test clustering for UL results
        database.selectExpResults(analysisIDs='ATLAS-SUSY-2013-02',datasetIDs=None)
        dataset = database.expResultList[0]
        dataset = dataset.getDataset(None)
        tx = [t for t in dataset.txnameList if str(t) == 'T1'][0]
        smsList = [tx.hasSMSas(sms) for sms in smsList[:]]

        for isms,sms in enumerate(smsList):
            sms.eff = 1.0
            sms.weight = 1.5*(isms+1)*1.0*fb
            sms.txname = tx
        weights = [sms.weight.asNumber(fb) for sms in smsList]

        #Check clustering with distinct SMS
        clusters = clusterTools.clusterSMS(smsList,maxDist=0.2,dataset=dataset)
        self.assertEqual(len(clusters),2)
        avgMass1 = (1000.*GeV*weights[0]+1020.*GeV*weights[1])/(weights[0]+weights[1])
        avgMass1 = [avgMass1,avgMass1,100*GeV,100*GeV]
        avgMass1 = [roundValue(m.asNumber(GeV),5)for m in avgMass1]
        avgMass2 = [node.mass for node in sms3.nodes if not node.isSM]
        avgMass2 = [roundValue(m.asNumber(GeV),5) for m in avgMass2]
        averageMasses = [avgMass1, avgMass2]

        totalwidth1 = (n2.totalwidth*weights[0]+n3.totalwidth*weights[1])/(weights[0]+weights[1])
        totalwidth1 = [totalwidth1,totalwidth1,0*GeV,0*GeV]
        totalwidth1 = [roundValue(m.asNumber(GeV),5)for m in totalwidth1]
        totalwidth2 = [node.totalwidth for node in sms3.nodes if not node.isSM]
        totalwidth2 = [roundValue(m.asNumber(GeV),5) for m in totalwidth2]
        averageWidths = [totalwidth1, totalwidth2]


        smsClusters = [smsList[:2], [smsList[2]]]
        for ic,cluster in enumerate(clusters):
            avgSMS = cluster.averageSMS
            self.assertEqual(sorted(cluster.smsList),sorted(smsClusters[ic]))
            avgMass = [node.mass for node in avgSMS.nodes if not node.isSM]
            avgWidth = [node.totalwidth for node in avgSMS.nodes if not node.isSM]
            for im,m in enumerate(avgMass):
                m = roundValue(m.asNumber(GeV),5)
                self.assertAlmostEqual(m,averageMasses[ic][im],2)
            for iw,w in enumerate(avgWidth):
                    w = roundValue(w.asNumber(GeV),5)
                    self.assertAlmostEqual(w,averageWidths[ic][iw],2)
            

        #Check clustering with distinct elements but no maxDist limit
        clusters = clusterTools.clusterSMS(smsList,maxDist=10.,dataset=dataset)
        self.assertEqual(len(clusters),1)
        cluster = clusters[0]
        self.assertEqual(sorted(cluster.smsList),sorted(smsList))
        avgSMS = cluster.averageSMS
        avgMass = [node.mass for node in avgSMS.nodes if not node.isSM]
        avgWidth = [node.totalwidth for node in avgSMS.nodes if not node.isSM]

        averageMass = (1000.*GeV*weights[0]+1020.*GeV*weights[1] + 500.*GeV*weights[2])/sum(weights)
        averageMass = [averageMass,averageMass,100.*GeV,100.*GeV]
        averageMass = [roundValue(m,5) for m in averageMass]

        averageWidths = (n2.totalwidth*weights[0]+n3.totalwidth*weights[1] + n4.totalwidth*weights[2])/sum(weights)
        averageWidths = [averageWidths,averageWidths,0*GeV,0*GeV]
        averageWidths = [roundValue(m,5)for m in averageWidths]
        for im,m in enumerate(avgMass):
            m = roundValue(m,5)
            self.assertAlmostEqual(m.asNumber(GeV),averageMass[im].asNumber(GeV),2)
        for iw,w in enumerate(avgWidth):
            w = roundValue(w,5)
            self.assertAlmostEqual(w.asNumber(GeV),averageWidths[iw].asNumber(GeV),2)


        #Check clustering where elements have same upper limits, but not the average element:
        smsList[0]._upperLimit = 1.*fb
        smsList[1]._upperLimit = 1.*fb
        smsList[2]._upperLimit = 1.*fb
        clusters = clusterTools.clusterSMS(smsList,maxDist=0.1,dataset=dataset)
        self.assertEqual(len(clusters),2)

    def testClustererLifeTimes(self):
        """ test the clustering with distinct lifetimes"""


        data = [[ [[ 674.99, 199.999], [ 674.99, 199.999] ],.03 ],
               [ [[ 725.0001,200.], [ 725.0001,200.] ], .06 ] ,
               [ [[ 750.,250.], [ 750.,250.] ], .03 ] ]

        xvalues = [np.array(pt[0]).flatten() for pt in data]
        yvalues = [pt[1] for pt in data]

        info = Info("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/data/dataInfo.txt")
        globalInfo = Info("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/globalInfo.txt")
        txnameData=TxNameData(x=xvalues, y=yvalues, txdataId=5)
        txname=TxName("./database/8TeV/ATLAS/ATLAS-SUSY-2013-05/data/T2bb.txt",
                        globalInfo,info,finalStates)
        txname.txnameData = txnameData
        dataset = DataSet(info = globalInfo, createInfo = False)
        dataset.dataInfo = info
        dataset.txnameList = [txname]

        slhafile = 'testFiles/slha/lightEWinos.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(slhafile)

        sms1 = fromString("(PV > sb_1(1),sb_1(2)), (sb_1(1) > b,N1), (sb_1(2) > b,N1)",
                         model=model)
        sms2 = fromString("(PV > sb_2(1),sb_2(2)), (sb_2(1) > b,N1), (sb_2(2) > b,N1)",
                         model=model)

        smsList = [sms1,sms2]
        n1 = model.getParticle(label='N1')
        sb1 = model.getParticle(label='sb_1')
        sb2 = model.getParticle(label='sb_2')
        n1.mass = 200*GeV
        n1.totalwidth = 0.*GeV
        sb1.mass = 675.*GeV
        sb1.totalwidth = 1e-15*GeV
        sb2.mass = 675.*GeV
        sb2.totalwidth = 0.9e-15*GeV

        w1 = XSectionList()
        w1.xSections.append(XSection())
        w1.xSections[0].info = XSectionInfo()
        w1.xSections[0].info.sqrts = 8.*TeV
        w1.xSections[0].info.label = '8 TeV'
        w1.xSections[0].info.order = 0
        w1.xSections[0].value = 10.*fb

        smsList = [txname.hasSMSas(sms) for sms in smsList]
        for sms in smsList:
            sms.weight = w1.getMaxXsec()
            sms.txname = txname
            sms.eff = 1. #(Used in clustering)

        clusters = clusterTools.clusterSMS(smsList, 5., dataset)
        self.assertEqual(len(clusters),1)
        cluster = clusters[0]
        avgSMS = cluster.averageSMS
        self.assertEqual(avgSMS.mass[1],675.*GeV)
        self.assertAlmostEqual(avgSMS.totalwidth[1].asNumber(GeV)*1e15,0.95,4)

        clusters = clusterTools.clusterSMS(smsList, .5, dataset)
        #in this example the distance is in maxdist, so we cluster
        self.assertEqual(len(clusters),1)

        clusters=clusterTools.clusterSMS(smsList, .1, dataset)
        #in this example the distance is not in maxdist, so we dont cluster
        self.assertEqual(len(clusters),2)

    def testComplexCluster(self):
        """ test the mass clusterer """

        slhafile = 'testFiles/slha/416126634.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(slhafile)
        sigmacut = 0.03*fb
        mingap = 5.*GeV
        topDict = decomposer.decompose(model, sigmacut, massCompress=True, invisibleCompress=True, minmassgap=mingap)


        #Test clustering for UL results
        database.selectExpResults(analysisIDs='CMS-SUS-16-039',dataTypes='upperLimit')
        predictions = theoryPredictionsFor(database, topDict, combinedResults=False)
        clusterSizes = sorted([len(p.smsList) for p in predictions])
        self.assertEqual(clusterSizes, [1,16,16])

    def testClusterFails(self):
        """ test a case where the average cluster SMS does not have a well defined UL"""
        from smodels.experiment.databaseObj import Database
        db = Database ( "https://smodels.github.io/database/unittest_PAS12022" )
        er = db.getExpResults ( analysisIDs='CMS-PAS-SUS-12-022', txnames='TChiWZ' )[0]
        dataset = er.datasets[0]
        txname = dataset.txnameList[0]

        slhafile = 'testFiles/clusterfails/failed_cluster.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(slhafile,ignorePromptQNumbers=['spin','eCharge','colordim'])
        sigmacut = 0.001*fb
        mingap = 5.*GeV
        smsDict = decomposer.decompose(model, sigmacut, massCompress=True, invisibleCompress=True, minmassgap=mingap)

        # Select SMS:
        smsFilter = ['(PV > C1-(1),N2(2)), (C1-(1) > N1~,W-), (N2(2) > N1,Z)', 
                     '(PV > C1-(1),N3(2)), (C1-(1) > N1~,W-), (N3(2) > N1,Z)', 
                     '(PV > N2(1),C1+(2)), (N2(1) > N1,Z), (C1+(2) > N1,W+)',                     
                     '(PV > C1+(1),N2(2)), (C1+(1) > N1,W+), (N2(2) > N1,Z)', 
                     '(PV > N3(1),C1+(2)), (N3(1) > N1,Z), (C1+(2) > N1,W+)',
                     '(PV > C1+(1),N3(2)), (C1+(1) > N1,W+), (N3(2) > N1,Z)'
                     ]
        # The SMS have to be in the correct order for clustering of all SMS to fail!
        smsList = [sms for sms in smsDict.getSMSList() if str(sms) in smsFilter]
        for sms in smsList:
            sms.weight = 1*fb
            sms.txname = txname
        smsList = sorted(smsList, key = lambda sms: smsFilter.index(str(sms)))

        clusters = clusterTools.clusterSMS(smsList, 0.2, dataset)
        self.assertEqual(len(clusters),2)


if __name__ == "__main__":
    unittest.main()
