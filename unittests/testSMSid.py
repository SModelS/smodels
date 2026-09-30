#!/usr/bin/env python3

"""
.. module:: testElementID
   :synopsis: ??

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.decomposition import decomposer
from smodels.base.physicsUnits import GeV,fb
from databaseLoader import database
from smodels.matching.theoryPrediction import theoryPredictionsFor
from smodels.share.models.SMparticles import SMList
from smodels.share.models.mssm import BSMList
from smodels.base.model import Model

class SMSIdTest(unittest.TestCase):
    def testGoodFile(self):

        listOfIDs = {'ATLAS-CONF-2013-037':  [31+8, 26+8, 33+8, 34+8, 27+8, 28+8, 30+8, 32+8],
                     'ATLAS-SUSY-2013-05' : [29+8]}
        filename = "./testFiles/slha/higgsinoStop.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(filename,promptWidth = 1e-12*GeV,
                              ignorePromptQNumbers=['spin','eCharge','colordim']) #Make sure C1 and N2 are treated as prompt
        ## model.describe()

        topDict = decomposer.decompose(model, sigmacut= 0.1*fb, massCompress=True, invisibleCompress=True, 
                                       minmassgap= 5*GeV, minmassgapISR= 5*GeV)
        database.selectExpResults(analysisIDs=['*:8*TeV','CMS-PAS-SUS-15-002','CMS-PAS-SUS-16-024'])
        allTheorypredictions = [tp for tp in theoryPredictionsFor(database, topDict) if 'EXO' not in tp.analysisId()]
        tpDict = {tp.analysisId() : [] for tp in allTheorypredictions}
        for tp in allTheorypredictions:
            tpDict[tp.analysisId()].append(tp)
                
        for anaID,tpList in tpDict.items():
            self.assertEqual(len(tpList),1)
            if not tpList:
                continue
            tp = tpList[0]
            tpIDs = [sms.smsID for sms in tp.smsList]
            self.assertEqual(sorted(tpIDs),sorted(listOfIDs[anaID]))            

if __name__ == "__main__":
    unittest.main()
