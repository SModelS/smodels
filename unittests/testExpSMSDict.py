#!/usr/bin/env python3

"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys

sys.path.insert(0,"../")
import unittest

from databaseLoader import database

from smodels.experiment.databaseObj import Database


class ExpSMSDictTest(unittest.TestCase):

    def testDict(self):
        
        database.selectExpResults(useNonValidated=True)
        expDict = database.expSMSDict
        self.assertEqual(len(list(expDict.getSMS())),80)
        self.assertEqual(len(list(expDict.getTx())),644)

    def testFilter(self):

        database.selectExpResults(analysisIDs='CMS*',useNonValidated=True)
        expDict = database.expSMSDict

        self.assertEqual(len(list(expDict.getSMS())),36)
        self.assertEqual(len(list(expDict.getTx())),439)

        # Reset database:
        database.selectExpResults(useNonValidated=True)
        expDict = database.expSMSDict
        self.assertEqual(len(list(expDict.getSMS())),80)
        self.assertEqual(len(list(expDict.getTx())),644)

    def testNodeDict(self):

        database.selectExpResults(txnames='TChiWH',analysisIDs='ATLAS-SUSY-2013-12')
        expDict = database.expSMSDict
        self.assertEqual(len(list(expDict.getSMS())),1)
        self.assertEqual(len(list(expDict.getTx())),1)

        tx = list(expDict.getTx())[0]
        txSMS = list(tx.smsMap.keys())[0]
        uniqueSMS = list(expDict.getSMS())[0]
        nodesDict = expDict._nodesDict[tx]['sms_1']

        # Test if nodeDict is correctly translating the nodes
        for node in uniqueSMS.nodeIndices:
            nA = uniqueSMS.indexToNode(node)
            nB = txSMS.indexToNode(nodesDict[node])
            self.assertEqual(str(nA),str(nB))

        # Test if relabeling is working
        sms = uniqueSMS.copy()
        smsR = expDict.setTxNodeOrdering(sms,tx,'sms_1')
        for nodeIndex,nA in zip(smsR.nodeIndices,smsR.nodes):
            nB = txSMS.indexToNode(nodeIndex)
            self.assertEqual(str(nA),str(nB))
        

if __name__ == "__main__":
    unittest.main()
