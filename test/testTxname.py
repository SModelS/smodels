#!/usr/bin/env python3

"""
.. module:: testTx
   :synopsis: Tests with Tx slha input files.
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""
import sys
sys.path.insert(0,"../")
from smodels.share.models import mssm
from smodels.theory.element import Element
from smodels.experiment import infoObj
from smodels.experiment.txnameObj import TxName
from smodels.tools.physicsUnits import GeV
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
import unittest

c1 = mssm.c1
n1 = mssm.n1

class TxTest(unittest.TestCase):

    def testTxnameElements(self):
        
        f = './databaseBroken/13TeV/CMS/CMS-PAS-EXO-16-036-eff/c000/THSCPM2.txt'
        gInfo = infoObj.Info('./databaseBroken/13TeV/CMS/CMS-PAS-EXO-16-036-eff/globalInfo.txt')
        gInfo.addInfo('dataId','c000')
        tx = TxName(f,gInfo,gInfo)
        
        el = Element(info="[[*],[]]",finalState = ['MET','HSCP'])
        
        self.assertTrue(len(tx._topologyList.getElements()), 1)
        self.assertEqual(tx._topologyList.getElements()[0], el)
        
        
    def testBrokenFinalState(self):
        
        f = './databaseBroken/13TeV/CMS/CMS-PAS-EXO-16-036-eff/c000/THSCPM2broken.txt'
        gInfo = infoObj.Info('./databaseBroken/13TeV/CMS/CMS-PAS-EXO-16-036-eff/globalInfo.txt')
        gInfo.addInfo('dataId','c000')
        try:
            print (gInfo)
            tx = TxName(f,gInfo,gInfo)
            print (tx)
            gotError = False
        except SModelSError as e:         
            gotError = e
        
        errstr = 'Final state non-MET has not been defined in finalStateParticles.py'
        self.assertEqual(gotError.args[0], errstr)
        
        
    def testgetEffFor(self):
        
        f = './databaseBroken/13TeV/CMS/CMS-PAS-EXO-16-036-eff/c000/THSCPM2.txt'
        gInfo = infoObj.Info('./databaseBroken/13TeV/CMS/CMS-PAS-EXO-16-036-eff/globalInfo.txt')
        gInfo.addInfo('dataId','c000')
        tx = TxName(f,gInfo,gInfo)                        
        
        self.assertFalse(hasattr(tx._topologyList.getElements()[0], 'mass'))
        
        el = Element(info="[[[e+]],[]]",finalState = ['HSCP','MET'])
        setattr(n1, 'mass', 200*GeV)
        setattr(c1, 'mass', 150*GeV)
        el.branches[0].oddParticles = [n1]
        el.branches[1].oddParticles = [c1]        
        
        #test getting efficiency with mass only
        self.assertEqual(tx.getEfficiencyFor(el.mass), 0.12396)

        #test getting efficiency with element and reweighted efficiency
        setattr(n1, 'totalwidth', 0.*GeV)
        setattr(c1, 'totalwidth', 10**(-17)*GeV)
        self.assertAlmostEqual(tx.getEfficiencyFor(el), 0.08697,5)


    def testGetValueFor(self):

        f = './database/13TeV/CMS/CMS-PAS-EXO-16-036-eff/c000/THSCPM1.txt'
        gInfo = infoObj.Info('./database/13TeV/CMS/CMS-PAS-EXO-16-036-eff/globalInfo.txt')
        gInfo.addInfo('dataId','c000')
        tx = TxName(f,gInfo,gInfo)

        el = Element(info="[[],[]]",finalState = ['HSCP','HSCP'])
        setattr(c1, 'mass', 150*GeV)
        el.branches[0].oddParticles = [c1]
        el.branches[1].oddParticles = [c1]

        #test getting UL with mass only
        self.assertEqual(tx.txnameData.getValueFor(el.mass), 0.21496)

        #test getting UL with element and reweighted efficiency
        setattr(c1, 'totalwidth', 10**(-17)*GeV)
        self.assertAlmostEqual(tx.txnameData.getValueFor(el),0.49*0.21496,3)

      
         



if __name__ == "__main__":
    unittest.main()
