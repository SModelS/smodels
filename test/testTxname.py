#!/usr/bin/env python

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
        tx = TxName(f,gInfo,gInfo,True)
        
        el = Element(info="[[*],[]]",finalState = ['MET','HSCP'])
        
        self.assertTrue(len(tx._topologyList.getElements()), 1)
        self.assertEqual(tx._topologyList.getElements()[0], el)
        
        
    def testBrokenFinalState(self):
        
        f = './databaseBroken/13TeV/CMS/CMS-PAS-EXO-16-036-eff/c000/THSCPM2broken.txt'
        gInfo = infoObj.Info('./databaseBroken/13TeV/CMS/CMS-PAS-EXO-16-036-eff/globalInfo.txt')
        gInfo.addInfo('dataId','c000')
        try:
            print (gInfo)
            tx = TxName(f,gInfo,gInfo,True)
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
        tx = TxName(f,gInfo,gInfo,True)                        
        
        self.assertFalse(hasattr(tx._topologyList.getElements()[0], 'mass'))
        
        el = Element(info="[[[e+]],[]]",finalState = ['HSCP','MET'])
        setattr(n1, 'mass', 200*GeV)
        setattr(c1, 'mass', 150*GeV)
        el.branches[0].BSMparticles = [n1]
        el.branches[1].BSMparticles = [c1]        
        
        self.assertEqual(tx.getEfficiencyFor(el), 0.12396)
         



if __name__ == "__main__":
    unittest.main()
