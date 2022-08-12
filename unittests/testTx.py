#!/usr/bin/env python3

"""
.. module:: testTx
   :synopsis: Tests with Tx slha input files.
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""
import sys
sys.path.insert(0,"../")
from smodels.decomposition import decomposer
from smodels.base.model import Model
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.physicsUnits import GeV, fb
import unittest

class TxTest(unittest.TestCase):

    def testT1(self):

        """ test with the T1 slha input file """
        slhafile="./testFiles/slha/simplyGluino.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(inputFile=slhafile)
        topos = decomposer.decompose(model, .1*fb, False, False, 5.*GeV)
        self.assertEqual(len(topos),1)
        self.assertEqual(len(topos.getSMSList()),1)
        sms = topos.getSMSList()[0]
        masses = [node.mass for node in sms.nodes if node.isSM is False]
        self.assertEqual(masses,[675.*GeV,675.*GeV,200.*GeV,200.*GeV])
        evenParticles = sms.treeToBrackets()[0]
        evenParticles = str(evenParticles).replace("'","").replace(' ', '')
        self.assertEqual(evenParticles, "[[[q,q]],[[q,q]]]" )

if __name__ == "__main__":
    unittest.main()
