#!/usr/bin/env python3

"""
.. module:: testSplitSMSstr
   :synopsis: Test the splitting of SMS process strings

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from __future__ import print_function
import sys
sys.path.insert(0,"../")
import unittest
from smodels.experiment.expAuxiliaryFuncs import smsInStr, bracketToProcessStr
from unitTestHelpers import getNodesIndices, getEdges



class TestSplitStr(unittest.TestCase):
    
    
    def test_split(self):
        cString = "71.*([[['mu+','mu-']],[['l','nu']]] + [[['e+','e-']],[['l','nu']]])"
        procStr = "71.*({(PV > anyBSM(1),anyBSM(2)),(anyBSM(1) > mu+,mu-,MET),(anyBSM(2) > l,nu,MET)} + {(PV > anyBSM(1),anyBSM(2)),(anyBSM(1) > e+,e-,MET),(anyBSM(2) > l,nu,MET)})"

        cStr2 = "cGtr([[['mu+','mu-']],[['l','nu']]],[[['e+','e-']],[['l','nu']]])"
        procStr2 = "cGtr({(PV > anyBSM(1),anyBSM(2)),(anyBSM(1) > mu+,mu-,MET),(anyBSM(2) > l,nu,MET)},{(PV > anyBSM(1),anyBSM(2)),(anyBSM(1) > e+,e-,MET),(anyBSM(2) > l,nu,MET)})"


        procStr3 = "{(PV > anyBSM(1),anyBSM(2)),(anyBSM(1) > mu+,mu-,MET),(anyBSM(2) > l,nu,MET)}"

        smsSplit = smsInStr(cString)
        self.assertEqual(len(smsSplit),2)
        self.assertEqual(smsSplit,["[[['mu+','mu-']],[['l','nu']]]", 
                                   "[[['e+','e-']],[['l','nu']]]"])
        
        smsSplit = smsInStr(procStr)
        self.assertEqual(len(smsSplit),2)
        self.assertEqual(smsSplit,['(PV>anyBSM(1),anyBSM(2)),(anyBSM(1)>mu+,mu-,MET),(anyBSM(2)>l,nu,MET)', 
                                   '(PV>anyBSM(1),anyBSM(2)),(anyBSM(1)>e+,e-,MET),(anyBSM(2)>l,nu,MET)'])
        
        smsSplit = smsInStr(cStr2)
        self.assertEqual(len(smsSplit),2)
        self.assertEqual(smsSplit,["[[['mu+','mu-']],[['l','nu']]]", 
                                   "[[['e+','e-']],[['l','nu']]]"])
        
        smsSplit = smsInStr(procStr2)
        self.assertEqual(len(smsSplit),2)
        self.assertEqual(smsSplit,['(PV>anyBSM(1),anyBSM(2)),(anyBSM(1)>mu+,mu-,MET),(anyBSM(2)>l,nu,MET)', 
                                   '(PV>anyBSM(1),anyBSM(2)),(anyBSM(1)>e+,e-,MET),(anyBSM(2)>l,nu,MET)'])
        
        smsSplit = smsInStr(procStr3)
        self.assertEqual(len(smsSplit),1)
        self.assertEqual(smsSplit,['(PV>anyBSM(1),anyBSM(2)),(anyBSM(1)>mu+,mu-,MET),(anyBSM(2)>l,nu,MET)'])

        cStr4 = "cGtr([[[mu+,mu-]],[['l',nu]]],[[['e+',e+]],[['l',nu]]])"
        smsSplit = smsInStr(cStr4)
        self.assertEqual(len(smsSplit),2)
        self.assertEqual(smsSplit,["[[['mu+','mu-']],[['l','nu']]]", 
                                   "[[['e+','e+']],[['l','nu']]]"])
        
        cStr5 = "[[[jet,jet]],[[jet,jet]]]"
        smsSplit = smsInStr(cStr5)
        self.assertEqual(len(smsSplit),1)
        self.assertEqual(smsSplit,["[[['jet','jet']],[['jet','jet']]]"])


        cStr6 = "[[['*']],[]]+[[['*','*']],[]]+[[['*','*','*']],[]]"
        smsSplit = smsInStr(cStr6)
        self.assertEqual(len(smsSplit),3)
        self.assertEqual(smsSplit,["[[['*']],[]]", 
                                   "[[['*','*']],[]]", 
                                   "[[['*','*','*']],[]]"])
        
        cStr7 = "[[[jet],[Z]],[[jet],[W]]]+[[[jet],[higgs]],[[jet],[W]]]+[[[jet],[Z]],[[jet],[Z]]]+[[[jet],[higgs]],[[jet],[higgs]]]+[[[jet],[Z]],[[jet],[higgs]]]"
        smsSplit = smsInStr(cStr7)
        self.assertEqual(len(smsSplit),5)
        self.assertEqual(smsSplit,["[[['jet'],['Z']],[['jet'],['W']]]", 
                                   "[[['jet'],['higgs']],[['jet'],['W']]]", 
                                   "[[['jet'],['Z']],[['jet'],['Z']]]", 
                                   "[[['jet'],['higgs']],[['jet'],['higgs']]]", 
                                   "[[['jet'],['Z']],[['jet'],['higgs']]]"])
        
        cStr8 = "[[[jet],[Z]],[[jet],[W]]]+[[[jet],[higgs]],[[jet],[W]]]+[[[jet],[Z]],[[jet],[Z]]]"
        smsSplit = smsInStr(cStr8)
        self.assertEqual(len(smsSplit),3)
        self.assertEqual(smsSplit,["[[['jet'],['Z']],[['jet'],['W']]]", 
                                   "[[['jet'],['higgs']],[['jet'],['W']]]", 
                                   "[[['jet'],['Z']],[['jet'],['Z']]]"])
        
        cStr9 = "2.*([[['L+'],['L-']],[['L'],['nu']]] + [[['L+'],['L-']],[['nu'],['L']]] + [[['L-'],['L+']],[['L'],['nu']]] + [[['L-'],['L+']],[['nu'],['L']]])"
        smsSplit = smsInStr(cStr9)
        self.assertEqual(len(smsSplit),4)
        self.assertEqual(smsSplit,["[[['L+'],['L-']],[['L'],['nu']]]", 
                                   "[[['L+'],['L-']],[['nu'],['L']]]", 
                                   "[[['L-'],['L+']],[['L'],['nu']]]", 
                                   "[[['L-'],['L+']],[['nu'],['L']]]"])


if __name__ == "__main__":
    unittest.main()                         