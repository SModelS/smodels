#!/usr/bin/env python

"""
.. module:: testFlongCalc
   :synopsis: Tests the use of internal and external Flong calculators
    
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
    
"""
import sys
sys.path.insert(0,"../")
from smodels.tools import runtime
from smodels.installation import installDirectory
from smodels.theory.slhaDecomposer import  _getPromptDecays,_getDictionariesFromSLHA
import unittest
# from smodels.tools.smodelsLogging import setLogLevel
# 
# setLogLevel("info")


def myFlong(pdg,mass,width):
    
    return {'Fprompt' : 0.1, 'Flong' : 0.5}


class TxTest(unittest.TestCase):

    def testFlongDefault(self):
        runtime.setFlongCalc(path= "%ssmodels/tools/flongCalc.py" %(installDirectory()),
                                method = "FlongCalculator")
        slhafile = "%sinputFiles/slha/longLived.slha" %(installDirectory())
        # Get BRs and masses from file
        brDic,_ = _getDictionariesFromSLHA(slhafile)
        
        
        brDic = _getPromptDecays(slhafile,brDic)
        stauPrompt = 0.
        stauLong = 0.
        for decay in brDic[1000015]:
            if decay.ids:
                stauPrompt += decay.br
            else:
                stauLong += decay.br

        self.assertAlmostEqual(stauLong,0.9916,4)
        self.assertAlmostEqual(stauPrompt,5.068e-7,4)
        
        runtime.setFlongCalc()


    def testFlongFallback(self):
        runtime.setFlongCalc(path= "Wrong Path", method = "FlongCalculator")
        slhafile = "%sinputFiles/slha/longLived.slha" %(installDirectory())
        # Get BRs and masses from file
        brDic,_ = _getDictionariesFromSLHA(slhafile)
        
        
        brDic = _getPromptDecays(slhafile,brDic)
        stauPrompt = 0.
        stauLong = 0.
        for decay in brDic[1000015]:
            if decay.ids:
                stauPrompt += decay.br
            else:
                stauLong += decay.br

        self.assertAlmostEqual(stauLong,0.9916,4)
        self.assertAlmostEqual(stauPrompt,5.068e-7,4)

        runtime.setFlongCalc()

    def testFlongNew(self):
        
        
        runtime.setFlongCalc(path="./testFlongCalc.py",method="myFlong")

        slhafile = "%sinputFiles/slha/longLived.slha" %(installDirectory())
        # Get BRs and masses from file
        brDic,_ = _getDictionariesFromSLHA(slhafile)
        
        
        brDic = _getPromptDecays(slhafile,brDic)
        stauPrompt = 0.
        stauLong = 0.
        for decay in brDic[1000015]:
            if decay.ids:
                stauPrompt += decay.br
            else:
                stauLong += decay.br

        self.assertAlmostEqual(stauLong,0.5,4)
        self.assertAlmostEqual(stauPrompt,0.1,4)
        
        runtime.setFlongCalc()
        
        

if __name__ == "__main__":
    unittest.main()
