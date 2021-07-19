#!/usr/bin/env python3
 
"""
.. module:: testInteractivePlots
   :synopsis: Tests the interactive plot tool
 
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
 
"""
 
import sys,os,shutil
sys.path.insert(0,"../")
import unittest
from smodels.tools.interactivePlots import main
try:
    from types import SimpleNamespace
except ImportError: ## doesnt exist in python2
    class SimpleNamespace:
        pass


class RunInteractivePlotSTest(unittest.TestCase):
    
    # this test corresponds to calling
    # ../smodelsTools.py interactive-plots -f ./testFiles/scanExample/smodels-output.tar.gz -s testFiles/scanExample/slhas.tar.gz -p iplots_parameters.py
    def testInteractivePlots(self):
        slhaFolder = './testFiles/scanExample/slhas.tar.gz'
        smodelsFolder = './testFiles/scanExample/smodels-output.tar.gz'
        parametersFile = './iplots_parameters.py'
        outFolder = './plots_test'
        
        defaultFolder = './plots_test_default'
        
        if os.path.isdir(outFolder):
            shutil.rmtree(outFolder)

        parser = SimpleNamespace()
        parser.smodelsFolder = smodelsFolder
        parser.slhaFolder  =  slhaFolder
        parser.parameters = parametersFile
        parser.outputFolder = outFolder
        parser.verbosity = 'error'
        parser.npoints = -1
        parser.modelFile='../smodels/share/models/mssm.py'

        run = main(parser)
        
        self.assertEqual(run,outFolder)        
        self.assertEqual(sorted(os.listdir(outFolder)), sorted(os.listdir(defaultFolder)))

        if os.path.isdir(outFolder):
            shutil.rmtree(outFolder)




    # this test corresponds to calling
    # ../smodelsTools.py interactive-plots -f ./testFiles/scanExampleIDM/smodels-output.tar.gz -s testFiles/scanExampleIDM/slhas.tar.gz -p iplots_parameters_IDM.py

    def testInteractivePlotsIDM(self):
        slhaFolder = './testFiles/scanExampleIDM/slhas.tar.gz'
        smodelsFolder = './testFiles/scanExampleIDM/smodels-output.tar.gz'
        parametersFile = './iplots_parameters_IDM.py'
        outFolder = './plots_test_idm'
        
        defaultFolder = './plots_test_default_idm'
        
        if os.path.isdir(outFolder):
            shutil.rmtree(outFolder)

        parser = SimpleNamespace()
        parser.smodelsFolder = smodelsFolder
        parser.slhaFolder  =  slhaFolder
        parser.parameters = parametersFile
        parser.outputFolder = outFolder
        parser.verbosity = 'error'
        parser.npoints = -1
        parser.modelFile='../smodels/share/models/idm.py'

        run = main(parser)
        
        self.assertEqual(run,outFolder)
        self.assertEqual(sorted(os.listdir(outFolder)), sorted(os.listdir(defaultFolder)))

        if os.path.isdir(outFolder):
            shutil.rmtree(outFolder)



if __name__ == "__main__":
    unittest.main()
