#!/usr/bin/env python3
 
"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS
 
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
 
"""
 
import sys,os
import importlib
sys.path.insert(0,"../")
import unittest
from unitTestHelpers import equalObjs, runMain
from smodels.tools.smodelsLogging import setLogLevel
from smodels.tools import runtime
from smodels import particlesLoader
from imp import reload
import subprocess

setLogLevel('debug')
 
 
class ModelsTest(unittest.TestCase):
    definingRun = False ## meant only to adapt to changes in output format
    ## use with super great care!!
  
    def testRuntimeImport(self):
        filename = "./testFiles/slha/idm_example.slha"
        runtime.modelFile = 'idm'
        reload(particlesLoader)
        outputfile = runMain(filename,inifile='testParameters_noModel.ini',suppressStdout=True)
        if self.definingRun:
            logger.error ( "This is a definition run! Know what youre doing!" )
            default = "idm_example_defaultB.py"
            cmd = "cat %s | sed -e 's/smodelsOutput/smodelsOutputDefault/' > %s" % ( outputfile, default )
            a = subprocess.getoutput ( cmd )
        spec = importlib.util.spec_from_file_location( "output", outputfile)
        output_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(output_module)
        smodelsOutput = output_module.smodelsOutput
        from idm_example_defaultB import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'database version']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.02,
                           ignore=ignoreFields)
        self.assertTrue(equals)
        self.removeOutputs(outputfile)
         
    def testParameterFile(self):
        filename = "./testFiles/slha/idm_example.slha"
        outputfile = runMain(filename,inifile='testParameters_idm.ini',suppressStdout=True)        
        if self.definingRun:
            logger.error ( "This is a definition run! Know what youre doing!" )
            default = "idm_example_default.py"
            cmd = "cat %s | sed -e 's/smodelsOutput/smodelsOutputDefault/' > %s" % ( outputfile, default )
            a = subprocess.getoutput ( cmd )
        spec = importlib.util.spec_from_file_location( "output", outputfile)
        output_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(output_module)
        smodelsOutput = output_module.smodelsOutput
        from idm_example_default import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'database version']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.02,
                           ignore=ignoreFields)
        self.assertTrue(equals)
        self.removeOutputs(outputfile)
        
    def testWrongModel(self):
        runtime.modelFile = 'mssm'
        reload(particlesLoader)
        filename = "./testFiles/slha/idm_example.slha"
        outputfile = runMain(filename,suppressStdout=True)
        spec = importlib.util.spec_from_file_location( "output", outputfile)
        output_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(output_module)
        smodelsOutput = output_module.smodelsOutput
        self.assertTrue(smodelsOutput['OutputStatus']['decomposition status'] < 0)  
        self.removeOutputs(outputfile)


    def removeOutputs(self, f):
        """ remove cruft outputfiles """
        for i in [ f, f.replace(".py",".pyc") ]:
            if os.path.exists ( i ): os.remove ( i )
  
      

if __name__ == "__main__":
    unittest.main()
