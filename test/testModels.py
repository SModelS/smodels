#!/usr/bin/env python3

"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys,os
sys.path.insert(0,"../")
import unittest
from unitTestHelpers import equalObjs, runMain, importModule
from smodels.tools.smodelsLogging import setLogLevel
from smodels.tools import runtime
from smodels import particlesLoader
from imp import reload
import subprocess
from smodels.tools.smodelsLogging import logger
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
        smodelsOutput = importModule ( outputfile )
        from idm_example_defaultB import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'Element',
                    'database version', 'Total missed xsec',
                    'Missed xsec long-lived', 'Missed xsec displaced',
                    'Missed xsec MET', 'Total outside grid xsec',
                    'Total xsec for missing topologies (fb)',
                    'Total xsec for missing topologies with displaced decays (fb)',
                    'Total xsec for missing topologies with prompt decays (fb)',
                    'Total xsec for topologies outside the grid (fb)']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.1,
                           ignore=ignoreFields, fname = outputfile )
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
        smodelsOutput = importModule ( outputfile )
        from idm_example_default import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'Element', 'database version', 'Total missed xsec',
                            'Missed xsec long-lived', 'Missed xsec displaced', 'Missed xsec MET', 'Total outside grid xsec',
                            'Total xsec for missing topologies (fb)','Total xsec for missing topologies with displaced decays (fb)',
                            'Total xsec for missing topologies with prompt decays (fb)',
                            'Total xsec for topologies outside the grid (fb)']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.1,
                           ignore=ignoreFields, fname = outputfile )
        self.assertTrue(equals)
        self.removeOutputs(outputfile)

    def testModelFromSLHA(self):
        filename = "./testFiles/slha/idm_example.slha"
        #Test the case where the BSM particles are defined by the SLHA file:
        outputfile = runMain(filename,inifile='testParameters_idmB.ini',suppressStdout=True)
        smodelsOutput = importModule ( outputfile )
        from idm_example_default import smodelsOutputDefault
        ignoreFields = ['input file','smodels version', 'ncpus', 'Element', 'database version', 'Total missed xsec',
                            'Missed xsec long-lived', 'Missed xsec displaced', 'Missed xsec MET', 'Total outside grid xsec',
                            'Total xsec for missing topologies (fb)','Total xsec for missing topologies with displaced decays (fb)',
                            'Total xsec for missing topologies with prompt decays (fb)',
                            'Total xsec for topologies outside the grid (fb)']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.1,
                           ignore=ignoreFields, fname = outputfile )

        self.assertTrue(equals)
        self.removeOutputs(outputfile)

    def testParticlesFromSLHA(self):

        from smodels.particlesLoader import getParticlesFromSLHA
        from smodels.theory.particle import Particle

        BSMList = getParticlesFromSLHA("./testFiles/slha/idm_example.slha")
        BSMList = sorted(BSMList, key = lambda p: p.pdg)

        h2 = Particle(Z2parity=-1,label='h2',pdg=35,eCharge=0,colordim=1,spin=0)
        h3 = Particle(Z2parity=-1,label='h3',pdg=36,eCharge=0,colordim=1,spin=0)
        hp = Particle(Z2parity=-1,label='h+',pdg=37,eCharge=1,colordim=1,spin=0)
        hm = Particle(Z2parity=-1,label='h-',pdg=-37,eCharge=-1,colordim=1,spin=0)
        BSMListDefault = [hm,h2,h3,hp]
        for i,p in enumerate(BSMListDefault):
            self.assertEqual(p.__dict__,BSMList[i].__dict__)

    def testWrongModel(self):
        runtime.modelFile = 'mssm'
        reload(particlesLoader)
        filename = "./testFiles/slha/idm_example.slha"
        outputfile = runMain(filename,suppressStdout=True)
        smodelsOutput = importModule ( outputfile )
        self.assertTrue(smodelsOutput['OutputStatus']['decomposition status'] < 0)
        self.removeOutputs(outputfile)

    def removeOutputs(self, f):
        """ remove cruft outputfiles """
        for i in [ f, f.replace(".py",".pyc" ), f.replace(".py",".smodels") ]:
            if os.path.exists ( i ): os.remove ( i )



if __name__ == "__main__":
    unittest.main()
