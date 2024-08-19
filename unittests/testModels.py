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
from smodels.base.smodelsLogging import setLogLevel
from smodels.base import runtime
from smodels.tools.particlesLoader import load
import subprocess
from smodels.base.smodelsLogging import logger
from smodels.base.model import Model
from smodels.share.models.SMparticles import SMList

setLogLevel('error')


class ModelsTest(unittest.TestCase):
    definingRun = False ## meant only to adapt to changes in output format
    ## use with super great care!!

    def testRuntimeImport(self):
        filename = "./testFiles/slha/idm_example.slha"
        runtime.modelFile = 'idm'
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
                    'Total xsec for topologies outside the grid (fb)',  'model', 'promptwidth', 'stablewidth', 'checkinput',
                    'doinvisible', 'docompress', 'computestatistics', 'testcoverage', 'combinesrs']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedRelDiff=0.1,
                           ignore=ignoreFields, fname = outputfile )
        self.assertTrue(equals)
        self.removeOutputs(outputfile)
        self.removeOutputs("./idm_example.slha" )

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
                            'Total xsec for topologies outside the grid (fb)',  'model', 'promptwidth', 'stablewidth', 'checkinput',
                            'doinvisible', 'docompress', 'computestatistics', 'testcoverage', 'combinesrs']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedRelDiff=0.1,
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
                            'Total xsec for topologies outside the grid (fb)',  'model', 'promptwidth', 'stablewidth', 'checkinput',
                            'doinvisible', 'docompress', 'computestatistics', 'testcoverage', 'combinesrs']
        smodelsOutputDefault['ExptRes'] = sorted(smodelsOutputDefault['ExptRes'],
                    key=lambda res: res['r'], reverse=True)
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedRelDiff=0.1,
                           ignore=ignoreFields, fname = outputfile )

        self.assertTrue(equals)
        self.removeOutputs(outputfile)

    def testParticlesFromSLHA(self):

        from smodels.tools.particlesLoader import getParticlesFromSLHA
        from smodels.base.particle import Particle

        BSMList = getParticlesFromSLHA("./testFiles/slha/idm_example.slha")
        BSMList = sorted(BSMList, key = lambda p: p.pdg)

        h2 = Particle(isSM=False,label='h2',pdg=35,eCharge=0.0,colordim=1,spin=0.0)
        h3 = Particle(isSM=False,label='h3',pdg=36,eCharge=0.0,colordim=1,spin=0.0)
        hp = Particle(isSM=False,label='h+',pdg=37,eCharge=1.0,colordim=1,spin=0.0)
        hm = Particle(isSM=False,label='h-',pdg=-37,eCharge=-1.0,colordim=1,spin=0.0)
        BSMListDefault = [hm,h2,h3,hp]

        self.assertEqual(len(BSMList),len(BSMListDefault))
        for i,p in enumerate(BSMListDefault):
            pDict = p.__dict__
            pDict.pop('_comp')
            pDict.pop('_id')
            pBDict = BSMList[i].__dict__
            pBDict.pop('_comp')
            pBDict.pop('_id')
            self.assertEqual(pDict,pBDict)

    def testWrongModel(self):
        filename = "./testFiles/slha/simplyGluino.slha"
        outputfile = runMain(filename,inifile='testParameters_idm.ini',suppressStdout=True)
        smodelsOutput = importModule ( outputfile )
        self.assertTrue(smodelsOutput['OutputStatus']['decomposition status'] < 0)
        self.removeOutputs(outputfile)

    def testUpdateParameters(self):
        from smodels.share.models.mssm import BSMList
        from smodels.base.physicsUnits import GeV

        slhafile = "./testFiles/slha/hscpTest_short.slha"
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile,minMass=1.0*GeV)
        p = model.getParticlesWith(pdg=1000012)[0]
        self.assertEqual(p.mass,1.0*GeV)

        model.updateParticles(inputFile=slhafile,minMass=5.0*GeV)
        p = model.getParticlesWith(pdg=1000012)[0]
        self.assertEqual(p.mass,5.0*GeV)

    # Check if models loaded through python module or through SLHA QNUMBERS give the same model
    def testRuntimeImport(self):

        filename = "./testFiles/slha/gluino_squarks.slha"
        runtime.modelFile = 'mssm'
        BSMList = load()
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(filename)

        runtime.modelFile = 'mssmQNumbers.slha'
        BSMList = load()
        modelB = Model(BSMparticles=BSMList, SMparticles=SMList)
        modelB.updateParticles(filename)

        for ptc in model.BSMparticles:
            ptcB = modelB.getParticlesWith(pdg = ptc.pdg)
            if not ptcB: #If particule is its own anti-particle, it should not appear in modelB
                ptcB = modelB.getParticlesWith(pdg = -ptc.pdg)
            self.assertEqual(len(ptcB),1)
            ptcB = ptcB[0]
            for attr,val in ptc.__dict__.items():
                if attr in ['_id','label','_comp','pdg','_isInvisible']:
                    continue
                valB = getattr(ptcB,attr)
                if attr == 'decays':
                    self.assertEqual(len(val),len(valB))
                    continue
                self.assertEqual(val,valB)


    def removeOutputs(self, f):
        """ remove cruft outputfiles """
        for i in [ f, f.replace(".py",".pyc" ), f.replace(".py",".smodels") ]:
            if os.path.exists ( i ): os.remove ( i )



if __name__ == "__main__":
    unittest.main()
