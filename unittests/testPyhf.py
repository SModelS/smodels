#!/usr/bin/env python3

"""
.. module:: testPyhf
   :synopsis: Test the pyhfInterface module

.. moduleauthor:: Gael Alguero <gael.alguero@lpsc.in2p3.fr>

"""

import sys
sys.path.insert(0,"../")
import unittest
import json
import jsonpatch
from smodels.statistics.pyhfInterface import PyhfData, PyhfUpperLimitComputer, pyhf
from smodels.base.smodelsLogging import logger
import warnings

class PyhfTest(unittest.TestCase):

    def simpleJson(self, bkg, obs):
        """
        Define a simple likelihood model under the json format
        :param bkg: list of bkg numbers
        :param obs: list of ebserved numbers

        :return: a simple likelihood specification under the json dictionary
        """
        #Defining the channels
        modifiers = []
        modifiers.append(dict(data=None,
                              type='lumi',
                              name='lumi'))
        samples = [dict(name='bkg',
                        data=bkg,
                        modifiers=modifiers)]
        channels = [dict(name='SR1',
                         samples=samples)]
        # Defining the measurements
        config = dict(poi='mu_SIG',
                      parameters=[dict(auxdata=[1],
                                       bounds=[[0.915, 1.085]],
                                       inits=[1],
                                       sigmas=[0.017],
                                       name='lumi')])
        measurements = [dict(name='BasicMeasurement',
                             config=config)]
        # Defining the observations
        observations = [dict(name='SR1',
                             data=obs)]
        ws = dict(channels=channels,
                  measurements=measurements,
                  observations=observations,
                  version='1.0.0')
        return ws

    def testCorruptJson1Signal(self):
        """
        Tests how the module handles corrupted json files
        Maybe it is needed to test different types of corruptions
        """
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        from smodels.base.smodelsLogging import setLogLevel
        setLogLevel("fatal")
        #Defining the channels
        modifiers = []
        modifiers.append(dict(data=None,
                              type='lumi',
                              name='lumi'))
        samples = [dict(name='bkg',
                        data=[10],
                        modifiers=modifiers)]
        channels = [dict(name='SR1',
                         samples=samples)]
        # Defining the measurements
        config = dict(poi='mu_SIG',
                      parameters=[dict(auxdata=[1],
                                       bounds=[[0.915, 1.085]],
                                       inits=[1],
                                       sigmas=[0.017],
                                       name='lumi')])
        measurements = [dict(name='BasicMeasurement',
                             config=config)]
        # Defining the observations
        observations = [dict(name='SR1',
                             data=[0.9])]
        # Missing channels
        ws = dict(#channels=channels,
                  measurements=measurements,
                  observations=observations,
                  version='1.0.0'
                  )
        signal = { "dummy0": { "SR1": 0.1 } }
        data = PyhfData( signal, [ws] )
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.getUpperLimitOnMu()
        self.assertEqual(ulcomputer.workspaces, None)
        self.assertEqual(ul, None)
        # Missing measurements
        ws = dict(channels=channels,
                  #measurements=measurements,
                  observations=observations,
                  version='1.0.0'
                  )
        data = PyhfData( signal, [ws] )
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.getUpperLimitOnMu()
        self.assertEqual(ulcomputer.workspaces, None)
        self.assertEqual(ul, None)
        # Missing observations
        ws = dict(channels=channels,
                  measurements=measurements,
                  #observations=observations,
                  version='1.0.0'
                  )
        data = PyhfData( signal, [ws])
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.getUpperLimitOnMu()
        self.assertEqual(ulcomputer.workspaces, None)
        self.assertEqual(ul, None)
        # Missing version
        ws = dict(channels=channels,
                  measurements=measurements,
                  observations=observations,
                  #version='1.0.0'
                  )
        data = PyhfData(signal, [ws])
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.getUpperLimitOnMu()
        self.assertIsNone(ulcomputer.workspaces)
        self.assertIsNone(ul)

    def testCorruptJson2Signal(self):
        """
        Tests how the module handles corrupted json files
        Maybe it is needed to test different types of corruptions
        """
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        #Defining the channels
        modifiers = []
        modifiers.append(dict(data=None,
                              type='lumi',
                              name='lumi'))
        samples = [dict(name='bkg',
                        data=[10, 9],
                        modifiers=modifiers)]
        channels = [dict(name='SR1',
                         samples=samples)]
        # Defining the measurements
        config = dict(poi='mu_SIG',
                      parameters=[dict(auxdata=[1],
                                       bounds=[[0.915, 1.085]],
                                       inits=[1],
                                       sigmas=[0.017],
                                       name='lumi')])
        measurements = [dict(name='BasicMeasurement',
                             config=config)]
        # Defining the observations
        observations = [dict(name='SR1',
                             data=[0.9, 0.8])]
        # Missing channels
        ws = dict(#channels=channels,
                  measurements=measurements,
                  observations=observations,
                  version='1.0.0'
                  ) 
        signal = { "dummy0": { "SR1": 0.1, "SR2": 0.2 } }
        data = PyhfData(signal, [ws])
        #data = PyhfData([[0.1, 0.2]], [ws])
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.getUpperLimitOnMu()
        self.assertEqual(ulcomputer.workspaces, None)
        self.assertEqual(ul, None)
        # Missing measurements
        ws = dict(channels=channels,
                  #measurements=measurements,
                  observations=observations,
                  version='1.0.0'
                  )
        data = PyhfData(signal, [ws] )
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.getUpperLimitOnMu()
        self.assertEqual(ulcomputer.workspaces, None)
        self.assertEqual(ul, None)
        # Missing observations
        ws = dict(channels=channels,
                  measurements=measurements,
                  #observations=observations,
                  version='1.0.0'
                  )
        data = PyhfData(signal, [ws] )
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.getUpperLimitOnMu()
        self.assertEqual(ulcomputer.workspaces, None)
        self.assertEqual(ul, None)
        # Missing version
        ws = dict(channels=channels,
                  measurements=measurements,
                  observations=observations,
                  #version='1.0.0'
                  )
        data = PyhfData(signal, [ws] )
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.getUpperLimitOnMu()
        self.assertIsNone(ulcomputer.workspaces)
        self.assertIsNone(ul)

    def testNoSignal(self):
        """
        Tests the case where all SRs are empty
        """
        ws = self.simpleJson([0.9], [10])
        signal = { "dummy0": { "SR1": 0. } }
        data = PyhfData(signal, [ws])
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.getUpperLimitOnMu()
        self.assertIsNone(ul)

    def testWrongNbOfSignals(self):
        """
        Tests the case where the number of SRs feeded to the module doesn't match the number of SRs in the json
        """
        # One single json but too much signals
        ws = self.simpleJson([0.9], [10])
        nsignals = [[0.9, 0.5]]
        nsignals = { "dummy0": { "SR1": 0.9, "SR2": 0.5 } }
        data = PyhfData( nsignals, [ws] )
        ulcomputer = PyhfUpperLimitComputer(data)
        ul1 = ulcomputer.getUpperLimitOnMu()
        # Two jsons but only one signal
        ws = [self.simpleJson([0.9], [10]), self.simpleJson([0.8], [9])]
        signal = { "dummy0": { "SR1": 0.5 } }
        data = PyhfData( signal, ws )
        ulcomputer = PyhfUpperLimitComputer(data)
        ul2 = ulcomputer.getUpperLimitOnMu(workspace_index=0)
        self.assertIsNone(ul1)
        self.assertIsNone(ul2)

    def testWSindex(self):
        """
        Tests how the module reacts when giving several jsons but not specifying for
        which the UL should be computed
        """
        bg = [ .9, .8 ]
        obs = [ 10, 9 ]
        # nsig = [ .1, .2 ]
        ws = [ self.simpleJson([x], [y]) for x,y in zip (bg,obs) ]
        # nsignals = [ [x] for x in nsig ]
        nsignals = { "dummy0": { "SR1" : .1 }, "dummy1": { "SR1": .2 } }
        data = PyhfData( nsignals, ws)
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.getUpperLimitOnMu()
        self.assertAlmostEqual ( ul, 70.44942696708914, 1 )
        """ compare with:
        m = Data( observed=[10,9], backgrounds=[.9,.8],
        covariance=[[1e-6,0],[0,1e-6]], nsignal=[.1,.2])
        ulComp = UpperLimitComputer()
        ul = ulComp.getUpperLimitOnMu(m )
        """

    def testFullPyhfModule1(self):
        """
        Computes the UL using the pyhfInterface module and checks if, outside of the module, this UL still gives a 95% CLs
        """
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        bkg = self.simpleJson([0.8], [10])
        signals = [ 0.4 ]
        # Make the patch by hand
        patch = [dict(
            op='add',
            path='/channels/0/samples/0',
            value=dict(
                name='SR1',
                data=signals,
                modifiers=[
                    dict(
                        name='lumi',
                        type='lumi',
                        data=None
                    ),
                    dict(
                        name='mu_SIG',
                        type='normfactor',
                        data=None
                    )
                ]
            )
        )]
        signals = { "dummy0": { "SR1": 0.4 } }
        llhdSpec = jsonpatch.apply_patch(bkg, patch)
        # Computing the upper limit with the SModelS/pyhf interface
        data = PyhfData(signals, [bkg])
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.getUpperLimitOnMu()
        # ul = ul * data.totalYield
        # Computing the cls outside of SModelS with POI = ul, should give 0.95
        msettings = {'normsys': {'interpcode': 'code4'}, 'histosys': {'interpcode': 'code4p'}}
        workspace = pyhf.Workspace(llhdSpec)
        model = workspace.model(modifier_settings=msettings)
        bounds = model.config.suggested_bounds()
        bounds[model.config.poi_index] = [0,100]
        args = { "return_expected": False }
        pver = float ( pyhf.__version__[:3] )
        if pver < 0.6:
            args["qtilde"]=True
        else:
            args["test_stat"]="qtilde"
        result = pyhf.infer.hypotest(ul, workspace.data(model), model, par_bounds=bounds, **args )
        try:
            CLs = float(result[0])
        except IndexError:
            CLs = float(result)
        self.assertAlmostEqual(CLs, 0.05, 2)

    def testFullPyhfModule2(self):
        """
        Same as previous but with two SRs
        """
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        bkg = self.simpleJson([0.8, 0.9], [10, 11])
        signals = [0.4, 0.2]
        # Make the patch by hand
        patch = [dict(
            op='add',
            path='/channels/0/samples/0',
            value=dict(
                name='SR1',
                data=signals,
                modifiers=[
                    dict(
                        name='lumi',
                        type='lumi',
                        data=None
                    ),
                    dict(
                        name='mu_SIG',
                        type='normfactor',
                        data=None
                    )
                ]
            )
        )]
        llhdSpec = jsonpatch.apply_patch(bkg, patch)
        # Computing the upper limit with the SModelS/pyhf interface
        signals = { "dummy0": { "SR1[0]": 0.4, "SR1[1]": 0.2 } }
        data = PyhfData(signals, [bkg] )
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.getUpperLimitOnMu()
        # Computing the cls outside of SModelS with POI = ul, should give 0.95
        msettings = {'normsys': {'interpcode': 'code4'}, 'histosys': {'interpcode': 'code4p'}}
        workspace = pyhf.Workspace(llhdSpec)
        model = workspace.model(modifier_settings=msettings)
        bounds = model.config.suggested_bounds()
        bounds[model.config.poi_index] = [0,100]
        args = { "return_expected": False }
        pver = float ( pyhf.__version__[:3] )
        if pver < 0.6:
            args["qtilde"]=True
        else:
            args["test_stat"]="qtilde"
        result = pyhf.infer.hypotest(ul, workspace.data(model), model, par_bounds=bounds, **args )
        try:
            CLs = float(result[0])
        except IndexError:
            CLs = float(result)
        self.assertAlmostEqual(CLs, 0.05, 2)

    def testPatchUncertaintyAndCRs(self):
        """
        Test if signal uncertainty and signal in the CRs are pathed, with ATLAS-SUSY-2018-16 as a test case
        """
        from smodels.base.physicsUnits import GeV
        from smodels.experiment.databaseObj import Database
        from smodels.base import runtime
        from smodels.tools.particlesLoader import load
        from smodels.base.model import Model
        from smodels.share.models.SMparticles import SMList
        import os
        from smodels.decomposition import decomposer
        from smodels.matching.theoryPrediction import _getDataSetPredictions, _getCombinedResultFor, TheoryPredictionList
        warnings.filterwarnings("ignore", category=DeprecationWarning)

        database = Database('./database/')
        runtime.modelFile = "smodels.share.models.mssm"
        BSMList = load()

        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        slhafile = os.path.abspath('./testFiles/slha/TChiWZoff_150_125_150_125.slha')
        model.updateParticles(inputFile=slhafile,ignorePromptQNumbers = ['eCharge','colordim','spin'])

        # Set main options for decomposition
        topDict = decomposer.decompose(model, sigmacut=0.001,
                               massCompress=True, invisibleCompress=True,
                               minmassgap=5*GeV)

        expResult = database.getExpResults(analysisIDs=['ATLAS-SUSY-2018-16'],
                                            datasetIDs=['all'],
                                            txnames=['TChiWZoff'],
                                            dataTypes=['efficiencyMap'])[0]

        expSMSDict = database.expSMSDict
        smsMatch = expSMSDict.getMatchesFrom(topDict)
        dataSetResults = []
        for dataset in expResult.datasets:
            predList = _getDataSetPredictions(dataset, smsMatch, expSMSDict, maxMassDist=0.2)
            if predList:
                dataSetResults.append(predList)

        combinedRes = _getCombinedResultFor(dataSetResults,expResult)
        combinedRes.setStatsComputer()

        patch = combinedRes._statsComputer.likelihoodComputer.patches[0]

        with open('patch_2018-16.json','r') as f:
            patch_ref = json.load(f)

        equals = patch == patch_ref
        if not equals:
            logger.error ( f"patch_2018-16.json != debug.json" )
            with open('debug.json','w') as f:
                json.dump(patch,f,indent=4)
            
        self.assertEqual(patch,patch_ref)

if __name__ == "__main__":
    unittest.main()
