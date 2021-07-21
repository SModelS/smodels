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
from smodels.tools.pyhfInterface import PyhfData, PyhfUpperLimitComputer, pyhf

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
        data = PyhfData([[0.1]], [ws])
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.ulSigma()
        self.assertEqual(ulcomputer.workspaces, None)
        self.assertEqual(ul, None)
        # Missing measurements
        ws = dict(channels=channels,
                  #measurements=measurements,
                  observations=observations,
                  version='1.0.0'
                  )
        data = PyhfData([[0.1]], [ws])
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.ulSigma()
        self.assertEqual(ulcomputer.workspaces, None)
        self.assertEqual(ul, None)
        # Missing observations
        ws = dict(channels=channels,
                  measurements=measurements,
                  #observations=observations,
                  version='1.0.0'
                  )
        data = PyhfData([[0.1]], [ws])
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.ulSigma()
        self.assertEqual(ulcomputer.workspaces, None)
        self.assertEqual(ul, None)
        # Missing version
        ws = dict(channels=channels,
                  measurements=measurements,
                  observations=observations,
                  #version='1.0.0'
                  )
        data = PyhfData([[0.1]], [ws])
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.ulSigma()
        self.assertIsNone(ulcomputer.workspaces)
        self.assertIsNone(ul)

    def testCorruptJson2Signal(self):
        """
        Tests how the module handles corrupted json files
        Maybe it is needed to test different types of corruptions
        """
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
        data = PyhfData([[0.1, 0.2]], [ws])
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.ulSigma()
        self.assertEqual(ulcomputer.workspaces, None)
        self.assertEqual(ul, None)
        # Missing measurements
        ws = dict(channels=channels,
                  #measurements=measurements,
                  observations=observations,
                  version='1.0.0'
                  )
        data = PyhfData([[0.1, 0.2]], [ws])
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.ulSigma()
        self.assertEqual(ulcomputer.workspaces, None)
        self.assertEqual(ul, None)
        # Missing observations
        ws = dict(channels=channels,
                  measurements=measurements,
                  #observations=observations,
                  version='1.0.0'
                  )
        data = PyhfData([[0.1, 0.2]], [ws])
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.ulSigma()
        self.assertEqual(ulcomputer.workspaces, None)
        self.assertEqual(ul, None)
        # Missing version
        ws = dict(channels=channels,
                  measurements=measurements,
                  observations=observations,
                  #version='1.0.0'
                  )
        data = PyhfData([[0.1, 0.2]], [ws])
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.ulSigma()
        self.assertIsNone(ulcomputer.workspaces)
        self.assertIsNone(ul)

    def testNoSignal(self):
        """
        Tests the case where all SRs are empty
        """
        ws = self.simpleJson([0.9], [10])
        data = PyhfData([[0]], [ws])
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.ulSigma()
        self.assertIsNone(ul)

    def testWrongNbOfSignals(self):
        """
        Tests the case where the number of SRs feeded to the module doesn't match the number of SRs in the json
        """
        # One single json but too much signals
        ws = self.simpleJson([0.9], [10])
        data = PyhfData([[0.9, 0.5]], [ws])
        ulcomputer = PyhfUpperLimitComputer(data)
        ul1 = ulcomputer.ulSigma()
        # Two jsons but only one signal
        ws = [self.simpleJson([0.9], [10]), self.simpleJson([0.8], [9])]
        data = PyhfData([[0.5]], ws)
        ulcomputer = PyhfUpperLimitComputer(data)
        ul2 = ulcomputer.ulSigma(workspace_index=0)
        self.assertIsNone(ul1)
        self.assertIsNone(ul2)

    def testWSindex(self):
        """
	      Tests how the module reacts when giving several jsons but not specifying for
        which the UL should be computed
        """
        ws = [self.simpleJson([0.9], [10]), self.simpleJson([0.8], [9])]
        data = PyhfData([[0.1], [0.2]], ws)
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.ulSigma()
        self.assertAlmostEqual ( ul, 70.44942596708917, 1 )
        # self.assertIsNone(ul)

    def testFullPyhfModule1(self):
        """
        Computes the UL using the pyhfInterface module and checks if, outside of the module, this UL still gives a 95% CLs
        """
        bkg = self.simpleJson([0.8], [10])
        signals = [0.4]
        # Make the patch by hand
        patch = [dict(
            op='add',
            path='/channels/0/samples/0',
            value=dict(
                name='sig',
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
        data = PyhfData([signals], [bkg])
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.ulSigma()
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
        bkg = self.simpleJson([0.8, 0.9], [10, 11])
        signals = [0.4, 0.2]
        # Make the patch by hand
        patch = [dict(
            op='add',
            path='/channels/0/samples/0',
            value=dict(
                name='sig',
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
        data = PyhfData([signals], [bkg])
        ulcomputer = PyhfUpperLimitComputer(data)
        ul = ulcomputer.ulSigma()
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


if __name__ == "__main__":
    unittest.main()
