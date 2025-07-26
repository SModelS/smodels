#!/usr/bin/env python3

"""
.. module:: testMassCompress
   :synopsis: Tests mass compression.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
sys.path.insert(0,"../")
import unittest
from smodels.base.smodelsLogging import setLogLevel
from smodels.tools.particlesLoader import load
from smodels.base.model import Model
from smodels.share.models.SMparticles import SMList
from smodels.decomposition.theorySMS import TheorySMS
from smodels.experiment.expSMS import ExpSMS
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model
from smodels.base.physicsUnits import fb, GeV,MeV

setLogLevel('error')


class MassCompressTest(unittest.TestCase):
    definingRun = False ## meant only to adapt to changes in output format
    ## use with super great care!!

    def testInvCompA(self):
        
        slhafile = './testFiles/slha/lightEWinos.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile, promptWidth=1e-10*GeV,
                              ignorePromptQNumbers=['spin','eCharge','colordim'])
        
        stringEl = "(PV > H0(1),H-(2)), (H0(1) > H+(3),H-(4)), (H+(3) > nu,higgs), (H-(4) > nu,nu), (H-(2) > N1,N2)"
        # Hack to create a theory element from a string:
        expSMS = ExpSMS.from_string(stringEl, model=model)
        treeA = TheorySMS()
        treeA.add_nodes_from(expSMS.nodes)
        treeA.add_edges_from(expSMS.edgeIndices)
        treeA.prodXSec = 1.0*fb
        treeA.maxWeight = 1.0*fb
        treeA.setGlobalProperties()

        self.assertEqual([str(n) for n in treeA.nodes],
                         ['PV', 'H-', 'H0', 'N1', 'N2', 'H+', 'H-', 'higgs', 'nu', 'nu', 'nu'])
        self.assertEqual(treeA.nodeIndices,[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

        smsComp = treeA.invisibleCompress()
        self.assertEqual([str(n) for n in smsComp.nodes],
                         ['PV', 'inv', 'H0', 'inv', 'H+', 'higgs', 'nu'])
        self.assertEqual(smsComp.nodeIndices,[0, 1, 2, 3, 4, 5, 6])

    def testInvCompB(self):
        
        slhafile = './testFiles/slha/lightEWinos.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile, promptWidth=1e-10*GeV,
                              ignorePromptQNumbers=['spin','eCharge','colordim'])
        
        stringEl = "(PV > H0(1),H-(2)), (H0(1) > H+(3),H-(4)), (H+(3) > nu,nu), (H-(4) > nu,nu), (H-(2) > N1,N2)"
        # Hack to create a theory element from a string:
        expSMS = ExpSMS.from_string(stringEl, model=model)
        treeA = TheorySMS()
        treeA.add_nodes_from(expSMS.nodes)
        treeA.add_edges_from(expSMS.edgeIndices)
        treeA.prodXSec = 1.0*fb
        treeA.maxWeight = 1.0*fb
        treeA.setGlobalProperties()

        self.assertEqual([str(n) for n in treeA.nodes],
                         ['PV', 'H-', 'H0', 'N1', 'N2', 'H+', 'H-', 'nu', 'nu', 'nu', 'nu'])
        self.assertEqual(treeA.nodeIndices,[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

        smsComp = treeA.invisibleCompress()
        self.assertEqual([str(n) for n in smsComp.nodes],
                         ['PV', 'inv', 'inv'])
        self.assertEqual(smsComp.nodeIndices,[0, 1, 2])


    def testInvCompC(self):
        
        slhafile = './testFiles/slha/lightEWinos.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile, promptWidth=1e-10*GeV,
                              ignorePromptQNumbers=['spin','eCharge','colordim'])
        
        stringEl = "(PV > H0(1),H-(2)), (H0(1) > H+(3),H-(4)), (H+(3) > nu,nu,e-), (H-(4) > nu,nu), (H-(2) > N1,N2)"
        # Hack to create a theory element from a string:
        expSMS = ExpSMS.from_string(stringEl, model=model)
        treeA = TheorySMS()
        treeA.add_nodes_from(expSMS.nodes)
        treeA.add_edges_from(expSMS.edgeIndices)
        treeA.prodXSec = 1.0*fb
        treeA.maxWeight = 1.0*fb
        treeA.setGlobalProperties()

        self.assertEqual([str(n) for n in treeA.nodes],
                         ['PV', 'H-', 'H0', 'N1', 'N2', 'H-', 'H+', 'nu', 'nu', 'e-', 'nu', 'nu'])
        self.assertEqual(treeA.nodeIndices,[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])

        smsComp = treeA.invisibleCompress()
        self.assertEqual([str(n) for n in smsComp.nodes],
                         ['PV', 'inv', 'H0', 'inv', 'H+', 'e-', 'nu', 'nu'])
        self.assertEqual(smsComp.nodeIndices,[0, 1, 2, 3, 4, 5, 6, 7])

    def testInvCompD(self):
        
        slhafile = './testFiles/slha/lightEWinos.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile, promptWidth=1e-10*GeV,
                              ignorePromptQNumbers=['spin','eCharge','colordim'])
        
        stringEl = "(PV > C2+(1),N3(2)), (N3(2) > C1-(3),W+), (C1-(3) > N1(4),e-,nu), (C2+(1) > C1+(5),Z), (C1+(5) > N2(6),W+), (N2(6) > N1(7),nu,nu), (N1(7) > sta_1(8),ta+)"
        # Hack to create a theory element from a string:
        expSMS = ExpSMS.from_string(stringEl, model=model)
        treeA = TheorySMS()
        treeA.add_nodes_from(expSMS.nodes)
        treeA.add_edges_from(expSMS.edgeIndices)
        treeA.prodXSec = 1.0*fb
        treeA.maxWeight = 1.0*fb
        treeA.setGlobalProperties()

        self.assertEqual([str(n) for n in treeA.nodes],
                         ['PV', 'N3', 'C2+', 'W+', 'C1-', 'Z', 'C1+', 'N1', 'e-', 'nu', 'W+', 'N2', 'nu', 'nu', 'N1', 'sta_1', 'ta+'])
        self.assertEqual(treeA.nodeIndices,[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])

        smsComp = treeA.invisibleCompress()
        self.assertEqual(smsComp,None)

if __name__ == "__main__":
    unittest.main()
