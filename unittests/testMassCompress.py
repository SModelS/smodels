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

    def testMassCompA(self):
        
        slhafile = './testFiles/slha/higgsino_dm_4.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile, promptWidth=1e-10*GeV,
                              ignorePromptQNumbers=['spin','eCharge','colordim'])
        
        stringEl = "(PV > C2+(1),N3(2)), (N3(2) > C1-(3),W+), (C1-(3) > N1,e-,nu), (C2+(1) > N2(4),W+), (N2(4) > N1(5),nu,nu), (N1(5) > sta_1,ta+)"
        # Hack to create a theory element from a string:
        expSMS = ExpSMS.from_string(stringEl, model=model)
        treeA = TheorySMS()
        treeA.add_nodes_from(expSMS.nodes)
        treeA.add_edges_from(expSMS.edgeIndices)
        treeA.prodXSec = 1.0*fb
        treeA.maxWeight = 1.0*fb
        treeA.setGlobalProperties()

        self.assertEqual([str(n) for n in treeA.nodes],
                         ['PV', 'N3', 'C2+', 'W+', 'C1-', 'W+', 'N2', 'N1', 'e-', 'nu', 'nu', 'nu', 'N1', 'sta_1', 'ta+'])
        self.assertEqual(treeA.nodeIndices,[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])

        masses = [None, 5.676E+02*GeV, 5.68E+02*GeV, 8.00E+01*GeV, 1.232E+02*GeV, 8.00E+01*GeV, 
                  1.232E+02*GeV, 1.192E+02*GeV, 5.00E-01*MeV, 0.00E+00*MeV, 0.00E+00*MeV, 
                  0.00E+00*MeV, 1.192E+02*GeV, 5.0057E+03*GeV, 1.777E+03*MeV]
        for mA,mB in zip(treeA.mass,masses):
            if mA == mB:
                continue
            self.assertAlmostEqual(mA.asNumber(GeV),mB.asNumber(GeV),2)


        newSMS = treeA.massCompress(5*GeV,minmassgapISR=0*GeV)
        self.assertEqual([str(n) for n in newSMS.nodes],
                         ['PV', 'N3', 'C2+', 'N1', 'W+', 'W+', 'N1', 'sta_1', 'ta+'])
        masses = [None, 5.676E+02*GeV, 5.68E+02*GeV, 1.192E+02*GeV, 8.00E+01*GeV, 
                  8.00E+01*GeV, 1.192E+02*GeV, 5.0057E+03*GeV, 1.78E+03*MeV]
        for mA,mB in zip(newSMS.mass,masses):
            if mA == mB:
                continue
            self.assertAlmostEqual(mA.asNumber(GeV),mB.asNumber(GeV),2)

        newSMS = treeA.massCompress(3*GeV,0*GeV)
        self.assertEqual(newSMS,None)

    def testMassCompB(self):
        
        slhafile = './testFiles/slha/higgsino_dm_4.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile, promptWidth=1e-10*GeV,
                              ignorePromptQNumbers=['spin','eCharge','colordim'])
        
        stringEl = "(PV > C2+(1),N3(2)), (N3(2) > C1-(3),W+), (C1-(3) > N1,e-,nu), (C2+(1) > N2(4),W+), (N2(4) > N1(5),sta_1,nu), (N1(5) > sta_1,ta+)"
        # Hack to create a theory element from a string:
        expSMS = ExpSMS.from_string(stringEl, model=model)
        treeA = TheorySMS()
        treeA.add_nodes_from(expSMS.nodes)
        treeA.add_edges_from(expSMS.edgeIndices)
        treeA.prodXSec = 1.0*fb
        treeA.maxWeight = 1.0*fb
        treeA.setGlobalProperties()

        self.assertEqual([str(n) for n in treeA.nodes],
                         ['PV', 'N3', 'C2+', 'W+', 'C1-', 'W+', 'N2', 'N1', 'e-', 'nu', 'sta_1', 'nu', 'N1', 'sta_1', 'ta+'])
        self.assertEqual(treeA.nodeIndices,[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])


        newSMS = treeA.massCompress(10*GeV,0*GeV)
        self.assertEqual([str(n) for n in newSMS.nodes],
                         ['PV', 'N3', 'C2+', 'N1', 'W+', 'W+', 'N2', 'sta_1', 'nu', 'N1', 'sta_1', 'ta+'])
        masses = [None, 5.676E+02*GeV, 5.68E+02*GeV, 1.192E+02*GeV, 8.00E+01*GeV, 8.00E+01*GeV, 
                  1.232E+02*GeV, 5.0057E+03*GeV, 0.00E+00*MeV, 1.192E+02*GeV, 5.0057E+03*GeV, 1.78E+03*MeV]

        for mA,mB in zip(newSMS.mass,masses):
            if mA == mB:
                continue
            self.assertAlmostEqual(mA.asNumber(GeV),mB.asNumber(GeV),2)


    def testMassCompC(self):
        
        slhafile = './testFiles/slha/higgsino_dm_4.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile, promptWidth=1e-10*GeV,
                              ignorePromptQNumbers=['spin','eCharge','colordim'])
        
        stringEl = "(PV > C2+(1),N3(2)), (N3(2) > C1-(3),W+), (C1-(3) > N1,e-,nu), (C2+(1) > C1+(4),Z), (C1+(4) > N2(5),W+), (N2(5) > N1(6),nu,nu), (N1(6) > sta_1(7),ta+)"
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


        newSMS = treeA.massCompress(10*GeV,0*GeV)
        self.assertEqual([str(n) for n in newSMS.nodes],
                         ['PV', 'N3', 'C2+', 'N1', 'W+', 'Z', 'N1', 'sta_1', 'ta+'])
        masses = [None, 5.676E+02*GeV, 5.68E+02*GeV, 1.192E+02*GeV, 8.00E+01*GeV, 9.10E+01*GeV, 
                  1.192E+02*GeV, 5.0057E+03*GeV, 1.78E+03*MeV]


        for mA,mB in zip(newSMS.mass,masses):
            if mA == mB:
                continue
            self.assertAlmostEqual(mA.asNumber(GeV),mB.asNumber(GeV),2)



if __name__ == "__main__":
    unittest.main()
