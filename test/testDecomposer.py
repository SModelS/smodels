#!/usr/bin/env python3

"""
.. module:: testDecomposer
   :synopsis: Tests the ascii grapher.
              Depends also on lheReader, lheDecomposer.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import unittest
import sys
sys.path.insert(0,"../")
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model
from smodels.theory import decomposer
from smodels.theory.element import Element
from smodels.tools.physicsUnits import GeV,pb,TeV,fb
from smodels.experiment.defaultFinalStates import finalStates

class DecomposerTest(unittest.TestCase):

    def testDecomposerLHE(self):

        filename = "./testFiles/lhe/simplyGluino.lhe"
        model = Model(BSMList,SMList)
        model.updateParticles(filename)

        topList = decomposer.decompose(model, sigmacut=0*fb)
        self.assertTrue(len(topList.getElements()) == 1)
        element = topList.getElements()[0]
        el = Element("[[[q,q]],[[q,q]]]",finalState=['MET','MET'], model=finalStates)
        self.assertTrue(el == element)
        bsmLabels = [[bsm.label for bsm in branch] for branch in element.oddParticles]
        self.assertEqual(bsmLabels,[['gluino','N1']]*2)
        self.assertAlmostEqual(element.mass,[[675.*GeV,200.*GeV]]*2)
        xsec = [xsec for xsec in element.weight if xsec.info.sqrts == 8.*TeV][0]
        xsec = xsec.value.asNumber(pb)
        self.assertAlmostEqual(xsec,0.262,3)


    def testDecomposerSLHA(self):

        filename = "./testFiles/slha/simplyGluino.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(filename)

        topList = decomposer.decompose(model, sigmacut=0*fb)
        self.assertTrue(len(topList.getElements()) == 1)
        element = topList.getElements()[0]
        el = Element("[[[q,q]],[[q,q]]]",finalState=['MET','MET'], model=finalStates)
        self.assertTrue(el == element)
        bsmLabels = [[bsm.label for bsm in branch] for branch in element.oddParticles]
        self.assertEqual(bsmLabels,[['gluino','N1']]*2)
        self.assertAlmostEqual(element.mass,[[675.*GeV,200.*GeV]]*2)
        xsec = [xsec for xsec in element.weight if xsec.info.sqrts == 8.*TeV][0]
        xsec = xsec.value.asNumber(pb)
        self.assertAlmostEqual(element.weight[0].value.asNumber(pb),0.572,3)


    def testDecomposerLongLived(self):

        filename = "./testFiles/slha/longLived.slha"
        #Consider a simpler model
        newModel = [ptc for ptc in BSMList if not isinstance(ptc.pdg,list) and abs(ptc.pdg) in [1000015,1000022]]
        model = Model(newModel,SMList)
        model.updateParticles(filename)

        topList = decomposer.decompose(model, sigmacut=0*fb)
        self.assertTrue(len(topList.getElements()) == 10)
        expectedWeights = {str(sorted([['N1'],['N1']])).replace(' ','') : 0.020,
                           str(sorted([['sta_1'],['sta_1~']])).replace(' ','') : 0.26,
                           str(sorted([['sta_1'],['sta_1~','N1~']])).replace(' ','') : 0.13,
                           str(sorted([['sta_1~'],['sta_1','N1']])).replace(' ','') : 0.13,
                           str(sorted([['sta_1~','N1~'],['sta_1','N1']])).replace(' ','') : 0.065}

        for el in topList.getElements():
            bsmLabels = str(sorted([[bsm.label for bsm in branch] for branch in el.oddParticles]))
            bsmLabels = bsmLabels.replace(' ','')
            xsec = el.weight.getXsecsFor(8.*TeV)[0].value.asNumber(fb)
            self.assertAlmostEqual(expectedWeights[bsmLabels], xsec,2)

    def testCompression(self):

        filename = "./testFiles/slha/higgsinoStop.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(filename,promptWidth=1e-12*GeV) #Force charginos/neutralinos to be considered as prompt


        tested = False
        topos = decomposer.decompose(model, sigmacut=0.1*fb, doCompress=False, doInvisible=False, minmassgap=5.*GeV )
        toposExpected = {'[][1]' : 1,'[][2]' : 14,'[1][1]' : 1,'[1][2]' : 25,'[2][2]' : 72,
                         '[][2,2]' : 44,'[1][1,1]' : 2,'[1][1,2]' : 22,'[2][1,2]' : 48,
                         '[2][2,2]' : 284,'[1,1][1,1]' : 5,'[1,1][1,2]' : 22,'[1,2][1,2]' : 120,
                         '[1][1,1,1]' : 2,'[1][1,2,2]' : 64,'[1,1][1,1,1]' : 12,'[1,1][1,1,2]' : 16,
                         '[1,2][1,2,2]' : 240,'[1,1,1][1,1,2]' : 56,'[1,1,2][1,1,2]' : 16,
                         '[1][1,1,1,2]' : 4,'[1,1][1,1,1,2]' : 56,'[1,1,2][1,1,1,2]' : 176}


        for topo in topos:
            self.assertEqual(len(topo.elementList),toposExpected[str(topo)])
            if str(topo)!="[1,1][1,1]":
                continue
            for element in topo.elementList:
                if str(element)!="[[[t+],[t-]],[[q],[W+]]]":
                    continue
                tested = True
                self.assertEqual(len(element.motherElements),1)
        self.assertTrue(tested) #Make sure the test was performed

        tested = False
        topos = decomposer.decompose(model, sigmacut=0.1*fb, doCompress=False, doInvisible=True, minmassgap=5.*GeV )
        toposExpected = {"[][]" : 1,"[][1]" : 2,"[][2]" : 22,"[1][1]" : 4,"[1][2]" : 25,"[2][2]" : 72,
                         "[][2,2]" : 44,"[1][1,1]" : 4,"[1][1,2]" : 42,"[2][1,2]" : 48,"[2][2,2]" : 284,
                         "[1,1][1,1]" : 5,"[1,1][1,2]" : 22,"[1,2][1,2]" : 120,"[1][1,1,1]" : 2,"[1][1,2,2]" : 72,
                         "[1,1][1,1,1]" : 20,"[1,1][1,1,2]" : 16,"[1,2][1,2,2]" : 240,"[1,1,1][1,1,2]" : 56,
                         "[1,1,2][1,1,2]" : 16,"[1][1,1,1,2]" : 4,"[1,1][1,1,1,2]" : 56,"[1,1,2][1,1,1,2]" : 176}
        for topo in topos:
            self.assertEqual(len(topo.elementList),toposExpected[str(topo)])
            if str(topo)!="[][]":
                continue
            for element in topo.elementList:
                if str(element) != "[[],[]]":
                    continue
                tested = True
                self.assertEqual(str(element.motherElements[0]),"[[],[[nu,nu]]]")
                bsmLabels = [[bsm.label for bsm in br] for br in element.oddParticles]
                self.assertEqual(bsmLabels,[['N1'],['inv']])
                ## all neutrinos are considered as equal, so there should be a single mother:
                self.assertEqual(len(element.motherElements), 1)
        self.assertTrue(tested) #Make sure the test was performed


        tested = False
        topos = decomposer.decompose(model, sigmacut=0.1*fb, doCompress=True, doInvisible=False, minmassgap=5.*GeV)
        toposExpected = {"[][]" : 1,"[][1]" : 8,"[][2]" : 14,"[1][1]" : 4,"[1][2]" : 29,"[2][2]" : 72,
                        "[][1,2]" : 2,"[][2,2]" : 44,"[1][1,1]" : 4,"[1][1,2]" : 34,"[2][1,2]" : 48,
                        "[2][2,2]" : 284,"[1,1][1,1]" : 17,"[1,1][1,2]" : 22,"[1,2][1,2]" : 120,
                        "[1][1,1,1]" : 4,"[1][1,2,2]" : 64,"[1,1][1,1,1]" : 48,"[1,1][1,1,2]" : 16,
                        "[1,2][1,2,2]" : 240,"[1,1,1][1,1,2]" : 64,"[1,1,2][1,1,2]" : 16,"[1][1,1,1,2]" : 4,
                        "[1,1][1,1,1,2]" : 64,"[1,1,2][1,1,1,2]" : 176}
        for topo in topos:
            self.assertEqual(len(topo.elementList),toposExpected[str(topo)])
            if str(topo)!="[1][1]":
                continue
            for element in topo.elementList:
                if str(element)!="[[[b]],[[b]]]":
                    continue
                masses = element.motherElements[0].mass
                dm = abs(masses[0][1]-masses[0][2])/GeV
                tested = True
                self.assertEqual(len(element.motherElements),24)
                self.assertTrue(dm < 5.0)
        self.assertTrue(tested) #Make sure the test was performed

        tested = False
        topos = decomposer.decompose(model, sigmacut=0.1*fb, doCompress=True, doInvisible=True, minmassgap=5.*GeV )
        elIDs = {29+8 : Element("[[[b]],[[b]]]",finalState=['MET','MET'], model=finalStates),
                 30+8 : Element("[[[b]],[[t+]]]",finalState=['MET','MET'], model=finalStates),
                 32+8 : Element("[[[b]],[[t+]]]",finalState=['MET','MET'], model=finalStates),
                 26+8 : Element("[[[t-]],[[b]]]",finalState=['MET','MET'], model=finalStates),
                 31+8 : Element("[[[b]],[[t-]]]",finalState=['MET','MET'], model=finalStates),
                 33+8 : Element("[[[t+]],[[t-]]]",finalState=['MET','MET'], model=finalStates),
                 28+8 : Element("[[[t-]],[[t+]]]",finalState=['MET','MET'], model=finalStates),
                 34+8 : Element("[[[t-]],[[t+]]]",finalState=['MET','MET'], model=finalStates),
                 27+8 : Element("[[[t-]],[[t+]]]",finalState=['MET','MET'], model=finalStates)}

        toposExpected = {"[][]" : 2,"[][1]" : 9,"[][2]" : 22,"[1][1]" : 9,"[1][2]" : 29,
                         "[2][2]" : 72,"[][1,2]" : 2,"[][2,2]" : 44,"[1][1,1]" : 6,"[1][1,2]" : 54,
                         "[2][1,2]" : 48,"[2][2,2]" : 284,"[1,1][1,1]" : 17,"[1,1][1,2]" : 22,
                         "[1,2][1,2]" : 120,"[1][1,1,1]" : 4,"[1][1,2,2]" : 72,"[1,1][1,1,1]" : 56,
                         "[1,1][1,1,2]" : 16,"[1,2][1,2,2]" : 240,"[1,1,1][1,1,2]" : 64,
                         "[1,1,2][1,1,2]" : 16,"[1][1,1,1,2]" : 4,"[1,1][1,1,1,2]" : 64,
                         "[1,1,2][1,1,1,2]" : 176}
        for topo in topos:
            self.assertEqual(len(topo.elementList),toposExpected[str(topo)])
            if str(topo)!="[1][1]":
                continue
            for element in topo.elementList:
                if element.elID in elIDs:
                    tested = True
                    self.assertEqual(element,elIDs[element.elID])
        self.assertTrue(tested) #Make sure the test was performed

if __name__ == "__main__":
    unittest.main()
