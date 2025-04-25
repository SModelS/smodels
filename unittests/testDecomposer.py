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
from smodels.experiment.expSMS import ExpSMS
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model
from smodels.decomposition import decomposer
from smodels.base.physicsUnits import GeV,pb,TeV,fb
from smodels.experiment.defaultFinalStates import finalStates
from unitTestHelpers import theorySMSFromString as fromString
from unitTestHelpers import canonNameToVertNumb


class DecomposerTest(unittest.TestCase):

    def testDecomposerLHE(self):

        filename = "./testFiles/lhe/simplyGluino.lhe"
        model = Model(BSMList,SMList)
        model.updateParticles(filename,ignorePromptQNumbers=['spin','eCharge','colordim'])

        topList = decomposer.decompose(model, sigmacut=0*fb)
        self.assertTrue(len(topList.getSMSList()) == 1)
        sms = topList.getSMSList()[0]
        smsB = ExpSMS.from_string("[[['q','q']],[['q','q']]]",finalState=['MET','MET'], model=finalStates)
        self.assertTrue(smsB == sms)
        bsmLabels = [str(node) for node in sms.nodes
                     if node.isSM is False]
        self.assertEqual(bsmLabels,['gluino','gluino','N1', 'N1'])
        mass = [node.mass for node in sms.nodes if node.isSM is False]
        self.assertAlmostEqual(mass,[675.*GeV,675.*GeV,200.*GeV,200.*GeV])
        self.assertAlmostEqual(sms.maxWeight/1000.,0.262,3)


    def testDecomposerSLHA(self):

        filename = "./testFiles/slha/simplyGluino.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(filename,ignorePromptQNumbers=['spin','eCharge','colordim'])

        topList = decomposer.decompose(model, sigmacut=0*fb)
        self.assertTrue(len(topList.getSMSList()) == 1)
        sms = topList.getSMSList()[0]
        smsB = ExpSMS.from_string("[[['q','q']],[['q','q']]]",finalState=['MET','MET'], model=finalStates)
        self.assertTrue(smsB == sms)
        bsmLabels = [str(node) for node in sms.nodes
                     if node.isSM is False]
        self.assertEqual(bsmLabels,['gluino','gluino','N1', 'N1'])
        mass = [node.mass for node in sms.nodes if node.isSM is False]
        self.assertAlmostEqual(mass,[675.*GeV,675.*GeV,200.*GeV,200.*GeV])
        xsec = [xsec for xsec in sms.weightList if xsec.info.sqrts == 8.*TeV][0]
        xsec = xsec.value.asNumber(pb)
        self.assertAlmostEqual(xsec,0.572,3)


    def testDecomposerLongLived(self):

        filename = "./testFiles/slha/longLived.slha"
        #Consider a simpler model
        newModel = [ptc for ptc in BSMList if not isinstance(ptc.pdg,list) and abs(ptc.pdg) in [1000015,1000022]]
        model = Model(newModel,SMList)
        model.updateParticles(filename,ignorePromptQNumbers=['spin','eCharge','colordim'])

        topList = decomposer.decompose(model, sigmacut=0*fb)
        self.assertTrue(len(topList.getSMSList()) == 10)
        expectedWeights = {str(sorted(['N1','N1'])).replace(' ','') : 0.020,
                           str(sorted(['sta_1','sta_1~'])).replace(' ','') : 0.26,
                           str(sorted(['sta_1','sta_1~','N1~'])).replace(' ','') : 0.13,
                           str(sorted(['sta_1~','sta_1','N1'])).replace(' ','') : 0.13,
                           str(sorted(['sta_1~','sta_1','N1~','N1'])).replace(' ','') : 0.065}

        for sms in topList.getSMSList():
            bsmLabels = str(sorted([str(node) for node in sms.nodes if node.isSM is False]))
            bsmLabels = bsmLabels.replace(' ','')
            xsec = sms.weightList.getXsecsFor(8.*TeV)[0].value.asNumber(fb)
            self.assertAlmostEqual(expectedWeights[bsmLabels], xsec,2)

    def testCompression(self):

        filename = "./testFiles/slha/higgsinoStop.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(filename,promptWidth=1e-12*GeV,
                            ignorePromptQNumbers=['spin','eCharge','colordim']) #Force charginos/neutralinos to be considered as prompt


        tested = False
        topos = decomposer.decompose(model, sigmacut=0.1*fb, massCompress=False, invisibleCompress=False, 
                                     minmassgap=5.*GeV, minmassgapISR=5.*GeV )
        toposExpected = {'[][1]' : 1,'[][2]' : 14,'[1][1]' : 1,'[1][2]' : 25,'[2][2]' : 72,
                         '[][2,2]' : 44,'[1][1,1]' : 2,'[1][1,2]' : 22,'[2][1,2]' : 48,
                         '[2][2,2]' : 284,'[1,1][1,1]' : 5,'[1,1][1,2]' : 22,'[1,2][1,2]' : 120,
                         '[1][1,1,1]' : 2,'[1][1,2,2]' : 64,'[1,1][1,1,1]' : 12,'[1,1][1,1,2]' : 16,
                         '[1,2][1,2,2]' : 240,'[1,1,1][1,1,2]' : 56,'[1,1,2][1,1,2]' : 16,
                         '[1][1,1,1,2]' : 4,'[1,1][1,1,1,2]' : 56,'[1,1,2][1,1,1,2]' : 176}


        for topo in topos:
            vertnumb = canonNameToVertNumb(topos,topo)
            self.assertEqual(len(topos[topo]),toposExpected[vertnumb])
            if vertnumb != "[1,1][1,1]":
                continue
            sms = topos[topo][4]
            evenParticles = sms.treeToBrackets()[0]
            evenParticles = str(evenParticles).replace("'","").replace(' ', '')
            self.assertEqual(evenParticles,"[[[t+],[t-]],[[q],[W+]]]")
            self.assertEqual(len(sms.ancestors),1)

        topos = decomposer.decompose(model, sigmacut=0.1*fb, massCompress=False, invisibleCompress=True, 
                                     minmassgap=5.*GeV, minmassgapISR=5*GeV )
        toposExpected = {"[][]" : 1,"[][1]" : 2,"[][2]" : 22,"[1][1]" : 4,"[1][2]" : 25,"[2][2]" : 72,
                         "[][2,2]" : 44,"[1][1,1]" : 4,"[1][1,2]" : 42,"[2][1,2]" : 48,"[2][2,2]" : 284,
                         "[1,1][1,1]" : 5,"[1,1][1,2]" : 22,"[1,2][1,2]" : 120,"[1][1,1,1]" : 2,"[1][1,2,2]" : 72,
                         "[1,1][1,1,1]" : 20,"[1,1][1,1,2]" : 16,"[1,2][1,2,2]" : 240,"[1,1,1][1,1,2]" : 56,
                         "[1,1,2][1,1,2]" : 16,"[1][1,1,1,2]" : 4,"[1,1][1,1,1,2]" : 56,"[1,1,2][1,1,1,2]" : 176}
        for topo in topos:
            # Convert to old notation
            vertnumb = canonNameToVertNumb(topos,topo)
            self.assertEqual(len(topos[topo]),toposExpected[vertnumb])
            if vertnumb != "[][]":
                continue
            sms = topos[topo][0]
            evenParticles = sms.treeToBrackets()[0]
            evenParticles = str(evenParticles).replace("'","").replace(' ', '')
            smsMom = sms.ancestors[0]
            evenParticles = smsMom.treeToBrackets()[0]
            evenParticles = str(evenParticles).replace("'","").replace(' ', '')
            self.assertEqual(evenParticles,"[[],[[nu,nu]]]")
            bsmLabels = sorted([str(node) for node in sms.nodes if node.isSM is False])
            self.assertEqual(bsmLabels,['N1','inv'])
            ## all neutrinos are considered as equal, so there should be a single mother:
            self.assertEqual(len(sms.ancestors), 1)


        topos = decomposer.decompose(model, sigmacut=0.1*fb, massCompress=True, invisibleCompress=False, 
                                     minmassgap=5.*GeV, minmassgapISR=5*GeV)
        toposExpected = {"[][]" : 1,"[][1]" : 8,"[][2]" : 14,"[1][1]" : 4,"[1][2]" : 29,"[2][2]" : 72,
                        "[][1,2]" : 2,"[][2,2]" : 44,"[1][1,1]" : 4,"[1][1,2]" : 34,"[2][1,2]" : 48,
                        "[2][2,2]" : 284,"[1,1][1,1]" : 17,"[1,1][1,2]" : 22,"[1,2][1,2]" : 120,
                        "[1][1,1,1]" : 4,"[1][1,2,2]" : 64,"[1,1][1,1,1]" : 48,"[1,1][1,1,2]" : 16,
                        "[1,2][1,2,2]" : 240,"[1,1,1][1,1,2]" : 64,"[1,1,2][1,1,2]" : 16,"[1][1,1,1,2]" : 4,
                        "[1,1][1,1,1,2]" : 64,"[1,1,2][1,1,1,2]" : 176}
        for topo in topos:
            vertnumb = canonNameToVertNumb(topos,topo)
            self.assertEqual(len(topos[topo]),toposExpected[vertnumb])
            if vertnumb !="[1][1]":
                continue
            sms = topos[topo][2]
            evenParticles = sms.treeToBrackets()[0]
            evenParticles = str(evenParticles).replace("'","").replace(' ', '')
            self.assertEqual(evenParticles,"[[[b]],[[b]]]")
            masses = [node.mass for node in sms.ancestors[0].nodes if node.isSM is False]
            dm = abs(masses[2]-masses[4])/GeV
            self.assertEqual(len(sms.ancestors),24)
            self.assertTrue(dm < 5.0)


        topos = decomposer.decompose(model, sigmacut=0.1*fb, massCompress=True, invisibleCompress=True, 
                                     minmassgap=5.*GeV, minmassgapISR=5*GeV )
        smsIDs = {29+8 : ExpSMS.from_string("[[['b']],[['b']]]",finalState=['MET','MET'], model=finalStates),
                 30+8 : ExpSMS.from_string("[[['b']],[['t+']]]",finalState=['MET','MET'], model=finalStates),
                 32+8 : ExpSMS.from_string("[[['b']],[['t+']]]",finalState=['MET','MET'], model=finalStates),
                 26+8 : ExpSMS.from_string("[[['t-']],[['b']]]",finalState=['MET','MET'], model=finalStates),
                 31+8 : ExpSMS.from_string("[[['b']],[['t-']]]",finalState=['MET','MET'], model=finalStates),
                 33+8 : ExpSMS.from_string("[[['t+']],[['t-']]]",finalState=['MET','MET'], model=finalStates),
                 28+8 : ExpSMS.from_string("[[['t-']],[['t+']]]",finalState=['MET','MET'], model=finalStates),
                 34+8 : ExpSMS.from_string("[[['t-']],[['t+']]]",finalState=['MET','MET'], model=finalStates),
                 27+8 : ExpSMS.from_string("[[['t-']],[['t+']]]",finalState=['MET','MET'], model=finalStates)}

        toposExpected = {"[][]" : 2,"[][1]" : 9,"[][2]" : 22,"[1][1]" : 9,"[1][2]" : 29,
                         "[2][2]" : 72,"[][1,2]" : 2,"[][2,2]" : 44,"[1][1,1]" : 6,"[1][1,2]" : 54,
                         "[2][1,2]" : 48,"[2][2,2]" : 284,"[1,1][1,1]" : 17,"[1,1][1,2]" : 22,
                         "[1,2][1,2]" : 120,"[1][1,1,1]" : 4,"[1][1,2,2]" : 72,"[1,1][1,1,1]" : 56,
                         "[1,1][1,1,2]" : 16,"[1,2][1,2,2]" : 240,"[1,1,1][1,1,2]" : 64,
                         "[1,1,2][1,1,2]" : 16,"[1][1,1,1,2]" : 4,"[1,1][1,1,1,2]" : 64,
                         "[1,1,2][1,1,1,2]" : 176}
        for topo in topos:
            vertnumb = canonNameToVertNumb(topos,topo)
            self.assertEqual(len(topos[topo]),toposExpected[vertnumb])
            if vertnumb !="[1][1]":
                continue
            for sms in topos[topo]:
                self.assertTrue(sms.smsID in smsIDs)
                self.assertEqual(smsIDs[sms.smsID],sms)


    def test_extra1(self):

        slhafile = './testFiles/slha/lightEWinos.slha'
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        model.updateParticles(inputFile=slhafile,
                                    ignorePromptQNumbers=['spin','eCharge','colordim'])
        
        sigmacut = 10*fb
        topDict = decomposer.decompose(model, sigmacut= sigmacut, 
                                       massCompress=False, 
                                       invisibleCompress=False)
        nUnique = 103 # uncompressed
        self.assertEqual(nUnique,len(topDict.getSMSList()))

        topSummary = []
        for c in sorted(topDict.keys()):
            total = 0.0*fb
            for sms in topDict[c]:
                total += sms.weightList.getMaxXsec()
            topSummary.append([len(topDict[c]),round(total.asNumber(fb),3)])

        expectedSummary = [[2, 107.528], [58, 6126.453], [1, 31.861], [7, 490.94], [4, 81.12], [27, 1237.813], [4, 85.382]]
        self.assertEqual(sorted(topSummary),sorted(expectedSummary))

        # Smaller sigmacut:
        sigmacut = 0.1*fb
        topDict = decomposer.decompose(model, sigmacut= sigmacut, 
                                       massCompress=False, 
                                       invisibleCompress=False)
        nUnique = 5435
        self.assertEqual(nUnique,len(topDict.getSMSList()))

        topSummary = []
        for c in sorted(topDict.keys()):
            total = 0.0*fb
            for sms in topDict[c]:
                total += sms.weightList.getMaxXsec()
        #     print(c,len(topDict[c]),total)
            topSummary.append([len(topDict[c]),round(total.asNumber(fb),3)])
        expectedSummary = [[1, 3.883], [2, 1.003], [3, 1.094], [4, 0.892], [4, 2.616], [6, 2.853], [8, 3.643], [10, 3.376], [10, 14.988], [11, 10.526], [13, 13.994], [14, 11.171], [15, 164.364], [16, 78.045], [17, 30.516], [25, 12.91], [28, 19.106], [49, 42.873], [56, 25.215], [62, 48.534], [78, 51.703], [87, 150.926], [90, 6350.441], [102, 154.363], [104, 94.307], [104, 153.527], [107, 404.754], [113, 403.824], [126, 47.333], [132, 1078.776], [144, 267.975], [183, 106.275], [295, 166.074], [304, 377.261], [411, 157.421], [504, 2817.347], [648, 455.245], [747, 2086.069], [802, 1330.765]]
        self.assertEqual(sorted(topSummary),sorted(expectedSummary))

        # Invisible compression
        sigmacut = 1*fb
        topDict = decomposer.decompose(model, sigmacut= sigmacut, 
                                       massCompress=False, 
                                       invisibleCompress=True)

        nUnique = 917
        self.assertEqual(nUnique,len(topDict.getSMSList()))

        # Mass compression
        topDict = decomposer.decompose(model, sigmacut= sigmacut, 
                                       massCompress=True, 
                                       invisibleCompress=False)
        nUnique = 899
        self.assertEqual(nUnique,len(topDict.getSMSList()))


    def test_extra2(self):

        slhafile="./testFiles/slha/higgsinoStop.slha"
        model = Model(BSMList,SMList)
        model.updateParticles(inputFile=slhafile,promptWidth = 1e-12*GeV,
                              
                              ignorePromptQNumbers=['spin','eCharge','colordim'])
        
        sigmacut = 1*fb
        topDict = decomposer.decompose(model, sigmacut= sigmacut, 
                                       massCompress=True, invisibleCompress=False, 
                                       minmassgap=5*GeV, minmassgapISR=5*GeV)
        nUnique = 489
        self.assertEqual(nUnique,len(topDict.getSMSList()))


        topSummary = []
        for c in sorted(topDict.keys()):
            total = 0.0*fb
            for sms in topDict[c]:
                total += sms.weightList.getMaxXsec()
            topSummary.append([len(topDict[c]),round(total.asNumber(fb),3)])

        expectedSummary = [[1, 46932.601], [1, 1043.828], [14, 68354.626], [4, 8885.261], [8, 694.931], [72, 42401.477], [36, 438.84], [4, 119.091], [34, 8107.982], [116, 792.104], [16, 111.128], [111, 8119.626], [8, 23.764], [16, 47.891], [16, 41.485], [8, 11.09], [8, 11.09], [16, 19.359]]
        self.assertEqual(sorted(topSummary),sorted(expectedSummary))

        # Mass and invisible compression
        sigmacut = 0.1*fb
        topDict = decomposer.decompose(model, sigmacut= sigmacut, 
                            massCompress=True, 
                            invisibleCompress=True, minmassgap=5*GeV, minmassgapISR=5*GeV)
        nUnique = 1452
        self.assertEqual(nUnique,len(topDict.getSMSList()))
        topSummary = []
        for c in sorted(topDict.keys()):
            total = 0.0*fb
            for sms in topDict[c]:
                total += sms.weightList.getMaxXsec()
        #     print(c,len(topDict[c]),total)
            topSummary.append([len(topDict[c]),round(total.asNumber(fb),3)])


        expectedSummary = [[2, 56785.133], [9, 1055.527], [22, 74853.168], [9, 10118.398], [29, 702.603], [72, 42401.477], [2, 0.923], [44, 443.058], [6, 126.609], [54, 9124.952], [48, 10.39], [284, 879.227], [17, 5.406], [22, 114.309], [120, 8125.92], [4, 1.076], [72, 48.575], [56, 99.821], [16, 2.469], [240, 119.334], [64, 31.64], [16, 1.667], [4, 0.684], [64, 31.64], [176, 59.448]]
        self.assertEqual(sorted(topSummary),sorted(expectedSummary))



if __name__ == "__main__":
    unittest.main()
