#!/usr/bin/env python

"""
.. module:: testSlhaDecomposition
   :synopsis: Checks slha decomposition, alongside with compression
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

import sys
sys.path.insert(0,"../")
from smodels.theory.slhaReader import getInputData
from smodels.theory import decomposer
from smodels.tools.physicsUnits import GeV, fb
from smodels.theory.element import createElementFromStr
import unittest

class SlhaDecompositionTest(unittest.TestCase):
    from smodels.tools.smodelsLogging import logger

    def testSimple(self):
        self.logger.info ( "test a simple decomposition, no compression" )
        """ test the decomposition with no compression """
         
        #Load particles
        f = open("particleDefinitions.pcl","rb")
        modelParticles = pickle.load(f)
        f.close()
        particlesDict = dict([[p._name,p] for p in modelParticles])
         
        slhafile="../inputFiles/slha/simplyGluino.slha"
        xSectionDict,particlesList = getInputData(slhafile,modelParticles)
        topos = decomposer.decompose(xSectionDict,particlesList, .5*fb, False, False, 5.*GeV)        
        self.assertEqual( len(topos), 1 )
        #print len(topos),"topologies."
        topo=topos[0]
        #print topo
        ellist=topo.elementList
        self.assertEqual( len(ellist), 3 )
        self.assertEqual( str(ellist[0]), "[[[u*,u]],[[u*,u]]]" )
        elStr = createElementFromStr("[[[jet,jet]],[[jet,jet]]]",particlesDict)
        for el in ellist:
            self.assertEqual(el,elStr)
              
        gluinoxsec = 572.1689*fb
        br = 0.5
        self.assertAlmostEqual(ellist[0].weight[0].value.asNumber(fb),
                               (gluinoxsec*br*br).asNumber(fb),3)
        self.assertAlmostEqual(ellist[1].weight[0].value.asNumber(fb),
                               2*(gluinoxsec*br*br).asNumber(fb),3)
        self.assertAlmostEqual(ellist[2].weight[0].value.asNumber(fb),
                               (gluinoxsec*br*br).asNumber(fb),3)
        
    def testComplex(self):
        self.logger.info ( "test a complex decomposition, no compression" )
        """ test the decomposition with no compression """
         
        #Load particles
        f = open("particleDefinitions.pcl","rb")
        modelParticles = pickle.load(f)
        f.close()
        particlesDict = dict([[p._name,p] for p in modelParticles])
         
        slhafile="../inputFiles/slha/lightEWinos.slha"
        xSectionDict,particlesList = getInputData(slhafile,modelParticles)
        topos = decomposer.decompose(xSectionDict,particlesList, .5*fb, False, False, 5.*GeV)        
        self.assertEqual(len(topos), 17)
        self.assertEqual(len(topos.getElements()), 364)
        #print len(topos),"topologies."
        el=topos[0].elementList[0]
        self.assertEqual(len(topos[0].elementList), 1)
        self.assertEqual(len(topos[-1].elementList), 64)
        elStr = createElementFromStr("[[],[[Z]]]",particlesDict)
        self.assertEqual(el, elStr)
        el = topos.getElements()[-1]
        elStr = createElementFromStr("[[[b+,t+],[c*,s]],[[g],[W+],[c*,s]]]",particlesDict)       
        self.assertEqual(el, elStr)
        w = 2.*0.86*fb
        self.assertAlmostEqual(el.weight[0].value.asNumber(fb), w.asNumber(fb), 2)
        
    def testComplexCompression(self):
        self.logger.info ( "test a complex decomposition, with compression" )
        """ test the decomposition with no compression """
         
        #Load particles
        f = open("particleDefinitions.pcl","rb")
        modelParticles = pickle.load(f)
        f.close()        
        #Path to input file name (either a SLHA or LHE file)
        slhafile = '../inputFiles/slha/lightEWinos.slha'
        xSectionDict,particlesList = getInputData(slhafile,modelParticles)        
    
        #Set main options for decomposition:
        sigmacut = 0.3 * fb
        mingap = 5. * GeV    
        """ Decompose model (use slhaDecomposer for SLHA input or lheDecomposer for LHE input) """
        smstoplist = decomposer.decompose(xSectionDict,particlesList, sigmacut, doCompress=True, 
                                          doInvisible=True, minmassgap=mingap)
                
        topweights = [0.45069449699999997, 0.8060860593513246, 303.2711905718159, 
                      0.3637955190752928, 15.339097274505018, 3358.5703119038644, 
                      2.160128693731249, 4.5031659201699235, 9.486866215694839, 
                      17.980093334558738, 2.164988614601289, 147.823524958894, 
                      0.8102285793970079, 39.61814354617472, 16.436088764209956, 
                      9.291691094022376, 299.60965243794095, 46.34557917490008, 
                      15.658947303819541, 7.991672275037888, 158.71034751774343, 
                      56.78304329128415, 15.597472309381796]
       
        self.assertEqual(len(smstoplist), 23)        
        self.assertEqual(len(smstoplist.getElements()), 669)
        for itop,top in enumerate(smstoplist):
            self.assertAlmostEqual(top.getTotalWeight()[0].value.asNumber(fb), 
                                   topweights[itop],4)        
        

if __name__ == "__main__":
    unittest.main()
