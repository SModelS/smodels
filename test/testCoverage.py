#!/usr/bin/env python3

"""
.. module:: testReweighting
   :synopsis: Tests the function of coverage 
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>   
"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.theory.branch import Branch
from smodels.tools.coverage import Uncovered
from smodels.theory.topology import TopologyList
from smodels.tools.physicsUnits import GeV, TeV, fb
from smodels.theory.element import Element
from smodels.theory.crossSection import XSectionList, XSection, XSectionInfo
from smodels.share.models import SMparticles, mssm
from smodels.tools.reweighting import reweightFactorFor

n1 = mssm.n1
n1.totalwidth = 0.*GeV
st1 = mssm.st1
st1.totalwidth = 2.*GeV
st2 = mssm.st2
st2.totalwidth = 10**(-15)*GeV
gluino = mssm.gluino
gluino.totalwidth = 1.*10**(-30)*GeV
t = SMparticles.t
b = SMparticles.b


#Default definitions for the uncovered topology categories/groups:
##Element filters for each group:
##(it should be a function which takes an Element object as input
##and returns True if the element belongs to the group and False otherwise)
filters = {'missing MET (prompt)' : lambda el: (not ('prompt' in el.coveredBy)) and (el.branches[0].oddParticles[-1].isMET()
                                                                                     and el.branches[1].oddParticles[-1].isMET()),
            'missing non-MET (prompt)' : lambda el: (not ('prompt' in el.coveredBy)) and (not el.branches[0].oddParticles[-1].isMET()
                                                                                     or not el.branches[1].oddParticles[-1].isMET()),
                'missing (displaced)' : lambda el: not ('displaced' in el.coveredBy)}

##Weight factors for each group:
##(it should be a function which takes an Element object as input
##and returns the reweighting factor to be applied to the element weight. It is relevant if only
##the fraction of the weight going into prompt or displaced decays is required)
factors = {}
for key in filters:
    if 'prompt' in key.lower():
        factors[key] = lambda el: reweightFactorFor(el,'prompt')
    elif 'displaced' in key.lower():
        factors[key] = lambda el: reweightFactorFor(el,'displaced')
    else:
        factors[key] = lambda el: 1.


class CoverageTest(unittest.TestCase):     
   
    def testFill(self):
       
        topolist = TopologyList()
        
        # prompt
        b1 = Branch()
        b1.evenParticles = [[b,t]]
        b1.oddParticles = [st1,n1]
        b1.setInfo()
        el1 = Element()
        el1.branches=[b1,b1]
        
        # long-lived
        b3 = Branch()
        b3.evenParticles = []
        b3.oddParticles = [gluino]
        b4 = Branch()
        b4.evenParticles = [[b,t]]
        b4.oddParticles = [st1,n1]
        b3.setInfo()
        b4.setInfo()
        el2 = Element()
        el2.branches=[b3,b4]
        
        # prompt and displaced
        b5 = Branch()
        b5.evenParticles = [[t],[b,t]]
        b5.oddParticles = [st2,st1,n1]
        b6 = Branch()
        b6.evenParticles = [[b,t]]
        b6.oddParticles = [st1,n1]
        b5.setInfo()
        b6.setInfo()
        el3 = Element()
        el3.branches=[b5,b6]
        
        
        w1 = XSectionList()
        w1.xSections.append(XSection())
        w1.xSections[0].info = XSectionInfo()
        w1.xSections[0].info.sqrts = 8.*TeV
        w1.xSections[0].info.label = '8 TeV'
        w1.xSections[0].info.order = 0
        w1.xSections[0].value = 10.*fb    
        el1.weight = w1
        el2.weight = w1
        el3.weight = w1
        
        topolist.addElement(el1)
        topolist.addElement(el2)
        topolist.addElement(el3)
        topolist._setElementIds()
        uncovered = Uncovered(topolist,groupFilters=filters,groupFactors=factors)
        longLived = uncovered.getGroup('missing non-MET (prompt)')
       
        
        MET = uncovered.getGroup('missing MET (prompt)')
        displaced = uncovered.getGroup('missing (displaced)')
               
        self.assertEqual(len(longLived.generalElements), 1)
        self.assertEqual(len(displaced.generalElements), 1)
        self.assertEqual(len(MET.generalElements), 2)
       
        self.assertAlmostEqual(longLived.generalElements[0].missingX, 10.)
        self.assertAlmostEqual(displaced.generalElements[0].missingX, 9.96109334317542,places=4)
        self.assertAlmostEqual(MET.generalElements[0].missingX, 10.)
        self.assertAlmostEqual(MET.generalElements[1].missingX, 0.03890665682,places=6)

        
        
    
    
if __name__ == "__main__":
    unittest.main()       
