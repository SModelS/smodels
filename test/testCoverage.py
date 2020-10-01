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
from smodels.theory.topology import TopologyList, Topology
from smodels.tools.physicsUnits import GeV, TeV, fb
from smodels.theory.element import Element
from smodels.theory.crossSection import XSectionList, XSection, XSectionInfo
from smodels.share.models import SMparticles, mssm
from smodels.tools.reweighting import reweightFactorFor
from smodels.theory.particle import Particle
import random

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
        factors[key] = lambda el: reweightFactorFor(el,'prompt')+reweightFactorFor(el,'displaced')


class CoverageTest(unittest.TestCase):

    def testUncovered(self):

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
        self.assertAlmostEqual(MET.generalElements[1].missingX, 0.03890665682,places=4)


    def testUncoveredTree(self):

        p1 = Particle(Z2parity=-1, label='p1', pdg=10001, mass=100.*GeV,
              eCharge=0., colordim=1, totalwidth=0*GeV)
        p2 = Particle(Z2parity=-1, label='p2', pdg=10002, mass=200.*GeV,
              eCharge=0., colordim=1, totalwidth=1e-15*GeV)
        xsecA = XSection()
        xsecA.value = 10.*fb
        xsecA.info.sqrts = 13.*TeV
        xsecA.info.label = 'wA'
        xsecA.info.order = 0
        xsecB = XSection()
        xsecB.value = 5.*fb
        xsecB.info.sqrts = 13.*TeV
        xsecB.info.label = 'wb'
        xsecB.info.order = 0

        #Element family-tree: A0->A1+B1, A1->A2
        #Mother A
        elA = Element()
        elA.label = 'A0'
        elA.testedBy = ['prompt','displaced']
        elA.coveredBy = ['prompt','displaced']
        elA.weight.add(xsecA)
        elA.motherElements = [elA]

        #Daughters:
        elA1 = Element()
        elA1.label = 'A1'
        elA1.testedBy = []
        elA1.coveredBy = []
        elA1.weight.add(xsecA)
        elA1.motherElements = [elA1,elA]
        elB1 = Element()
        elB1.label = 'B1'
        elB1.testedBy = []
        elB1.coveredBy = []
        elB1.weight.add(xsecA)
        elB1.motherElements = [elB1,elA]
        elA2 = Element()
        elA2.label = 'A2'
        elA2.testedBy = []
        elA2.coveredBy = []
        elA2.weight.add(xsecA)
        elA2.motherElements = [elA2,elA1]

        #Element family-tree: a0->a1+b1
        #Mother B
        ela = Element()
        ela.label = 'a0'
        ela.testedBy = []
        ela.coveredBy = []
        ela.weight.add(xsecB)
        ela.motherElements = [ela]
        #Daughters:
        ela1 = Element()
        ela1.label = 'a1'
        ela1.testedBy = []
        ela1.coveredBy = []
        ela1.weight.add(xsecB)
        ela1.motherElements = [ela1,ela]
        elb1 = Element()
        elb1.label = 'b1'
        elb1.testedBy = []
        elb1.coveredBy = []
        elb1.weight.add(xsecB)
        elb1.motherElements = [elb1,ela]

        #Merged element = (A2+b1)
        elCombined = Element()
        xsec = XSection()
        xsec.value = 15.*fb
        xsec.info.sqrts = 13.*TeV
        xsec.info.label = 'wA+wB'
        xsec.info.order = 0
        elCombined.label = 'A2+b1'
        elCombined.testedBy = []
        elCombined.coveredBy = []
        elCombined.weight.add(xsec)
        elCombined.motherElements = [elCombined,elA2,elb1]

        topoList = TopologyList()
        topoList.topos.append(Topology())
        elList = [elA,elA1,elA2,elB1,ela,elb1,ela1,elCombined]
        #Set odd particles (important for sorting the elements)
        for el in elList:
            for branch in el.branches:
                branch.oddParticles = [p1]
        elB1.branches[0].oddParticles = [p2,p1]
        ela1.branches[1].oddParticles = [p2,p1]
        elb1.branches[1].oddParticles = [p2,p1]
        ela.branches[0].oddParticles = [p2,p2,p1]

        #make sure the ordering in elList is not important:
        random.shuffle(elList)
        topoList.topos[0].elementList = elList[:]
        topoList._setElementIds()

        #Test if the family tree is being retrieved correctly and in the correct ordering (mother before grandmother,...):
        elListAncestors = {'A0' : [], 'A1' : ['A0'], 'A2' : ['A1','A0'], 'B1' : ['A0'],'A2+b1' : ['A2','b1','A1','a0','A0'],
                           'b1' : ['a0'], 'a1' : ['a0'] , 'a0' : []}
        for el in elList:
            ancestors = [mom.label for mom in el.getAncestors()]
            self.assertEqual(ancestors,elListAncestors[el.label])
        # A2+b1--> b1 is not tested, A2 is tested because its grandmother is tested
        missingTopos = Uncovered(topoList).getGroup('missing (all)')
        self.assertEqual(len(missingTopos.generalElements),1) #Only elCombined should appear (it is smaller than a1)
        self.assertAlmostEqual(missingTopos.generalElements[0].missingX,5.) #Only the fraction of the cross-section from elB is not missing
        self.assertEqual(len(missingTopos.generalElements[0]._contributingElements),1) #Only elCombined should appear
        self.assertTrue(missingTopos.generalElements[0]._contributingElements[0] is elCombined)

if __name__ == "__main__":
    unittest.main()
