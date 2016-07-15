#!/usr/bin/env python

"""
.. module:: tools.missingTopologies
   :synopsis: Definitions of classes used to find, format missing topologies
    
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>    
.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>

"""
from smodels.tools.physicsUnits import fb
import logging, sys
log = logging.getLogger(__name__)

class Uncovered(object):
    """
    Object collecting all information of non-tested/covered elements
    :ivar topoList: sms topology list
    :ivar sumL: if true, sum up electron and muon to lepton, for missing topos
    :ivar sumJet: if true, sum up jets, for missing topos
    """
    def __init__(self, topoList, sumL=True, sumJet=True):
        self.sqrts = max([xsec.info.sqrts for xsec in topoList.getTotalWeight()])
        self.missingTopos = UncoveredList(sumL, sumJet, self.sqrts)
        self.outsideGrid = UncoveredList(sumL, sumJet, self.sqrts) # FIXME change this to derived objects for printout
        self.longCascade = UncoveredClassifier()
        self.asymmetricBranches = UncoveredClassifier()
        self.motherIDs = []
        self.getAllMothers(topoList)
        self.fill(topoList)

    def getAllMothers(self, topoList):
        """
        Find all IDs of mother elements, only most compressed element can be missing topology
        :ivar topoList: sms topology list
        """
        for el in topoList.getElements():
            for mEl in el.motherElements:
                motherID = mEl[-1].elID
                if not motherID in self.motherIDs: self.motherIDs.append(motherID)

    def fill(self, topoList):
        """
        Check all elements, categorise those not tested / missing, classify long cascade decays and asymmetric branches
        Fills all corresponding objects
        :ivar topoList: sms topology list
        """
        for el in topoList.getElements():
            missing = self.isMissingTopo(el) #missing topo only if not mother, covered, mother covered
            if not missing: # if not missing check if it is actually tested
                if not el.tested:
                    self.outsideGrid.addToTopos(el) # not missing but not tested means we are outside the mass grid
                continue
            self.missingTopos.addToTopos(el) #keep track of all missing topologies
            if self.hasLongCascade(el): self.longCascade.addToClasses(el)
            elif self.hasAsymmetricBranches(el): self.asymmetricBranches.addToClasses(el) # if no long cascade, check for asymmetric branches

    def hasLongCascade(self, el):
        """
        Return True if element has more than 3 particles in the decay chain
        :ivar el: Element
        """
        if el._getLength() > 3: return True
        return False

    def hasAsymmetricBranches(self, el):
        """
        Return True if Element branches are not equal
        :ivar el: Element
        """
        if el.branches[0] == el.branches[1]: return False
        return True

    def isMissingTopo(self, el):
        """
        A missing topology is not a mother element, not covered, and does not have mother which is covered
        :ivar el: Element
        """
        if el.elID in self.motherIDs: return False
        if el.covered: return False
        if self.motherCovered(el): return False
        return True

    def motherCovered(self, el):
        """
        Recursively check if a mother of the given Element is covered
        :ivar el: Element
        """
        mothers = el.motherElements
        while mothers:
            newmothers = []
            for mother in mothers:
                if mother[-1].covered: return True
                newmothers += mother[-1].motherElements
            mothers = newmothers
        return False

    def getMissingXsec(self, sqrts=None):
        """
        Calculate total missing topology cross section at sqrts. If no sqrts is given use self.sqrts
        :ivar sqrts: sqrts
        """
        xsec = 0.
        if not sqrts: sqrts = self.sqrts
        for topo in self.missingTopos.topos:
            for el in topo.contributingElements:
                if not el.weight.getXsecsFor(sqrts): continue
                xsec += el.weight.getXsecsFor(sqrts)[0].value.asNumber(fb)
        return xsec

    def getOutOfGridXsec(self, sqrts=None): #FIXME same as getMissingXsec but different object, should not be separate functions
        xsec = 0.
        if not sqrts: sqrts = self.sqrts
        for topo in self.outsideGrid.topos:
            for el in topo.contributingElements:
                if not el.weight.getXsecsFor(sqrts): continue
                xsec += el.weight.getXsecsFor(sqrts)[0].value.asNumber(fb)
        return xsec

    def getLongCascadeXsec(self, sqrts=None):
        xsec = 0.
        if not sqrts: sqrts = self.sqrts
        for uncovClass in self.longCascade.classes:
            for el in uncovClass.contributingElements:
                if not el.weight.getXsecsFor(sqrts): continue
                xsec += el.weight.getXsecsFor(sqrts)[0].value.asNumber(fb)
        return xsec

    def getAsymmetricXsec(self, sqrts=None):
        xsec = 0.
        if not sqrts: sqrts = self.sqrts
        for uncovClass in self.asymmetricBranches.classes:
            for el in uncovClass.contributingElements:
                if not el.weight.getXsecsFor(sqrts): continue
                xsec += el.weight.getXsecsFor(sqrts)[0].value.asNumber(fb)
        return xsec
          
class UncoveredClassifier(object):
    """
    Object collecting elements with long cascade decays or asymmetric branches.
    Objects are grouped according to the initially produced particle PID pair.
    """
    def __init__(self):
        self.classes = []

    def addToClasses(self, el):
        """
        Add Element in corresponding UncoveredClass, defined by mother PIDs.
        If no corresponding class in self.classes, add new UncoveredClass
        :ivar el: Element
        """
        motherPIDs = [abs(pid) for pid in el.getMothers()[0]]
        motherPIDs.sort()
        for entry in self.classes:
            if entry.add(motherPIDs, el): return
        self.classes.append(UncoveredClass(motherPIDs, el))

    def getSorted(self,sqrts):
        """
        Returns list of UncoveredClass objects in self.classes, sorted by weight
        :ivar sqrts: sqrts for weight lookup
        """
        return sorted(self.classes, key=lambda x: x.getWeight(sqrts).asNumber(fb), reverse=True)

class UncoveredClass(object):
    """
    Object collecting all elements contributing to the same uncovered class, defined by the mother PIDs.
    :ivar motherPIDs: PID of initially produces particles, sorted and without charge information
    :ivar el: Element
    """
    def __init__(self, motherPIDs, el):
        self.motherPIDs = motherPIDs # holds list of mother PIDs as given by element.getMothers
        self.contributingElements = [el] # collect all contributing elements, to keep track of weights as well
    def add(self, motherPIDs, el):
        """
        Add Element to this UncoveredClass object if motherPIDs match and return True, else return False
        :ivar motherPIDs: PID of initially produces particles, sorted and without charge information
        :ivar el: Element
        """
        if not motherPIDs == self.motherPIDs: return False
        self.contributingElements.append(el)
        return True
    def getWeight(self, sqrts):
        """
        Calculate weight at sqrts
        :ivar sqrts: sqrts
        """
        xsec = 0.*fb
        for el in self.contributingElements:
            elxsec = el.weight.getXsecsFor(sqrts)            
            if not elxsec: continue
            xsec += elxsec[0].value            
        return xsec

  
class UncoveredTopo(object):
    """
    Object to describe one missing topology result / one topology outside the mass grid
    :ivar topo: topology description
    :ivar weights: weights dictionary
    """
    def __init__(self, topo, weights, contributingElements=[]):
        self.topo = topo
        self.contributingElements = contributingElements
        self.value = 0. # weight for sqrts set in uncoveredList, only set this in printout

class UncoveredList(object):
    """
    Object to find and collect UncoveredTopo objects, plus printout functionality
    :ivar sumL: if true sum electrons and muons to leptons
    :ivar sumJet: if true, sum up jets
    :ivar sqrts: sqrts, for printout
    """
    def __init__(self, sumL, sumJet, sqrts):
        self.topos = []
        self.sumL = sumL
        self.sumJet = sumJet
        self.sqrts = sqrts

    def addToTopos(self, el):
        """
        adds an element to the list of missing topologies
        if the element contributes to a missing topology that is already
        in the list, add weight to topology
        :parameter el: element to be added
        """
        name = self.orderbranches(self.generalName(el.__str__()))
        for topo in self.topos:
            if name == topo.topo:
                topo.contributingElements.append(el)
                return
        self.topos.append(UncoveredTopo(name, el.weight, [el]))
        return

    def generalName(self, instr):
        """
        generalize by summing over charges
        e, mu are combined to l
        :parameter instr: element as string
        :returns: string of generalized element
        """
        from smodels.theory.particleNames import ptcDic
        if self.sumL: exch = ["W", "l", "t", "ta"]
        else: exch = ["W", "e", "mu", "t", "ta"]
        if self.sumJet: exch.append("jet")
        for pn in exch:
            for on in ptcDic[pn]:
                instr = instr.replace(on, pn).replace("hijetjets","higgs")
        return instr

    def orderbranches(self, instr):
        """
        unique ordering of branches
        :parameter instr: element as string
        :returns: string of ordered element
        """
        from smodels.theory.element import Element
        li = Element(instr).getParticles()
        for be in li:
            for ve in be:
                ve.sort()
        li.sort()
        return str(li).replace("'", "").replace(" ", "")

