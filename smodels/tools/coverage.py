#!/usr/bin/env python3

"""
.. module:: coverage
   :synopsis: Definitions of classes used to find, format missing topologies
    
.. moduleauthor:: Ursula Laa <ursula.laa@lpsc.in2p3.fr>    
.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>

"""
import copy
from smodels.tools.physicsUnits import fb

class Uncovered(object):
    """
    Object collecting all information of non-tested/covered elements
    :ivar topoList: sms topology list
    :ivar sumL: if true, sum up electron and muon to lepton, for missing topos
    :ivar sumJet: if true, sum up jets, for missing topos
    :ivar sqrts: Center of mass energy. If defined it will only consider cross-sections
    for this value. Otherwise the highest sqrts value will be used.
    """
    def __init__(self, topoList, sumL=True, sumJet=True, sqrts=None):
        
        
        if sqrts is None:
            self.sqrts = max([xsec.info.sqrts for xsec in topoList.getTotalWeight()])
        else:
            self.sqrts = sqrts
        self.missingTopos = UncoveredList(sumL, sumJet, self.sqrts)
        self.outsideGrid = UncoveredList(sumL, sumJet, self.sqrts) # FIXME change this to derived objects for printout
        self.longCascade = UncoveredClassifier()
        self.asymmetricBranches = UncoveredClassifier()
        self.motherIDs = []
        self.prevMothers = []
        self.outsideGridMothers = []
        self.getAllMothers(topoList)
        self.fill(topoList)
        self.asymmetricBranches.combine()
        self.longCascade.combine()

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
        for el in topoList.getElements(): # loop over all elements, by construction we start with the most compressed
            if self.inPrevMothers(el): missing = False # cannot be missing if element with same mothers has already appeared
            # this is because it can certainly be compressed further to the smaller element already seen in the loop
            else: #if not the case, we add mothers to previous mothers and test if the topology is missin
                self.addPrevMothers(el)
                missing = self.isMissingTopo(el) #missing topo only if not covered, and counting only weights of non-covered mothers
                # in addition, mother elements cannot be missing, we consider only the most compressed one
            if not missing: # any element that is not missing might be outside the grid
                # outside grid should be smalles covered but not tested in compression
                if el.covered and not el.tested: # verify first that element is covered but not tested
                    if not el.weight.getXsecsFor(self.sqrts): continue # remove elements that only have weight at higher sqrts
                    if self.inOutsideGridMothers(el): continue # if daughter element of current element is already counted skip this element
                    # in this way we do not double count, but always count the smallest compression that is outside the grid
                    outsideX = self.getOutsideX(el) # as for missing topo, recursively find untested cross section
                    if outsideX: # if all mothers are tested, this is no longer outside grid contribution, otherwise add to list
                        el.missingX =  outsideX # for combined printing function, call outside grid weight missingX as well
                        self.outsideGrid.addToTopos(el) # add to list of outsideGrid topos
                continue
            self.missingTopos.addToTopos(el) #keep track of all missing topologies
            if self.hasLongCascade(el): self.longCascade.addToClasses(el)
            elif self.hasAsymmetricBranches(el): self.asymmetricBranches.addToClasses(el) # if no long cascade, check for asymmetric branches

    def inPrevMothers(self, el): #check if smaller element with same mother has already been checked
        for mEl in el.motherElements:
            if mEl[-1].elID in self.prevMothers: return True
        return False

    def inOutsideGridMothers(self, el): #check if this element or smaller element with same mother has already been checked
        if el.elID in self.outsideGridMothers: return True
        for mEl in el.motherElements:
            if mEl[-1].elID in self.outsideGridMothers: return True
        return False

    def addPrevMothers(self, el): #add mother elements of currently tested element to previous mothers
        for mEl in el.motherElements:
            self.prevMothers.append(mEl[-1].elID)

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
        if el.elID in self.motherIDs: return False # mother element can not be missing
        if el.covered: return False # covered = not missing
        missingX = self.getMissingX(el) # find total missing cross section by checking if mothers are covered
        if not missingX: return False # if all mothers covered, element is not missing
        el.missingX = missingX # missing cross section is found by adding up cross section of mothers not covered
        return True

    def getMissingX(self,el):
        """
        Calculate total missing cross section of element, by recursively checking if mothers are covered
        :ivar el: Element
        :returns: missing cross section in fb as number
        """
        mothers = el.motherElements
        alreadyChecked = [] # for sanity check
        # if element has no mothers, the full cross section is missing
        if not el.weight.getXsecsFor(self.sqrts): return 0.
        missingX =  el.weight.getXsecsFor(self.sqrts)[0].value.asNumber(fb)
        if not mothers: return missingX
        while mothers: # recursive loop to check all mothers
            newmothers = []
            for mother in mothers:
                if mother[-1].elID in alreadyChecked: continue # sanity check, to avoid double counting
                alreadyChecked.append(mother[-1].elID) # now checking, so adding to alreadyChecked
                if mother[-1].covered:
                    if not mother[-1].weight.getXsecsFor(self.sqrts): continue
                    missingX -= mother[-1].weight.getXsecsFor(self.sqrts)[0].value.asNumber(fb)
                    continue # do not count cross section if mother is covered, do not continue recursion for this contribution
                if not mother[-1].motherElements: continue # end of recursion if element has no mothers, we keep its cross section in missingX
                else: newmothers += mother[-1].motherElements # if element has mother element, check also if those are covered before adding the cross section contribution
            mothers = newmothers # all new mothers will be checked until we reached the end of all recursions
        return missingX

    def getOutsideX(self,el):
        """
        Calculate total outside grid cross section of element, by recursively checking if mothers are covered
        :ivar el: Element
        :returns: missing cross section in fb as number
        """
        #same as getMissingX, but we also keep track of the mother elements of outsideGrid contributions
        # this is so we can find the smallest covered but not tested element in a chain of compressed elements
        mothers = el.motherElements
        alreadyChecked = []
        if not el.weight.getXsecsFor(self.sqrts): return 0.
        missingX =  el.weight.getXsecsFor(self.sqrts)[0].value.asNumber(fb)
        if not mothers: return missingX
        while mothers:
            newmothers = []
            for mother in mothers:
                if mother[-1].elID in alreadyChecked: continue # sanity check, to avoid double counting
                alreadyChecked.append(mother[-1].elID) # now checking, so adding to alreadyChecked
                if mother[-1].tested:
                    if not mother[-1].weight.getXsecsFor(self.sqrts): continue
                    missingX -= mother[-1].weight.getXsecsFor(self.sqrts)[0].value.asNumber(fb)
                    continue
                self.outsideGridMothers.append(mother[-1].elID) # mother element is not tested, but should no longer be considered as outside grid, because we already count its contribution here
                if not mother[-1].motherElements: continue
                else: newmothers += mother[-1].motherElements
            mothers = newmothers
        return missingX



    def getMissingXsec(self, sqrts=None):
        """
        Calculate total missing topology cross section at sqrts. If no sqrts is given use self.sqrts
        :ivar sqrts: sqrts
        """
        xsec = 0.
        if not sqrts: sqrts = self.sqrts
        for topo in self.missingTopos.topos:
            for el in topo.contributingElements:
                xsec += el.missingX
        return xsec

    def getOutOfGridXsec(self, sqrts=None): #FIXME same as getMissingXsec but different object, should not be separate functions
        xsec = 0.
        if not sqrts: sqrts = self.sqrts
        for topo in self.outsideGrid.topos:
            for el in topo.contributingElements:
                xsec += el.missingX
        return xsec

    def getLongCascadeXsec(self, sqrts=None):
        xsec = 0.
        if not sqrts: sqrts = self.sqrts
        for uncovClass in self.longCascade.classes:
            for el in uncovClass.contributingElements:
                xsec += el.missingX
        return xsec

    def getAsymmetricXsec(self, sqrts=None):
        xsec = 0.
        if not sqrts: sqrts = self.sqrts
        for uncovClass in self.asymmetricBranches.classes:
            for el in uncovClass.contributingElements:
                xsec += el.missingX
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
        motherPIDs = self.getMotherPIDs(el)
        for entry in self.classes:
            if entry.add(motherPIDs, el): return
        self.classes.append(UncoveredClass(motherPIDs, el))

    def getMotherPIDs(self, el):
        allPIDs = []
        for pids in el.getMothers():
            cPIDs = []
            for pid in pids:
                cPIDs.append(abs(pid))
            cPIDs.sort()
            if not cPIDs in allPIDs:
                allPIDs.append(cPIDs)
        allPIDs.sort()
        return allPIDs

    def combine(self):
        for ecopy in copy.deepcopy(self.classes):
            for e in self.classes:
                if e.isSubset(ecopy):
                    e.combine(ecopy)
                    self.remove(ecopy)

    def remove(self, cl):
        """
        Remove element where mother pids match exactly
        """
        for i, o in enumerate(self.classes):
            if o.motherPIDs == cl.motherPIDs:
                del self.classes[i]
                break

    def getSorted(self,sqrts):
        """
        Returns list of UncoveredClass objects in self.classes, sorted by weight
        :ivar sqrts: sqrts for weight lookup
        """
        return sorted(self.classes, key=lambda x: x.getWeight(), reverse=True)

class UncoveredClass(object):
    """
    Object collecting all elements contributing to the same uncovered class, defined by the mother PIDs.
    :ivar motherPIDs: PID of initially produces particles, sorted and without charge information
    :ivar el: Element
    """
    def __init__(self, motherPIDs, el):
        self.motherPIDs = motherPIDs # holds nested list of mother PIDs as given by element.getMothers
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

    def combine(self, other):
        for el in other.contributingElements:
            self.contributingElements.append(el)

    def getWeight(self):
        """
        Calculate weight at sqrts
        :ivar sqrts: sqrts
        """
        xsec = 0.
        for el in self.contributingElements:
            xsec += el.missingX     
        return xsec

    def isSubset(self, other):
        """
        True if motherPIDs of others are subset of the motherPIDs of this UncoveredClass
        """
        if len(other.motherPIDs) >= len(self.motherPIDs): return False
        for mothers in other.motherPIDs:
            if not mothers in self.motherPIDs: return False
        return True
  
class UncoveredTopo(object):
    """
    Object to describe one missing topology result / one topology outside the mass grid
    :ivar topo: topology description
    :ivar weights: weights dictionary
    """
    def __init__(self, topo, contributingElements=[]):
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
        name = name + ' (%s)'%(str(el.getFinalStates()).replace('[','').replace(']',''))
        for topo in self.topos:
            if name == topo.topo:
                topo.contributingElements.append(el)
                return
        self.topos.append(UncoveredTopo(name, [el]))
        return

    def generalName(self, instr):
        """
        generalize by summing over charges
        e, mu are combined to l
        :parameter instr: element as string
        :returns: string of generalized element
        """
        
        # 180318 mat: BUG? #############################
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

