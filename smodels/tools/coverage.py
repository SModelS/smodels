#!/usr/bin/env python3

"""
.. module:: coverage
   :synopsis: Definitions of classes used to find, format missing topologies
    
.. moduleauthor:: Ursula Laa <ursula.laa@lpsc.in2p3.fr>    
.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
"""

from smodels.tools.physicsUnits import fb
from smodels.tools.reweighting import addPromptAndDisplaced
from smodels.theory.element import Element
import smodels.experiment.finalStateParticles as fS

class Uncovered(object):
    """
    Object collecting all information of non-tested/covered elements
    :ivar topoList: sms topology list
    :ivar sigmacut: if defined, will only keep elements with weight above sigmacut.
    :ivar sumL: if true, sum up electron and muon to lepton, for missing topos
    :ivar sumJet: if true, sum up jets, for missing topos
    :ivar sqrts: Center of mass energy. If defined it will only consider cross-sections
                for this value. Otherwise the highest sqrts value will be used.
    """
    def __init__(self, topoList, sigmacut=None, sumL=True, sumJet=True, sqrts=None):
        
        
        if sqrts is None:
            self.sqrts = max([xsec.info.sqrts for xsec in topoList.getTotalWeight()])
        else:
            self.sqrts = sqrts
        self.missingTopos = UncoveredList(sumL, sumJet, self.sqrts)
        self.outsideGrid = UncoveredList(sumL, sumJet, self.sqrts) 
        self.longLived = UncoveredList(sumL, sumJet, self.sqrts)
        self.displaced = UncoveredList(sumL, sumJet, self.sqrts)
        self.MET = UncoveredList(sumL, sumJet, self.sqrts)
        
        self.motherIDs = []
        self.prevMothers = []
        self.outsideGridMothers = []
        self.getAllMothers(topoList)
        self.fill(topoList, sigmacut)

    def getAllMothers(self, topoList):
        """
        Find all IDs of mother elements, only most compressed element can be missing topology
        :ivar topoList: sms topology list
        """
        
        for el in topoList.getElements():
            for mEl in el.motherElements:
                motherID = mEl[-1].elID
                if not el.elID == motherID and not motherID in self.motherIDs:                     
                    self.motherIDs.append(motherID)
                

    def fill(self, topoList, sigmacut):
        """
        Check all elements, categorise those not tested / missing, classify long cascade decays and asymmetric branches
        Fills all corresponding objects
        :ivar topoList: sms topology list
        :ivar sigmacut: if defined, will only keep elements with weight above sigmacut.
        """
        
        for element in topoList.getElements(): # loop over all elements, by construction we start with the most compressed
            allElements = []       
            probabilities1, branches1 = addPromptAndDisplaced(element.branches[0])
            probabilities2, branches2 = addPromptAndDisplaced(element.branches[1])               
            for i,probability1 in enumerate(probabilities1):
                for j,probability2 in enumerate(probabilities2): 
                    newEl = element.copy()   
                    newEl.tested = element.tested
                    newEl.covered = element.covered
                    newEl.branches[0]._decayType = branches1[i]._decayType
                    newEl.branches[1]._decayType = branches2[j]._decayType      
                    newEl.weight *= (probability1*probability2)    
                    if (not sigmacut is None) and newEl.weight.getMaxXsec() < sigmacut:
                        continue 
                    allElements.append(newEl) 
            
            
            for el in allElements:
                                
                if self.inPrevMothers(el): 
                    missing = False # cannot be missing if element with same mothers has already appeared
                # this is because it can certainly be compressed further to the smaller element already seen in the loop
                else: #if not the case, we add mothers to previous mothers and test if the topology is missing             
                    self.addPrevMothers(el)
                    missing = self.isMissingTopo(el) #missing topo only if not covered, and counting only weights of non-covered mothers
                    # in addition, mother elements cannot be missing, we consider only the most compressed one                                  
                if not missing: # any element that is not missing might be outside the grid      
                    # outside grid should be smalles covered but not tested in compression                  
                    if el.covered and not el.tested: # verify first that element is covered but not tested                              
                        if not el.weight.getXsecsFor(self.sqrts): continue # remove elements that only have weight at higher sqrts                    
                        if self.inOutsideGridMothers(el): continue # if daughter element of current element is already counted skip this element                  
                        # in this way we do not double count, but always count the smallest compression that is outside the grid
                        outsideX = self.getMissingX(el, checkFor = 'tested') # as for missing topo, recursively find untested cross section

                        if outsideX: # if all mothers are tested, this is no longer outside grid contribution, otherwise add to list
                            el.missingX =  outsideX # for combined printing function, call outside grid weight missingX as well
                            self.outsideGrid.addToGeneralElements(el) # add to list of outsideGrid topos                        
                    continue         
                self.missingTopos.addToGeneralElements(el) #keep track of all missing topologies                
                if self.hasDisplaced(el):
                    self.displaced.addToGeneralElements(el)
                elif self.hasLongLived(el):
                    self.longLived.addToGeneralElements(el)
                else:
                    self.MET.addToGeneralElements(el)

    def inPrevMothers(self, el): #check if smaller element with same mother has already been checked
        for mEl in el.motherElements:
            if mEl[0] != 'original' and mEl[-1].elID in self.prevMothers:
                return True
        return False

    def inOutsideGridMothers(self, el): #check if this element or smaller element with same mother has already been checked
        if el.elID in self.outsideGridMothers:
            return True
        for mEl in el.motherElements:
            if mEl[0] != 'original' and mEl[-1].elID in self.outsideGridMothers:
                return True
        return False

    def addPrevMothers(self, el): #add mother elements of currently tested element to previous mothers
        for mEl in el.motherElements:
            if mEl[0] != 'original' and mEl[-1].elID!=0: 
                self.prevMothers.append(mEl[-1].elID)
        
    def hasDisplaced(self, el):
        """
        Return True if Element has at least one displaced branch 
        :ivar el: Element
        """  
        if 'displaced' in el.branches[0]._decayType or 'displaced' in el.branches[1]._decayType: 
            return True
        else:
            return False                    
        
    def hasLongLived(self, el):
        """
        Return True if Element has at least one long lived branch but no displaced branch
        :ivar el: Element
        """  
        if el.branches[0]._decayType == 'longlived' or el.branches[1]._decayType == 'longlived': 
            return True 
        else:
            return False     
        
    def isMissingTopo(self, el):
        """
        A missing topology is not a mother element, not covered, and does not have a mother which is covered
        :ivar el: Element
        """     
        if el.elID in self.motherIDs: return False # mother element can not be missing        
        if el.covered: return False # covered = not missing
        missingX = self.getMissingX(el) # find total missing cross section by checking if mothers are covered
        if not missingX: return False # if all mothers covered, element is not missing
        el.missingX = missingX # missing cross section is found by adding up cross section of mothers not covered
        return True

    def getMissingX(self,el, checkFor = 'covered'):
        """
        Calculate total missing cross section of element, by recursively checking if mothers are covered
        :ivar el: Element
        :returns: missing cross section in fb as number
        """
        mothers = el.motherElements
        alreadyChecked = [] # for sanity check
        # if element has no mothers, the full cross section is missing
        if not el.weight.getXsecsFor(self.sqrts):
            return 0.
        missingX =  el.weight.getXsecsFor(self.sqrts)[0].value
        if not mothers:
            return missingX.asNumber(fb)
        coveredX = 0.*fb
        while mothers: # recursive loop to check all mothers
            newmothers = []
            for mother in mothers:
                if mother[-1].elID in alreadyChecked: continue # sanity check, to avoid double counting
                alreadyChecked.append(mother[-1].elID) # now checking, so adding to alreadyChecked
                
                if getattr(mother[-1], checkFor):
                    if not mother[-1].weight.getXsecsFor(self.sqrts): continue
                    coveredX += mother[-1].weight.getXsecsFor(self.sqrts)[0].value
                if checkFor=='tested':
                    if mother[-1].elID!=0: self.outsideGridMothers.append(mother[-1].elID) # mother element is not tested, but should no longer be considered as outside grid, because we already count its contribution here
                if all(grandmother[0] == 'original' for grandmother in mother[-1].motherElements) : continue # end of recursion if element has no mothers, we keep its cross section in missingX
                else: newmothers += mother[-1].motherElements # if element has mother element, check also if those are covered before adding the cross section contribution
            mothers = newmothers # all new mothers will be checked until we reached the end of all recursions
        missingX = missingX - coveredX
        return missingX.asNumber(fb)
    



class UncoveredList(object):
    """
    Object to find and collect UncoveredTopo objects, plus printout functionality
    :ivar generalElements: missing elements, grouped by common general name using inclusive labels (e.g. jet)
    :ivar sumL: if True sum electrons and muons to leptons
    :ivar sumJet: if True, sum up jets
    :ivar sqrts: sqrts, for printout
    """
    def __init__(self, sumL, sumJet, sqrts):
        self.generalElements = []
        self.sumL = sumL
        self.sumJet = sumJet
        self.sqrts = sqrts
        
        
    def getTotalXsec(self, sqrts=None):
        """
        Calculate total missing topology cross section at sqrts. If no sqrts is given use self.sqrts
        :ivar sqrts: sqrts
        """
        xsec = 0.
        if not sqrts: sqrts = self.sqrts
        for genEl in self.generalElements:
            xsec += genEl.missingX
        return xsec
            

    def addToGeneralElements(self, el):
        """
        Adds an element to the list of missing topologies = general elements.
        If the element contributes to a missing topology that is already in the list, add element and weight to topology.
        :parameter el: element to be added
        """     

        newGenEl = self.generalElement(el)
        newGenEl.sortBranches()

        name = str(newGenEl) + '  (%s)'%(str( newGenEl.branches[0]._decayType) + ","+ str(newGenEl.branches[1]._decayType) )    

        # Check if an element with the same generalized name has already been added
        for genEl in self.generalElements:
            if name == genEl._outputDescription:
                genEl._contributingElements.append(el)
                genEl.missingX += el.missingX
                return
            
        newGenEl._outputDescription = name           
        newGenEl._contributingElements = [el]    
        self.generalElements.append(newGenEl)
        

    def generalElement(self, el):
        """
        Creates a new element with a generalized name according to the inclusive labels. 
        Generalization is done by summing over charges:
        sumL: e, mu are combined to l
        sumJet: quarks, g, pi are combined to jet
        :parameter instr: element as string
        :returns: string of generalized element
        """

        newEl = Element()
        newEl.branches[0]._decayType = el.branches[0]._decayType
        newEl.branches[1]._decayType = el.branches[1]._decayType
        for ib,branch in enumerate(el.branches):
            newEl.branches[ib].oddParticles = branch.oddParticles[:]
        newEl.missingX = el.missingX
                 
        if self.sumL: exch = [fS.WList, fS.lList, fS.tList, fS.taList, fS.nuList] 
        else: exch = [fS.WList, fS.eList, fS.muList,fS.tList, fS.taList, fS.nuList]
        if self.sumJet: exch.append(fS.jetList)
          
        for ib,branch in enumerate(el.branches):
            newParticles = []
            for vertex in branch.evenParticles:
                newVertex = vertex[:]
                for ip,particle in enumerate(vertex):
                    for particleList in exch:
                        if particle == particleList:
                            newVertex[ip] = particleList
                newParticles.append(newVertex)
            newEl.branches[ib].evenParticles = newParticles
            newEl.branches[ib].setInfo()
        
        return newEl 

