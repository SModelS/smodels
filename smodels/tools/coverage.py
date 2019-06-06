#!/usr/bin/env python3

"""
.. module:: coverage
   :synopsis: Definitions of classes used to find, group and format missing topologies
    
.. moduleauthor:: Ursula Laa <ursula.laa@lpsc.in2p3.fr>    
.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
"""

from smodels.tools.physicsUnits import fb
from smodels.tools.reweighting import reweightFactorFor
from smodels.theory.branch import Branch
from smodels.theory.particle import MultiParticle, ParticleList
from smodels.share.models.SMparticles import e,mu,ta,taC,eC,muC,W,WC,t,tC,q,c,g,pion,nu
from smodels.theory import particle
from smodels.theory.auxiliaryFunctions import index_bisect
from smodels.experiment.databaseParticles import MET,HSCP,RHadronG,RHadronQ

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

    def __init__(self,topoList,sumL=True, sumJet=True, sqrts=None,
                 groupFilters = {'outsideGrid (prompt)' : lambda el: ('prompt' in el.coveredBy) and not ('prompt' in el.testedBy),
                           'outsideGrid (displaced)' : lambda el: ('displaced' in el.coveredBy) and not ('displaced' in el.testedBy),
                            'missing (prompt)' : lambda el: not ('prompt' in el.coveredBy),
                            'missing (displaced)' : lambda el: not ('displaced' in el.coveredBy),
                            'missing (long cascade)' : lambda el: (not el.coveredBy) and el._getLength() > 3},
                 groupFactors = {'outsideGrid (prompt)' : lambda el: reweightFactorFor(el,'prompt'),
                           'outsideGrid (displaced)' : lambda el: reweightFactorFor(el,'displaced'),
                            'missing (prompt)' : lambda el: reweightFactorFor(el,'prompt'),
                            'missing (displaced)' : lambda el: reweightFactorFor(el,'displaced'),
                            'missing (long cascade)' : lambda el: 1.}):
        
        #Define multiparticles for conveniently grouping final states
        eList = MultiParticle('e' , [e,eC])
        muList = MultiParticle('mu', [mu,muC])
        taList = MultiParticle('ta', [ta,taC])
        lList = MultiParticle('l', [e,mu,eC,muC])
        WList = MultiParticle('W', [W,WC])
        tList = MultiParticle('t', [t,tC])
        jetList = MultiParticle('jet', [q,c,g,pion])
        nuList  = nu
        if sumL:
            particleGroups = [WList, lList, tList, taList, nuList]
        else:
            particleGroups = [WList, eList, muList,tList, taList, nuList]
        if sumJet:
            particleGroups.append(jetList)

        if sqrts is None:
            sqrts = max([xsec.info.sqrts for xsec in topoList.getTotalWeight()])
        else:
            sqrts = sqrts

        #Define global final states to label the General Elements:
        bsmFinalStates = [MET,HSCP,RHadronG,RHadronQ]
        
        self.topoList = topoList
        self.groups = []
        #Create each uncovered group and get the topologies from topoList
        for gLabel,gFilter in groupFilters.items():
            #Initialize the uncovered topology list:
            uncoveredTopos = UncoveredList(label=gLabel,elementFilter=gFilter,
                                           reweightFactor = groupFactors[gLabel],
                                           particleGroups=particleGroups,
                                           bsmFinalStates=bsmFinalStates,
                                           sqrts=sqrts)
            #Fill the list with the elements in topoList:
            uncoveredTopos.getToposFrom(topoList)
            self.groups.append(uncoveredTopos)

class UncoveredList(object):
    """
    Object to find and collect GeneralElement objects, plus printout functionality
    :ivar generalElements: missing elements, grouped by common general name using inclusive labels (e.g. jet)
    :ivar particleGroups: Lists of MultiParticles used to group final states.
    :ivar sqrts: sqrts, for printout
    """
    def __init__(self, label, elementFilter, reweightFactor, particleGroups, bsmFinalStates, sqrts):
        self.generalElements = []
        self.particleGroups = particleGroups
        self.bsmFinalStates = bsmFinalStates
        self.sqrts = sqrts
        self.label = label
        self.elementFilter = elementFilter
        self.reweightFactor = reweightFactor
        
    def __str__(self):
        return self.label

    def __repr__(self):
        return str(self)

    def getToposFrom(self,topoList):
        """
        Select the elements from topoList according to self.elementFilter
        and build GeneralElements from the selected elements.
        The GeneralElement weights corresponds to the missing cross-section
        with double counting from compressed elements already accounted for.
        """
        
        #First select all elements according to the filter (type of uncovered/missing topology):
        elementList = [el for el in topoList.getElements() if self.elementFilter(el)]
        
        #Get missing xsections including the reweight factor:
        missingXandEls = [[self.getMissingX(el)*self.reweightFactor(el),el] for el in elementList]
        #Sort according to largest missingX and then smallest size
        missingXandEls = sorted(missingXandEls, key = lambda pt: [pt[0],-pt[1]._getLength()], reverse=True)

        #Split lists of elements and missingX:
        missingXsecs = [pt[0] for pt in missingXandEls]
        elementList = [pt[1] for pt in missingXandEls]
        
        #Remove all elements which are related to each other in order to avoid double counting
        #(keep always the first appearance in the list, so we always keep the ones with largest missing xsec)
        elementListUnique = []
        missingXsecsUnique = []
        for i,element in enumerate(elementList):
            if any(element.isRelatedTo(el) for el in elementListUnique):
                continue
            elementListUnique.append(element)
            missingXsecsUnique.append(missingXsecs[i])


        #Now that we only have unique elements with their effective missing cross-sections
        #we create General Elements out of them
        for i,el in enumerate(elementListUnique):
            missingX = missingXsecsUnique[i]
            if not missingX:
                continue
            self.addToGeneralElements(el,missingX)

    def getMissingX(self,element):
        """
        Calculate total missing cross section of an element, by recursively checking if its
        mothers already appear in the list.
        :param element: Element object

        :returns: missing cross section without units (in fb)
        """

        ancestorList = element.getAncestors()
        alreadyChecked = [] # for sanity check
        #If the element has no xsec for the required sqrts, return 0.
        if not element.weight.getXsecsFor(self.sqrts):
            return 0.
        #Get the total element weight:
        missingX =  element.weight.getXsecsFor(self.sqrts)[0].value
        overlapXsec = 0.*fb
        for ancestor in ancestorList: # recursive loop to check all mothers
            #Skip entries which correspond to the element itself
            if ancestor is element:
                continue
            #Make sure we do not subtract the same mother twice
            if any(ancestor is el for el in alreadyChecked):
                continue
            alreadyChecked.append(ancestor)
            #Check if mother passes the group filter (if it is missing or not):
            if self.elementFilter(ancestor):
                continue
            #Subtract the weight of the ancestor and skip checking for all older ancestors
            #(avoid subtracting the same weight twice from mother and grandmother for instance)
            alreadyChecked += ancestor.getAncestors()
            ancestorXsec = ancestor.weight.getXsecsFor(self.sqrts)[0]
            if not ancestorXsec:
                continue
            overlapXsec += ancestorXsec.value

        missingX -= overlapXsec

        return missingX.asNumber(fb)
        
        
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
            

    def addToGeneralElements(self, el, missingX):
        """
        Adds an element to the list of missing topologies = general elements.
        If the element contributes to a missing topology that is already in the list, add element and weight to topology.
        :parameter el: element to be added
        :parameter missingX: missing cross-section for the element (in fb)
        """     

        newGenEl = GeneralElement(el,missingX,self.particleGroups,self.bsmFinalStates)

        index = index_bisect(self.generalElements,newGenEl)
        if index != len(self.generalElements) and self.generalElements[index] == newGenEl:
            self.generalElements[index]._contributingElements.append(el)
            self.generalElements[index].missingX += missingX
        else:
            self.generalElements.insert(index,newGenEl)


class GeneralElement(object):
    """
    This class represents a simplified (general) element which does
    only holds information about its even particles and decay type.
    The even particles are replaced/grouped by the particles defined in particleGroups.
    """

    def __init__(self,el,missingX,particleGroups,bsmFinalStates):

        self.branches = [Branch() for _ in el.branches]
        self.missingX = missingX
        self.finalBSMstates = []

        for ib,branch in enumerate(el.branches):
            newParticles = []
            for vertex in branch.evenParticles:
                newVertex = vertex[:]
                for ip,particle in enumerate(vertex):
                    for particleList in particleGroups:
                        if particleList.contains(particle):
                            newVertex[ip] = particleList
                newParticles.append(ParticleList(newVertex))
            self.branches[ib].evenParticles = newParticles

            finalBSM = branch.oddParticles[-1]
            for bsmFS in bsmFinalStates:
                if finalBSM == bsmFS:
                    finalBSM = bsmFS
                    break
            self.branches[ib].finalBSMstate = finalBSM
            self.branches[ib].setInfo()

        self.sortBranches()
        self.finalBSMstates = [br.finalBSMstate for br in self.branches]
        self._contributingElements = [el]
        self._outputDescription = str(self)

    def __cmp__(self,other):
        """
        Compares the element with other for any branch ordering.
        The comparison is made based on branches.
        OBS: The elements and the branches must be sorted!
        :param other:  element to be compared (GeneralElement object)
        :return: -1 if self < other, 0 if self == other, +1, if self > other.
        """

        if not isinstance(other,GeneralElement):
            return -1

        #Compare branches:
        if self.branches != other.branches:
            comp = (self.branches > other.branches)
            if comp:
                return 1
            else:
                return -1
        #Compare decay type
        if self.finalBSMstates != other.finalBSMstates:
            comp = (self.finalBSMstates > other.finalBSMstates)
            if comp:
                return 1
            else:
                return -1

        return 0

    def __eq__(self,other):
        return self.__cmp__(other)==0

    def __lt__(self,other):
        return self.__cmp__(other)<0

    def __str__(self):
        """
        Create the element bracket notation string, e.g. [[[jet]],[[jet]],
        including the decay type.

        :returns: string representation of the element (in bracket notation)
        """

        elStr = "["+",".join([str(br) for br in self.branches])+"]"
        elStr = elStr.replace(" ", "").replace("'", "")
        name = elStr.replace('~','') + '  (%s)'%(','.join([str(p) for p in self.finalBSMstates]))
        return name

    def __repr__(self):

        return self.__str__()

    def sortBranches(self):
        """
        Sort branches. The smallest branch is the first one.
        If branches are equal, sort accoding to decayType.
        """

        #Now sort branches
        self.branches = sorted(self.branches, key = lambda br: (br,br.finalBSMstate))

