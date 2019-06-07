#!/usr/bin/env python3

"""
.. module:: coverage
   :synopsis: Definitions of classes used to find, group and format missing topologies
    
.. moduleauthor:: Ursula Laa <ursula.laa@lpsc.in2p3.fr>    
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
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
from smodels.theory.exceptions import SModelSTheoryError as SModelSError


#Default definitions for the uncovered topology categories/groups:
##Element filters for each group:
##(it should be a function which takes an Element object as input
##and returns True if the element belongs to the group and False otherwise)
filtersDefault = {'outsideGrid (prompt)' : lambda el: ('prompt' in el.coveredBy) and not ('prompt' in el.testedBy),
                'outsideGrid (displaced)' : lambda el: ('displaced' in el.coveredBy) and not ('displaced' in el.testedBy),
                'missing (prompt)' : lambda el: not ('prompt' in el.coveredBy),
                'missing (displaced)' : lambda el: not ('displaced' in el.coveredBy),
                'missing (long cascade)' : lambda el: (not el.coveredBy) and el._getLength() > 3,
                'missing (all)' : lambda el: (not el.coveredBy),
                'outsideGrid (all)' : lambda el: (el.coveredBy and not el.testedBy)}

##Description for each group (optional and only used for printing):
##(if not defined, the group label will be used instead)
descriptionDefault = {'outsideGrid (prompt)' : 'topologies outside the mass grid with prompt decays',
                'outsideGrid (displaced)' : 'topologies outside the mass grid with displaced decays',
                'missing (prompt)' : 'missing topologies with prompt decays',
                'missing (displaced)' : 'missing topologies with displaced decays',
                'missing (long cascade)' : 'missing topologies with long cascade decays',
                'missing (all)' : 'missing topologies',
                'outsideGrid (all)' : 'topologies outside the mass grid'}

##Weight factors for each group:
##(it should be a function which takes an Element object as input
##and returns the reweighting factor to be applied to the element weight. It is relevant if only
##the fraction of the weight going into prompt or displaced decays is required)
factorsDefault = {}
for key in filtersDefault:
    if 'prompt' in key.lower():
        factorsDefault[key] = lambda el: reweightFactorFor(el,'prompt')
    elif 'displaced' in key.lower():
        factorsDefault[key] = lambda el: reweightFactorFor(el,'displaced')
    else:
        factorsDefault[key] = lambda el: 1.


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

    def __init__(self,topoList, sqrts=None,
                 groupFilters = filtersDefault,
                 groupFactors = factorsDefault,
                 groupdDescriptions = descriptionDefault,
                 smFinalStates=None,bsmFinalSates=None):

        #Sanity checks:
        if not isinstance(groupFilters,dict):
            raise SModelSError("groupFilters input should be a Dictionary and not %s" %type(groupFilters))
        if not isinstance(groupFactors,dict):
            raise SModelSError("groupFactors input should be a Dictionary and not %s" %type(groupFactors))
        if sorted(groupFilters.keys()) != sorted(groupFactors.keys()):
            raise SModelSError("Keys in groupFilters and groupFactors do not match")
        if any(not hasattr(gFilter,'__call__') for gFilter in groupFilters.values()):
            raise SModelSError("Group filters must be callable methods")


        if smFinalStates is None:
            #Define multiparticles for conveniently grouping SM final states
            taList = MultiParticle('ta', [ta,taC])
            lList = MultiParticle('l', [e,mu,eC,muC])
            WList = MultiParticle('W', [W,WC])
            tList = MultiParticle('t', [t,tC])
            jetList = MultiParticle('jet', [q,c,g,pion])
            nuList  = nu
            smFinalStates = [WList, lList, tList, taList, nuList, jetList]
        if bsmFinalSates is None:
            #Define inclusive BSM states to group/label the last BSM states:
            bsmFinalStates = [MET,HSCP,RHadronG,RHadronQ]
        
        if sqrts is None:
            sqrts = max([xsec.info.sqrts for xsec in topoList.getTotalWeight()])
        else:
            sqrts = sqrts
        
        self.groups = []
        #Create each uncovered group and get the topologies from topoList
        for gLabel,gFilter in groupFilters.items():
            #Initialize the uncovered topology list:
            uncoveredTopos = UncoveredList(label=gLabel,elementFilter=gFilter,
                                           reweightFactor = groupFactors[gLabel],
                                           smFinalStates=smFinalStates,
                                           bsmFinalStates=bsmFinalStates,
                                           sqrts=sqrts)
            if groupdDescriptions and gLabel in groupdDescriptions:
                uncoveredTopos.description = groupdDescriptions[gLabel]
            else:
                uncoveredTopos.description = gLabel
            #Fill the list with the elements in topoList:
            uncoveredTopos.getToposFrom(topoList)
            self.groups.append(uncoveredTopos)

    def getGroup(self,label):
        """
        Returns the group with the required label. If not found, returns None.

        :param label: String corresponding to the specific group label

        :return: UncoveredList object which matches the label
        """

        for group in self.groups:
            if group.label == label:
                return group

        return None

class UncoveredList(object):
    """
    Object to find and collect GeneralElement objects, plus printout functionality
    :ivar generalElements: missing elements, grouped by common general name using inclusive labels (e.g. jet)
    :ivar smFinalStates: Lists of MultiParticles used to group final states.
    :ivar sqrts: sqrts, for printout
    """
    def __init__(self, label, elementFilter, reweightFactor, smFinalStates, bsmFinalStates, sqrts):
        self.generalElements = []
        self.smFinalStates = smFinalStates
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
        missingXandEls = [[self.getMissingX(el)*self.reweightFactor(el),el] for el in elementList if self.reweightFactor(el)]
        #Sort according to largest missingX, smallest size and largest ID
        missingXandEls = sorted(missingXandEls, key = lambda pt: [pt[0],-pt[1]._getLength(),pt[1].elID], reverse=True)

        #Split lists of elements and missingX:
        missingXsecs = [pt[0] for pt in missingXandEls]
        elementList = [pt[1] for pt in missingXandEls]
        
        #Remove all elements which are related to each other in order to avoid double counting
        #(keep always the first appearance in the list, so we always keep the ones with largest missing xsec)
        elementListUnique = []
        missingXsecsUnique = []
        ancestors = set() #Keep track of all the ancestors of the elements in the unique list
        for i,element in enumerate(elementList):
            ancestorsIDs = set([el.elID for el in element.getAncestors() if el.elID != 0])
            #If the element has any common ancestor with any of the previous elements,
            #skip it to avoid double counting
            if ancestors.intersection(ancestorsIDs):
                continue
            elementListUnique.append(element)
            missingXsecsUnique.append(missingXsecs[i])
            ancestors = ancestors.union(ancestorsIDs)

        #Now that we only have unique elements with their effective missing cross-sections
        #we create General Elements out of them
        for i,el in enumerate(elementListUnique):
            missingX = missingXsecsUnique[i]
            if not missingX:
                continue
            self.addToGeneralElements(el,missingX)

        #Finally sort general elements by their missing cross-section:
        self.generalElements = sorted(self.generalElements[:], key = lambda genEl: genEl.missingX, reverse=True)

    def getMissingX(self,element):
        """
        Calculate total missing cross section of an element, by recursively checking if its
        mothers already appear in the list.
        :param element: Element object

        :returns: missing cross section without units (in fb)
        """

        ancestorList = element.getAncestors()
        alreadyChecked = [] # keep track of which elements have already been checked
        #If the element has no xsec for the required sqrts, return 0.
        if not element.weight.getXsecsFor(self.sqrts):
            return 0.
        #Get the total element weight:
        missingX =  element.weight.getXsecsFor(self.sqrts)[0].value
        overlapXsec = 0.*fb
        for ancestor in ancestorList: #check all ancestors (ancestorList is sorted by generation)
            #Skip entries which correspond to the element itself
            if ancestor is element:
                continue
            #Make sure we do not subtract the same mother twice
            if any(ancestor is el for el in alreadyChecked):
                continue
            alreadyChecked.append(ancestor)
            #Check if ancestor passes the group filter (if it has been covered/tested or not):
            if self.elementFilter(ancestor):
                continue
            #Subtract the weight of the ancestor and skip checking for all the older family tree of the ancestor
            #(avoid subtracting the same weight twice from mother and grandmother for instance)
            #(since the ancestorList is sorted by generation, the mother always
            #appears before the grandmother in the list)
            alreadyChecked += ancestor.getAncestors()
            ancestorXsec = ancestor.weight.getXsecsFor(self.sqrts)[0]
            if not ancestorXsec:
                continue
            overlapXsec += ancestorXsec.value

        missingX -= overlapXsec

        return missingX.asNumber(fb)
        
        
    def getTotalXSec(self, sqrts=None):
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

        newGenEl = GeneralElement(el,missingX,self.smFinalStates,self.bsmFinalStates)

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
    The even particles are replaced/grouped by the particles defined in smFinalStates.
    """

    def __init__(self,el,missingX,smFinalStates,bsmFinalStates):

        self.branches = [Branch() for _ in el.branches]
        self.missingX = missingX
        self.finalBSMstates = []

        for ib,branch in enumerate(el.branches):
            newParticles = []
            for vertex in branch.evenParticles:
                newVertex = vertex[:]
                for ip,particle in enumerate(vertex):
                    for particleList in smFinalStates:
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

