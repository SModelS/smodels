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
from smodels.theory.particle import MultiParticle, ParticleList, Particle
from smodels.share.models.SMparticles import e,mu,ta,taC,eC,muC,W,WC,t,tC,q,c,g,pion,nu
from smodels.theory.auxiliaryFunctions import index_bisect
from smodels.theory.exceptions import SModelSTheoryError as SModelSError


#Default definitions for the uncovered topology categories/groups:
##Element filters for each group:
##(it should be a function which takes an Element object as input
##and returns True if the element belongs to the group and False otherwise)
filtersDefault = {'missing (prompt)' : lambda el: not ('prompt' in el.coveredBy),
                'missing (displaced)' : lambda el: not ('displaced' in el.coveredBy),
#                 'missing (long cascade)' : lambda el: (not el.coveredBy) and el._getLength() > 3,
                'missing (all)' : lambda el: (not el.coveredBy),
                'outsideGrid (all)' : lambda el: (el.coveredBy and not el.testedBy)}

##Description for each group (optional and only used for printing):
##(if not defined, the group label will be used instead)
descriptionDefault = {'missing (prompt)' : 'missing topologies with prompt decays',
                'missing (displaced)' : 'missing topologies with displaced decays',
#                 'missing (long cascade)' : 'missing topologies with long cascade decays',
                'missing (all)' : 'missing topologies',
                'outsideGrid (all)' : 'topologies outside the grid'}

##Default final BSM states for grouping topologies:
#Used to construct BSM final states:
MET = Particle(label='MET', Z2parity = -1, eCharge = 0, colordim = 1)
HSCPp = Particle(label='HSCP+', Z2parity = -1, eCharge = +1, colordim = 1)
HSCPm = Particle(label='HSCP-', Z2parity = -1, eCharge = -1, colordim = 1)
HSCP = MultiParticle(label='HSCP', particles = [HSCPp,HSCPm])

RHadronG = Particle(label='RHadronG', Z2parity = -1, eCharge = 0, colordim = 8)
RHadronU = Particle(label='RHadronU', Z2parity = -1, eCharge = 2./3., colordim = 3)
RHadronD = Particle(label='RHadronD', Z2parity = -1, eCharge = -1./3., colordim = 3)
RHadronQ = MultiParticle(label='RHadronQ', particles = [RHadronU,RHadronU.chargeConjugate(),
                                                        RHadronD,RHadronD.chargeConjugate()])
bsmDefault = [MET,HSCP,RHadronG,RHadronQ]

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
        #If not specified assumed all fractions
        #(note that we can not include any long-lived fraction since this is already included in
        #the topologies where the meta-stable particle appears as a final state,
        #so the total is = (fraction of all decays being prompt)
        #+ (fraction of at least one displaced decay and no long-lived decays)
        factorsDefault[key] = lambda el: reweightFactorFor(el,'prompt') + reweightFactorFor(el,'displaced')


class Uncovered(object):
    """
    Wrapper object for defining and holding a list of coverage groups  (UncoveredGroup objects).

    The class builds a series of UncoveredGroup objects and stores them.
    """

    def __init__(self,topoList, sqrts=None, sigmacut=0*fb,
                 groupFilters = filtersDefault,
                 groupFactors = factorsDefault,
                 groupdDescriptions = descriptionDefault,
                 smFinalStates=None,bsmFinalSates=None):
        """
        Inititalize the object.

        :param topoList: TopologyList object used to select elements from.
        :param sqrts: Value (with units) for the center of mass energy used to compute the missing cross sections.
                     If not specified the largest value available will be used.
        :param sigmacut: Minimum cross-section/weight value (after applying the reweight factor)
                       for an element to be included. The value should in fb (unitless)
        :param groupFilters: Dictionary containing the groups' labels and the method for selecting
                            elements.
        :param groupFactors: Dictionary containing the groups' labels and the method for reweighting
                            cross sections.
        :param groupdDescriptions: Dictionary containing the groups' labels and strings describing the group
                                  (used for printout)
        :param smFinalStates: List of (inclusive) Particle or MultiParticle objects used for grouping Z2-even
                             particles when creating GeneralElements.
        :param bsmFinalSates: List of (inclusive) Particle or MultiParticle objects used for grouping Z2-odd
                             particles when creating GeneralElements.
        """

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
            bsmFinalStates = bsmDefault

        if sqrts is None:
            sqrts = max([xsec.info.sqrts for xsec in topoList.getTotalWeight()])
        else:
            sqrts = sqrts

        #Store the relevant element cross-sections to improve performance:
        for el in topoList.getElements():
            xsec = el.weight.getXsecsFor(sqrts)
            if xsec:
                el._totalXsec = xsec[0].value.asNumber(fb)
            else:
                el._totalXsec = 0.

        self.groups = []
        #Create each uncovered group and get the topologies from topoList
        for gLabel,gFilter in groupFilters.items():
            #Initialize the uncovered topology list:
            uncoveredTopos = UncoveredGroup(label=gLabel,elementFilter=gFilter,
                                           reweightFactor = groupFactors[gLabel],
                                           smFinalStates=smFinalStates,
                                           bsmFinalStates=bsmFinalStates,
                                           sqrts=sqrts, sigmacut=sigmacut.asNumber(fb))
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

        :return: UncoveredGroup object which matches the label
        """

        for group in self.groups:
            if group.label == label:
                return group

        return None

class UncoveredGroup(object):
    """
    Holds information about a single coverage group: criteria for selecting and grouping elements,
    function for reweighting cross sections, etc.
    """

    def __init__(self, label, elementFilter, reweightFactor,
                 smFinalStates, bsmFinalStates, sqrts, sigmacut=0.):
        """
        :param label: Group label
        :param elementFilter: Function which takes an element as argument and returns True (False) if
                             the element should (not) be selected.
        :param reweightFactor: Function which takes an element as argument and returns the reweighting
                              factor to be applied to the element weight.
        :param smFinalStates: List of Particle/MultiParticle objects used to group Z2-even particles appearing
                            in the final state
        :param bsmFinalStates: List of Particle/MultiParticle objects used to group Z2-odd particles appearing
                            in the final state
        :param sqrts: Value (with units) for the center of mass energy used to compute the missing cross sections.
                     If not specified the largest value available will be used.
        :param sigmacut: Minimum cross-section/weight value (after applying the reweight factor)
                       for an element to be included. The value should in fb (unitless)
        """

        self.generalElements = []
        self.smFinalStates = smFinalStates
        self.bsmFinalStates = bsmFinalStates
        self.sqrts = sqrts
        self.sigmacut = sigmacut
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
        #Only keep the ones elements sigmacut:
        if self.sigmacut:
            missingXandEls = [x for x in missingXandEls[:] if x[0] > self.sigmacut]
        else:
            missingXandEls = [x for x in missingXandEls[:] if x[0] > 0.]
        #Sort according to largest missingX, smallest size and largest ID
        missingXandEls = sorted(missingXandEls, key = lambda pt: [pt[0],-pt[1]._getLength(),-pt[1].elID], reverse=True)

        #Split lists of elements and missingX:
        missingXsecs = [pt[0] for pt in missingXandEls]
        elementList = [pt[1] for pt in missingXandEls]

        #Remove all elements which are related to each other in order to avoid double counting
        #(keep always the first appearance in the list, so we always keep the ones with largest missing xsec)
        elementListUnique = []
        missingXsecsUnique = []
        ancestors = set() #Keep track of all the ancestor ids of the elements in the unique list
        for i,element in enumerate(elementList):
            #Get ancestor object ids
            #(safer than element.elID since some elements can share elID = 0 if they
            #were never directly inserted in the TopologyList)
            ancestorsIDs = set([id(element)] + [id(el) for el in element.getAncestors()])
            #If the element has any common ancestor with any of the previous elements
            #or if it is an ancestor of any of the previous elements
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
        #Get the (pre-loaded) total element weight in fb:
        missingX =  element._totalXsec
        if not missingX:
            return 0.
        overlapXsec = 0.
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
            if not hasattr(ancestor,'_totalXsec'):
                xsec = ancestor.weight.getXsecsFor(self.sqrts)
                ancestor._totalXsec = xsec[0].value.asNumber(fb)
            overlapXsec += ancestor._totalXsec

        return missingX-overlapXsec


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
