"""
.. module:: element
   :synopsis: Module holding the Element class and its methods.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.theory.auxiliaryFunctions import elementsInStr
from smodels.theory.branch import Branch, InclusiveBranch
from smodels.theory import crossSection
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger
import itertools
from smodels.theory.particle import Particle

class Element(object):
    """
    An instance of this class represents an element.
    This class possesses a pair of branches and the element weight
    (cross-section * BR).
    """

    def __init__(self, info=None, finalState=None, intermediateState=None, model=None):
        """
        Initializes the element. If info is defined, tries to generate
        the element using it.

        :parameter info: string describing the element in bracket notation
                         (e.g. [[[e+],[jet]],[[e-],[jet]]])

        :parameter finalState: list containing the final state labels for each branch
                         (e.g. ['MET', 'HSCP'] or ['MET','MET'])
        :parameter intermediateState: nested list containing intermediate state labels
                                     for each branch  (e.g. [['gluino'], ['gluino']])
        :parameter model: The model (Model object) to be used when converting particle labels to
                          particle objects (only used if info, finalState or intermediateState != None).
        """

        self.branches = [Branch(), Branch()]
        self.weight = crossSection.XSectionList() # gives the weight for all decays promptly
        self.decayLabels = []
        self.motherElements = [self] #The motheElements includes self to keep track of merged elements
        self.elID = 0
        self.coveredBy = set()
        self.testedBy = set()

        if info:
            # Create element from particle string
            if type(info) == type(str()):
                elements = elementsInStr(info)
                if not elements or len(elements) > 1:
                    nel = 0
                    if elements:
                        nel = len(elements)
                    logger.error("Malformed input string. Number of elements "
                                  "is %d (expected 1) in: ``%s''", nel, info)
                    return None
                else:
                    el = elements[0]
                    branches = elementsInStr(el[1:-1])
                    if not branches or len(branches) != 2:
                        logger.error("Malformed input string. Number of "
                                      "branches is %d (expected 2) in: ``%s''",
                                      len(branches), info)
                        return None
                    self.branches = []
                    if intermediateState:
                        if len(intermediateState) != len(branches):
                            raise SModelSError("Number of intermediate states (%i) does not match number of branches (%i)"
                                               %(len(intermediateState),len(branches)))
                    else:
                        intermediateState = [None]*len(branches)

                    if finalState:
                        if len(finalState) != len(branches):
                            raise SModelSError("Number of final states (%i) does not match number of branches (%i)"
                                               %(len(finalState),len(branches)))
                    else:
                        finalState = [None]*len(branches)
                    for ibr,branch in enumerate(branches):
                        if branch == '[*]':
                            self.branches.append(InclusiveBranch(finalState[ibr],
                                                  intermediateState=intermediateState[ibr],model=model))
                        else:
                            self.branches.append(Branch(branch,finalState[ibr],
                                                    intermediateState[ibr],model=model))

            # Create element from branch pair
            elif isinstance(info,list) and all(isinstance(x,(Branch,InclusiveBranch)) for x in info):
                self.branches = [br.copy() for br in info]
            else:
                raise SModelSError("Can not create element from input type %s" %type(info))

        self.setEinfo()


    def __cmp__(self,other):
        """
        Compares the element with other for any branch ordering.
        The comparison is made based on branches.
        OBS: The elements and the branches must be sorted!
        :param other:  element to be compared (Element object)
        :return: -1 if self < other, 0 if self == other, +1, if self > other.
        """

        if not isinstance(other,Element):
            return -1

        #Compare branches:
        for branchA in itertools.permutations(self.branches):
            branchA = list(branchA)
            if branchA == other.branches:
                return 0

        comp = (self.branches > other.branches)
        if comp:
            return 1
        else:
            return -1

    def __eq__(self,other):
        return self.__cmp__(other)==0

    def __lt__(self,other):
        return self.__cmp__(other)<0

    def __hash__(self):
        return object.__hash__(self)

    def __getattr__(self, attr):
        """
        If the attribute has not been defined for the element
        try to fetch it from its branches.
        :param attr: Attribute name

        :return: Attribute value
        """

        #If calling another special method, return default (required for pickling)
        if attr.startswith('__') and attr.endswith('__'):
            return object.__getattr__(attr)

        try:
            val = [getattr(br,attr) for br in self.branches]
            return val
        except AttributeError:
            raise AttributeError("Neither element nor branch has attribute ``%s''" %attr)

    def __str__(self):
        """
        Create the element bracket notation string, e.g. [[[jet]],[[jet]].

        :returns: string representation of the element (in bracket notation)
        """

        elStr = "["+",".join([str(br) for br in self.branches])+"]"
        elStr = elStr.replace(" ", "").replace("'", "")
        return elStr

    def __repr__(self):

        return self.__str__()

    def __add__(self,other):
        """
        Adds two elements. Should only be used if the elements
        have the same topologies. The element weights are added and their
        odd and even particles are combined.
        """

        if not isinstance(other,Element):
            raise TypeError("Can not add an Element object to %s" %type(other))
        elif self.getEinfo() != other.getEinfo():
            raise SModelSError("Can not add elements with distinct topologies")

        newEl = self.__class__()
        newEl.motherElements = self.motherElements[:] + other.motherElements[:]
        newEl.weight = self.weight + other.weight
        newEl.branches = []
        for ibr,branch in enumerate(self.branches):
            newEl.branches.append(branch + other.branches[ibr])

        return newEl

    def __radd__(self,other):
        """
        Adds two elements. Only elements with the same
        topology can be combined.
        """

        return self.__add__(other)

    def __iadd__(self,other):
        """
        Combine two elements. Should only be used if the elements
        have the same topologies. The element weights are added and their
        odd and even particles are combined.
        """

        if not isinstance(other,Element):
            raise TypeError("Can not add an Element object to %s" %type(other))
        elif self.getEinfo() != other.getEinfo():
            raise SModelSError("Can not add elements with distinct topologies")

        self.motherElements += other.motherElements[:]
        self.weight += other.weight
        for ibr,_ in enumerate(self.branches):
            self.branches[ibr] += other.branches[ibr]

        return self

    def getAverage(self,attr):
        """
        Get the average value for a given attribute appearing in
        the odd particles of the element branches.
        """

        try:
            vals = [br.getAverage(attr) for br in self.branches]
        except (AttributeError,ZeroDivisionError):
            raise SModelSError("Could not compute average for %s" %attr)

        return vals

    def toStr(self):
        """
        Returns a string with the element represented in bracket notation,
        including the final states, e.g. [[[jet]],[[jet]] (MET,MET)
        """

        elStr = str(self)+' '+str(tuple(self.getFinalStates())).replace("'","")

        return elStr

    def sortBranches(self):
        """
        Sort branches. The smallest branch is the first one.
        See the Branch object for definition of branch size and comparison
        """

        #Now sort branches
        self.branches = sorted(self.branches)

    def copy(self):
        """
        Create a copy of self.
        Faster than deepcopy.

        :returns: copy of element (Element object)
        """

        #Allows for derived classes (like inclusive classes)
        newel = self.__class__()
        newel.branches = []
        for branch in self.branches:
            newel.branches.append(branch.copy())
        newel.weight = self.weight.copy()
        newel.motherElements = self.motherElements[:]
        newel.elID = self.elID
        return newel

    def switchBranches(self):
        """
        Switch branches, if the element contains a pair of them.

        :returns: element with switched branches (Element object)
        """

        newEl = self.copy()
        if len(self.branches) == 2:
            newEl.branches = [newEl.branches[1], newEl.branches[0]]
        return newEl

    def getFinalStates(self):
        """
        Get the array of final state (last BSM particle) particle objects in the element.

        :returns: list of Particle objects
        """

        fsparticles = [branch.oddParticles[-1] for branch in self.branches]

        return fsparticles

    def _getAncestorsDict(self,igen=0):
        """
        Returns a dictionary with all the ancestors
        of the element. The dictionary keys are integers
        labeling the generation (number of generations away from self)
        and the values are a list of Element objects (ancestors) for that generation.
        igen is used as the counter for the initial generation.
        The output is also stored in self._ancestorsDict for future use.

        :param igen: Auxiliary integer indicating to which generation self belongs.

        :return: Dictionary with generation index as key and ancestors as values
                 (e.g. {igen+1 : [mother1, mother2], igen+2 : [grandmother1,..],...})
        """

        ancestorsDict = {igen+1 : []}
        for mother in self.motherElements:
            if mother is self:
                continue
            ancestorsDict[igen+1].append(mother)
            for jgen,elList in mother._getAncestorsDict(igen+1).items():
                if not jgen in ancestorsDict:
                    ancestorsDict[jgen] = []
                ancestorsDict[jgen] += elList

        #Store the result
        self._ancestorsDict = dict([[key,val] for key,val in ancestorsDict.items()])

        return self._ancestorsDict

    def getAncestors(self):
        """
        Get a list of all the ancestors of the element.
        The list is ordered so the mothers appear first, then the grandmother,
        then the grandgrandmothers,...

        :return: A list of Element objects containing all the ancestors sorted by generation.
        """

        #Check if the ancestors have already been obtained (performance gain)
        if not hasattr(self,'_ancestorsDict'):
            self._getAncestorsDict()

        orderedAncestors = []
        for jgen in sorted(self._ancestorsDict.keys()):
            orderedAncestors += self._ancestorsDict[jgen]

        return orderedAncestors

    def isRelatedTo(self,other):
        """
        Checks if the element has any common ancestors with other or one
        is an ancestor of the other.
        Returns True if self and other have at least one ancestor in common
        or are the same element, otherwise returns False.

        :return: True/False
        """

        ancestorsA = set([id(self)] + [id(el) for el in self.getAncestors()])
        ancestorsB = set([id(other)] + [id(el) for el in other.getAncestors()])

        if ancestorsA.intersection(ancestorsB):
            return True
        else:
            return False

    def getEinfo(self):
        """
        Get element topology info from branch topology info.

        :returns: dictionary containing vertices and number of final states information
        """

        vertnumb = []
        vertparts = []
        for branch in self.branches:
            if branch.vertnumb is None:
                branch.setInfo()
            bInfo = branch.getInfo()
            vertparts.append(bInfo['vertparts'])
            vertnumb.append(bInfo['vertnumb'])

        return {"vertnumb" : vertnumb, "vertparts" : vertparts}

    def setEinfo(self):
        """
        Set topology info for each branch.
        """

        for branch in self.branches:
            branch.setInfo()

    def setTestedBy(self,resultType):
        """
        Tag the element, all its daughter and all its mothers
        as tested by the type of result specified.
        It also recursively tags all granddaughters, grandmothers,...

        :param resultType: String describing the type of result (e.g. 'prompt', 'displaced')
        """

        self.testedBy.add(resultType)
        for ancestor in self.getAncestors():
            ancestor.testedBy.add(resultType)

    def setCoveredBy(self,resultType):
        """
        Tag the element, all its daughter and all its mothers
        as covered by the type of result specified.
        It also recursively tags all granddaughters, grandmothers,...

        :param resultType: String describing the type of result (e.g. 'prompt', 'displaced')
        """

        self.coveredBy.add(resultType)
        for mother in self.getAncestors():
            mother.coveredBy.add(resultType)

    def _getLength(self):
        """
        Get the maximum of the two branch lengths.

        :returns: maximum length of the element branches (int)
        """
        return max(self.branches[0].getLength(), self.branches[1].getLength())

    def checkConsistency(self):
        """
        Check if the particles defined in the element are consistent
        with the element info.

        :returns: True if the element is consistent. Print error message
                  and exits otherwise.
        """
        info = self.getEinfo()
        for ib, branch in enumerate(self.branches):
            for iv, vertex in enumerate(branch.evenParticles):
                if len(vertex) != info['vertparts'][ib][iv]:
                    logger.error("Wrong syntax")
                    raise SModelSError()
        return True


    def compressElement(self, doCompress, doInvisible, minmassgap):
        """
        Keep compressing the original element and the derived ones till they
        can be compressed no more.

        :parameter doCompress: if True, perform mass compression
        :parameter doInvisible: if True, perform invisible compression
        :parameter minmassgap: value (in GeV) of the maximum
                               mass difference for compression
                               (if mass difference < minmassgap, perform mass compression)
        :returns: list with the compressed elements (Element objects)
        """

        if not doCompress and not doInvisible:
            return []

        added = True
        newElements = [self]
        # Keep compressing the new topologies generated so far until no new
        # compressions can happen:
        while added:
            added = False
            # Check for mass compressed topologies
            if doCompress:
                for element in newElements:
                    newel = element.massCompress(minmassgap)
                    # Avoids double counting
                    #(elements sharing the same parent are removed during clustering)
                    if newel and not any(newel == el for el in newElements[:]):
                        newElements.append(newel)
                        added = True

            # Check for invisible compressed topologies (look for effective
            # LSP, such as LSP + neutrino = LSP')
            if doInvisible:
                for element in newElements:
                    newel = element.invisibleCompress()
                    # Avoids double counting
                    #(elements sharing the same parent are removed during clustering)
                    if newel and not any(newel == el for el in newElements[:]):
                        newElements.append(newel)
                        added = True

        newElements.pop(0)  # Remove original element
        return newElements


    def removeVertex(self,ibr,iv):
        """
        Remove vertex iv located in branch ibr.
        The "vertex-mother" in BSMparticles and (SM) particles in the vertex
        are removed from the branch. The vertex index corresponds
        to the BSM decay (iv = 0 will remove the first BSM particle,...)

        :parameter ibr: Index of branch (int)
        :parameter iv: Index of vertex in branch ibr (int)

        """

        self.branches[ibr].removeVertex(iv)

    def massCompress(self, minmassgap):
        """
        Perform mass compression.

        :parameter minmassgap: value (in GeV) of the maximum
                               mass difference for compression
                               (if mass difference < minmassgap -> perform mass compression)
        :returns: compressed copy of the element, if two masses in this
                  element are degenerate; None, if compression is not possible;
        """

        newelement = self.copy()
        newelement.motherElements = [self]

        #Loop over branches and look for small mass differences
        for ibr,branch in enumerate(newelement.branches):
            #Get mass differences

            removeVertices = []
            for i,mom in enumerate(branch.oddParticles[:-1]):
                massDiff = mom.mass - branch.oddParticles[i+1].mass
                #Get vertices which have deltaM < minmassgap and the mother is prompt:
                if massDiff < minmassgap and mom.isPrompt():
                    removeVertices.append(i)
            #Remove vertices till they exist:
            while removeVertices:
                newelement.removeVertex(ibr,removeVertices[0])
                branch = newelement.branches[ibr]
                removeVertices = []
                for i,mom in enumerate(branch.oddParticles[:-1]):
                    massDiff = mom.mass - branch.oddParticles[i+1].mass
                    #Get vertices which have deltaM < minmassgap and the mother is prompt:
                    if massDiff < minmassgap and mom.isPrompt():
                        removeVertices.append(i)

        for ibr,branch in enumerate(newelement.branches):
            if branch.vertnumb != self.branches[ibr].vertnumb:
                newelement.sortBranches()
                return newelement

        #New element was not compressed, return None
        return None


    def invisibleCompress(self):
        """
        Perform invisible compression.

        :returns: compressed copy of the element, if element ends with invisible
                  particles; None, if compression is not possible
        """

        newelement = self.copy()
        newelement.motherElements = [self]

        # Loop over branches
        for branch in newelement.branches:
            if not branch.evenParticles:
                continue
            #Check if the last decay should be removed:
            neutralSM = all(ptc.isMET() for ptc in branch.evenParticles[-1])
            neutralDecay = neutralSM and branch.oddParticles[-1].isMET()
            #Check if the mother can be considered MET:
            neutralBSM = (branch.oddParticles[-2].isMET()
                            or branch.oddParticles[-2].isPrompt())
            if neutralBSM and neutralDecay:
                removeLastVertex = True
            else:
                removeLastVertex = False

            while len(branch.oddParticles) > 1 and removeLastVertex:
                bsmMom = branch.oddParticles[-2]
                effectiveDaughter = Particle(label='inv', mass = bsmMom.mass,
                                             eCharge = 0, colordim = 1,
                                             totalwidth = branch.oddParticles[-1].totalwidth,
                                             Z2parity = bsmMom.Z2parity, pdg = bsmMom.pdg)
                branch.removeVertex(len(branch.oddParticles)-2)
                #For invisible compression, keep an effective mother which corresponds to the invisible
                #daughter, but with the mass of the parent.
                branch.oddParticles[-1] = effectiveDaughter
                #Re-check if the last decay should be removed:
                if not branch.evenParticles:
                    continue
                neutralSM = all(ptc.isMET() for ptc in branch.evenParticles[-1])
                neutralDecay = neutralSM and branch.oddParticles[-1].isMET()
                neutralBSM = branch.oddParticles[-2].isMET()
                if neutralBSM and neutralDecay:
                    removeLastVertex = True
                else:
                    removeLastVertex = False

        for ibr,branch in enumerate(newelement.branches):
            if branch.vertnumb != self.branches[ibr].vertnumb:
                newelement.sortBranches()
                return newelement

        #New element was not compressed, return None
        return None


    def hasTopInList(self, elementList):
        """
        Check if the element topology matches any of the topologies in the
        element list.

        :parameter elementList: list of elements (Element objects)
        :returns: True, if element topology has a match in the list, False otherwise.
        """
        if type(elementList) != type([]) or len(elementList) == 0:
            return False
        for element in elementList:
            if type(element) != type(self):
                continue
            info1 = self.getEinfo()
            info2 = element.getEinfo()
            info2B = element.switchBranches().getEinfo()
            if info1 == info2 or info1 == info2B:
                return True
        return False
