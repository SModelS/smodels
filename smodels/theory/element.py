"""
.. module:: element
   :synopsis: Module holding the Element class and its methods.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.theory.graphTools import stringToTree, getCanonName, treeToString, compareNodes, drawTree, getTreeRoot, sortTree
from smodels.theory import crossSection
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from networkx import DiGraph


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
                         (e.g. [[[e+],[jet]],[[e-],[jet]]]) or process notation
                         (e.g (PV > gluino(1),gluino(2)), (gluino(1) > u,u~,N1), (gluino(2) > d,d~,N1))
                         OR a tree (DiGraph object)

        :parameter finalState: list containing the final state labels for each branch
                         (e.g. ['MET', 'HSCP'] or ['MET','MET']). Only needed for the bracket notation
                         input.
        :parameter intermediateState: nested list containing intermediate state labels
                                      for each branch  (e.g. [['gluino'], ['gluino']]). Only
                                      needed for the bracket notation input.

        :parameter model: The model (Model object) to be used when converting particle labels to
                          particle objects (only used if info, finalState or intermediateState != None).
        """
        self.tree = DiGraph()
        self.weight = crossSection.XSectionList()  # gives the weight for all decays promptly
        self.decayLabels = []
        self.motherElements = [self]  # The motheElements includes self to keep track of merged elements
        self.elID = 0
        self.coveredBy = set()
        self.testedBy = set()

        if info:
            if isinstance(info, str):
                try:
                    self.tree = stringToTree(info, finalState=finalState,
                                             intermediateState=intermediateState,
                                             model=model)
                except (SModelSError, TypeError):
                    raise SModelSError("Can not create element from input %s" % info)
            elif isinstance(info, DiGraph):
                self.tree = info.copy()  # Makes a shallow copy of the original tree
            else:
                raise SModelSError("Can not create element from input type %s" % type(info))

        self.setCanonName()
        self.sort()

    def setCanonName(self):
        """
        Compute and store the canonical name for the tree
        topology. The canonical name can be used to compare and sort topologies.
        The name is stored in self.tree.canonName
        """

        canonName = getCanonName(self.tree)
        self.tree.graph['canonName'] = canonName

    def sort(self):
        """
        Sort the tree according to the canonName and then particle. For each node,
        all the daughters are sorted according to their canonName.
        """

        # If canon names have not been defined, compute them:
        if 'canonName' not in self.tree.graph:
            self.setCanonName()

        self.tree = sortTree(self.tree)

    def compareTo(self, other, sortOther=False):
        """
        Compares the element with other.
        Uses the topology name (Tree canonincal name) to identify isomorphic topologies (trees).
        If trees are not isomorphic, compare topologies using the topology name.
        Else  check if nodes (particles) are equal. For nodes with the same canonName
        and the same parent, compare against all permutations.
        If the canonical names match, but the particle differs, returns the particle
        comparison.

        :param other:  element to be compared (Element object)
        :param sortOther: if True and the elements match,
                          also returns a copy of other, but with its tree sorted
                          according to the way it matched self.
                          If True, but elements do not match, return None.

        :return: -1 if self < other, 0 if self == other, +1, if self > other.
        """

        if not isinstance(other, Element):
            return -1

        # make sure the topology names have been computed:
        if 'canonName' not in self.tree.graph:
            self.setCanonName()
        if 'canonName' not in other.tree.graph:
            other.setCanonName()

        # Recursively compare the nodes:
        cmp, newTree = compareNodes(self.tree, other.tree,
                                    getTreeRoot(self.tree), getTreeRoot(other.tree))

        if sortOther:
            if cmp == 0:  # Elements matched, return copy of other with tree sorted
                otherNew = other.copy()
                otherNew.tree = newTree
            else:
                otherNew = None
            return cmp, otherNew
        else:
            return cmp

    def __eq__(self, other):
        return self.compareTo(other) == 0

    def __lt__(self, other):
        return self.compareTo(other) < 0

    def __gt__(self, other):
        return self.compareTo(other) > 0

    def __hash__(self):
        return object.__hash__(self)

    def __getattr__(self, attr):
        """
        If the attribute has not been defined for the element
        try to fetch it from its branches.
        :param attr: Attribute name

        :return: Attribute value
        """

        # If calling another special method, return default (required for pickling)
        if attr.startswith('__') and attr.endswith('__'):
            return object.__getattr__(attr)

        try:
            val = [getattr(node, attr) if str(node) != 'PV' else None
                   for node in self.tree.nodes()]
            return val
        except AttributeError:
            raise AttributeError("Neither element nor particles have attribute ``%s''" % attr)

    def __str__(self):
        """
        Create the element to a string.

        :returns: string representation of the element
        """

        return treeToString(self.tree)

    def __repr__(self):

        return self.__str__()

    # def __add__(self, other):
    #     """
    #     Adds two elements. Should only be used if the elements
    #     have the same topologies. The element weights are added and their
    #     odd and even particles are combined.
    #     """
    #
    #     if not isinstance(other, Element):
    #         raise TypeError("Can not add an Element object to %s" % type(other))
    #     elif self.getCanonName() != other.getCanonName():
    #         raise SModelSError("Can not add elements with distinct topologies")
    #
    #     newEl = self.__class__()
    #     newEl.motherElements = self.motherElements[:] + other.motherElements[:]
    #     newEl.weight = self.weight + other.weight
    #     newEl.branches = []
    #     for ibr, branch in enumerate(self.branches):
    #         newEl.branches.append(branch + other.branches[ibr])
    #
    #     return newEl
    #
    # def __radd__(self, other):
    #     """
    #     Adds two elements. Only elements with the same
    #     topology can be combined.
    #     """
    #
    #     return self.__add__(other)
    #
    # def __iadd__(self, other):
    #     """
    #     Combine two elements. Should only be used if the elements
    #     have the same topologies. The element weights are added and their
    #     odd and even particles are combined.
    #     """
    #
    #     if not isinstance(other, Element):
    #         raise TypeError("Can not add an Element object to %s" % type(other))
    #     elif self.getCanonName() != other.getCanonName():
    #         raise SModelSError("Can not add elements with distinct topologies")
    #
    #     self.motherElements += other.motherElements[:]
    #     self.weight += other.weight
    #     for ibr, _ in enumerate(self.branches):
    #         self.branches[ibr] += other.branches[ibr]
    #
    #     return self
    #
    # def getAverage(self, attr):
    #     """
    #     Get the average value for a given attribute appearing in
    #     the odd particles of the element branches.
    #     """
    #
    #     try:
    #         vals = [br.getAverage(attr) for br in self.branches]
    #     except (AttributeError, ZeroDivisionError):
    #         raise SModelSError("Could not compute average for %s" % attr)
    #
    #     return vals

    def drawTree(self, oddColor='lightcoral', evenColor='skyblue',
                 pvColor='darkgray', genericColor='violet',
                 nodeScale=4, labelAttr=None):
        """
        Draws Tree using matplotlib.

        :param tree: tree to be drawn
        :param oddColor: color for Z2-odd particles
        :param evenColor: color for Z2-even particles
        :param genericColor: color for particles without a defined Z2-parity
        :param pvColor: color for primary vertex
        :param nodeScale: scale size for nodes
        :param labelAttr: attribute to be used as label. If None, will use the string representation
                          of the node object.

        """

        return drawTree(self.tree, oddColor=oddColor,
                        evenColor=evenColor,
                        pvColor=pvColor, genericColor=genericColor,
                        nodeScale=nodeScale, labelAttr=labelAttr)

    def copy(self):
        """
        Create a copy of self.
        Faster than deepcopy.

        :returns: copy of element (Element object)
        """

        # Allows for derived classes (like inclusive classes)
        newel = self.__class__()
        newel.tree = self.tree.copy()
        newel.weight = self.weight.copy()
        newel.motherElements = self.motherElements[:]
        newel.elID = self.elID
        return newel

    def getFinalStates(self):
        """
        Get the list of particles which have not decayed (within the element)

        :returns: list of ParticleNode objects
        """

        finalStates = []
        for node in self.tree.nodes():
            if not list(self.tree.successors(node)):
                finalStates.append(node)

        return finalStates

    def _getAncestorsDict(self, igen=0):
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

        ancestorsDict = {igen+1: []}
        for mother in self.motherElements:
            if mother is self:
                continue
            ancestorsDict[igen+1].append(mother)
            for jgen, elList in mother._getAncestorsDict(igen+1).items():
                if jgen not in ancestorsDict:
                    ancestorsDict[jgen] = []
                ancestorsDict[jgen] += elList

        # Store the result
        self._ancestorsDict = dict([[key, val] for key, val in ancestorsDict.items()])

        return self._ancestorsDict

    def getAncestors(self):
        """
        Get a list of all the ancestors of the element.
        The list is ordered so the mothers appear first, then the grandmother,
        then the grandgrandmothers,...

        :return: A list of Element objects containing all the ancestors sorted by generation.
        """

        # Check if the ancestors have already been obtained (performance gain)
        if not hasattr(self, '_ancestorsDict'):
            self._getAncestorsDict()

        orderedAncestors = []
        for jgen in sorted(self._ancestorsDict.keys()):
            orderedAncestors += self._ancestorsDict[jgen]

        return orderedAncestors

    def isRelatedTo(self, other):
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

    def getCanonName(self):
        """
        Get element topology info from branch topology info.

        :returns: dictionary containing vertices and number of final states information
        """
        if 'canonName' not in self.tree.graph:
            self.setCanonName()
        return self.tree.graph['canonName']

    def setTestedBy(self, resultType):
        """
        Tag the element, all its daughter and all its mothers
        as tested by the type of result specified.
        It also recursively tags all granddaughters, grandmothers,...

        :param resultType: String describing the type of result (e.g. 'prompt', 'displaced')
        """

        self.testedBy.add(resultType)
        for ancestor in self.getAncestors():
            ancestor.testedBy.add(resultType)

    def setCoveredBy(self, resultType):
        """
        Tag the element, all its daughter and all its mothers
        as covered by the type of result specified.
        It also recursively tags all granddaughters, grandmothers,...

        :param resultType: String describing the type of result (e.g. 'prompt', 'displaced')
        """

        self.coveredBy.add(resultType)
        for mother in self.getAncestors():
            mother.coveredBy.add(resultType)

    # def compressElement(self, doCompress, doInvisible, minmassgap):
    #     """
    #     Keep compressing the original element and the derived ones till they
    #     can be compressed no more.
    #
    #     :parameter doCompress: if True, perform mass compression
    #     :parameter doInvisible: if True, perform invisible compression
    #     :parameter minmassgap: value (in GeV) of the maximum
    #                            mass difference for compression
    #                            (if mass difference < minmassgap, perform mass compression)
    #     :returns: list with the compressed elements (Element objects)
    #     """
    #
    #     if not doCompress and not doInvisible:
    #         return []
    #
    #     added = True
    #     newElements = [self]
    #     # Keep compressing the new topologies generated so far until no new
    #     # compressions can happen:
    #     while added:
    #         added = False
    #         # Check for mass compressed topologies
    #         if doCompress:
    #             for element in newElements:
    #                 newel = element.massCompress(minmassgap)
    #                 # Avoids double counting
    #                 #(elements sharing the same parent are removed during clustering)
    #                 if newel and not any(newel == el for el in newElements[:]):
    #                     newElements.append(newel)
    #                     added = True
    #
    #         # Check for invisible compressed topologies (look for effective
    #         # LSP, such as LSP + neutrino = LSP')
    #         if doInvisible:
    #             for element in newElements:
    #                 newel = element.invisibleCompress()
    #                 # Avoids double counting
    #                 #(elements sharing the same parent are removed during clustering)
    #                 if newel and not any(newel == el for el in newElements[:]):
    #                     newElements.append(newel)
    #                     added = True
    #
    #     newElements.pop(0)  # Remove original element
    #     return newElements
    #
    # def removeVertex(self, ibr, iv):
    #     """
    #     Remove vertex iv located in branch ibr.
    #     The "vertex-mother" in BSMparticles and (SM) particles in the vertex
    #     are removed from the branch. The vertex index corresponds
    #     to the BSM decay (iv = 0 will remove the first BSM particle,...)
    #
    #     :parameter ibr: Index of branch (int)
    #     :parameter iv: Index of vertex in branch ibr (int)
    #
    #     """
    #
    #     self.branches[ibr].removeVertex(iv)
    #
    # def massCompress(self, minmassgap):
    #     """
    #     Perform mass compression.
    #
    #     :parameter minmassgap: value (in GeV) of the maximum
    #                            mass difference for compression
    #                            (if mass difference < minmassgap -> perform mass compression)
    #     :returns: compressed copy of the element, if two masses in this
    #               element are degenerate; None, if compression is not possible;
    #     """
    #
    #     newelement = self.copy()
    #     newelement.motherElements = [self]
    #
    #     #Loop over branches and look for small mass differences
    #     for ibr, branch in enumerate(newelement.branches):
    #         #Get mass differences
    #
    #         removeVertices = []
    #         for i, mom in enumerate(branch.oddParticles[:-1]):
    #             massDiff = mom.mass - branch.oddParticles[i+1].mass
    #             #Get vertices which have deltaM < minmassgap and the mother is prompt:
    #             if massDiff < minmassgap and mom.isPrompt():
    #                 removeVertices.append(i)
    #         #Remove vertices till they exist:
    #         while removeVertices:
    #             newelement.removeVertex(ibr, removeVertices[0])
    #             branch = newelement.branches[ibr]
    #             removeVertices = []
    #             for i, mom in enumerate(branch.oddParticles[:-1]):
    #                 massDiff = mom.mass - branch.oddParticles[i+1].mass
    #                 #Get vertices which have deltaM < minmassgap and the mother is prompt:
    #                 if massDiff < minmassgap and mom.isPrompt():
    #                     removeVertices.append(i)
    #
    #     for ibr, branch in enumerate(newelement.branches):
    #         if branch.vertnumb != self.branches[ibr].vertnumb:
    #             newelement.sortBranches()
    #             return newelement
    #
    #     #New element was not compressed, return None
    #     return None
    #
    # def invisibleCompress(self):
    #     """
    #     Perform invisible compression.
    #
    #     :returns: compressed copy of the element, if element ends with invisible
    #               particles; None, if compression is not possible
    #     """
    #
    #     newelement = self.copy()
    #     newelement.motherElements = [self]
    #
    #     # Loop over branches
    #     for branch in newelement.branches:
    #         if not branch.evenParticles:
    #             continue
    #         #Check if the last decay should be removed:
    #         neutralSM = all(ptc.isMET() for ptc in branch.evenParticles[-1])
    #         neutralDecay = neutralSM and branch.oddParticles[-1].isMET()
    #         #Check if the mother can be considered MET:
    #         neutralBSM = (branch.oddParticles[-2].isMET()
    #                      or branch.oddParticles[-2].isPrompt())
    #         if neutralBSM and neutralDecay:
    #             removeLastVertex = True
    #         else:
    #             removeLastVertex = False
    #
    #         while len(branch.oddParticles) > 1 and removeLastVertex:
    #             bsmMom = branch.oddParticles[-2]
    #             effectiveDaughter = Particle(label='inv', mass=bsmMom.mass,
    #                                         eCharge=0, colordim=1,
    #                                         totalwidth=branch.oddParticles[-1].totalwidth,
    #                                         Z2parity=bsmMom.Z2parity, pdg=bsmMom.pdg)
    #             branch.removeVertex(len(branch.oddParticles)-2)
    #             #For invisible compression, keep an effective mother which corresponds to the invisible
    #             #daughter, but with the mass of the parent.
    #             branch.oddParticles[-1] = effectiveDaughter
    #             #Re-check if the last decay should be removed:
    #             if not branch.evenParticles:
    #                 continue
    #             neutralSM = all(ptc.isMET() for ptc in branch.evenParticles[-1])
    #             neutralDecay = neutralSM and branch.oddParticles[-1].isMET()
    #             neutralBSM = branch.oddParticles[-2].isMET()
    #             if neutralBSM and neutralDecay:
    #                 removeLastVertex = True
    #             else:
    #                 removeLastVertex = False
    #
    #     for ibr, branch in enumerate(newelement.branches):
    #         if branch.vertnumb != self.branches[ibr].vertnumb:
    #             newelement.sortBranches()
    #             return newelement
    #
    #     #New element was not compressed, return None
    #     return None
    #
    # def hasTopInList(self, elementList):
    #     """
    #     Check if the element topology matches any of the topologies in the
    #     element list.
    #
    #     :parameter elementList: list of elements (Element objects)
    #     :returns: True, if element topology has a match in the list, False otherwise.
    #     """
    #     if type(elementList) != type([]) or len(elementList) == 0:
    #         return False
    #     for element in elementList:
    #         if type(element) != type(self):
    #             continue
    #         info1 = self.getCanonName()
    #         info2 = element.getCanonName()
    #         info2B = element.switchBranches().getCanonName()
    #         if info1 == info2 or info1 == info2B:
    #             return True
    #     return False
