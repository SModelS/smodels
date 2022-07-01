"""
.. module:: element
   :synopsis: Module holding the Element class and its methods.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.theory.tree import Tree
from smodels.theory import crossSection
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.theory.particle import Particle


class Element(object):
    """
    An instance of this class represents an element.
    This class possesses a pair of branches and the element weightList
    (cross-section * BR).
    """

    def __init__(self, info=None, finalState=None,
                 intermediateState=None, model=None, sort=True):
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

        :parameter sort: If False, it will not sort the element.
        """

        if info is None:
            self.tree = None
        elif isinstance(info, Tree):
            self.tree = info.copyTree()
        else:
            self.tree = Tree(info=info, finalState=finalState,
                             intermediateState=intermediateState,
                             model=model)
        self.weightList = crossSection.XSectionList()
        self.motherElements = [self]  # The motheElements includes self to keep track of merged elements
        self.elID = 0
        self.coveredBy = set()
        self.testedBy = set()

        # Check graph consistency and sort it
        # (must be a rooted tree with each node having a single parent):
        if self.tree and self.tree.number_of_nodes():
            self.tree.setGlobalProperties()
            if sort:
                self.sort()
            self.tree.checkConsistency()

    def sort(self):
        """
        Sort the tree according to the canonName and then particle. For each node,
        all the daughters are sorted according to their canonName.
        """

        self.tree.sort()
        self.tree.numberNodes()

    def compareTo(self, other):
        """
        Compares the element with other.
        Uses the topology name (Tree canonincal name) to identify isomorphic topologies (trees).
        If trees are not isomorphic, compare topologies using the topology name.
        Else  check if nodes (particles) are equal. For nodes with the same canonName
        and the same parent, compare against all permutations.
        If the canonical names match, but the particle differs, returns the particle
        comparison.

        :param other:  element to be compared (Element object)

        :return: (cmp,otherSorted), where cmp = -1 if self < other, 0 if self == other, +1, if self > other and otherSorted is None if cmp != 0 or other sorted according to the way it matched self.
        """

        if not isinstance(other, Element):
            return -1, None

        # make sure the topology names have been computed:
        canonName = self.tree.canonName
        otherName = other.tree.canonName
        if canonName != otherName:
            if canonName > otherName:
                return 1, None
            else:
                return -1, None

        # Recursively compare the nodes:
        cmp, newTree = self.tree.compareTreeTo(other.tree)

        if cmp == 0:  # Elements matched, return copy of other with tree sorted
            otherNew = other.copy(emptyTree=True)
            otherNew.tree = newTree
        else:
            otherNew = None
        return cmp, otherNew

    def __eq__(self, other):
        cmp, otherSorted = self.compareTo(other)
        return (cmp == 0)

    def __lt__(self, other):
        cmp, otherSorted = self.compareTo(other)
        return (cmp < 0)

    def __gt__(self, other):
        return not (self < other)

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
                   for node in self.tree.nodes]
            return val
        except AttributeError:
            raise AttributeError("Neither element nor particles have attribute ``%s''" % attr)

    def __str__(self):
        """
        Create the element to a string.

        :returns: string representation of the element
        """

        return self.tree.treeToString()

    def __repr__(self):

        return self.__str__()

    def __add__(self, other):
        """
        Adds two elements. Should only be used if the elements
        have the same topologies. The element weights are added and their
        odd and even particles are combined.
        """

        if not isinstance(other, Element):
            raise TypeError("Can not add an Element object to %s" % type(other))
        elif self.canonName != other.canonName:
            raise SModelSError("Can not add elements with distinct topologies")

        newEl = self.__class__()
        newEl.motherElements = self.motherElements[:] + other.motherElements[:]
        newEl.weightList = self.weightList + other.weightList
        newEl.tree = self.tree + other.tree
        newEl.sort()

        return newEl

    def __radd__(self, other):
        """
        Adds two elements. Only elements with the same
        topology can be combined.
        """

        return self.__add__(other)

    def __iadd__(self, other):
        """
        Combine two elements. Should only be used if the elements
        have the same topologies. The element weights are added and their
        odd and even particles are combined.
        """

        if not isinstance(other, Element):
            raise TypeError("Can not add an Element object to %s" % type(other))
        elif self.canonName != other.canonName:
            raise SModelSError("Can not add elements with distinct topologies")

        self.motherElements += other.motherElements[:]
        self.weightList += other.weightList
        self.tree = self.tree + other.tree

        return self

    def drawTree(self, particleColor='lightcoral',
                 smColor='skyblue',
                 pvColor='darkgray',
                 nodeScale=4, labelAttr=None, attrUnit=None):
        """
        Draws Tree using matplotlib.

        :param tree: tree to be drawn
        :param particleColor: color for particle nodes
        :param smColor: color used for particles which have the isSM attribute set to True
        :param pvColor: color for primary vertex
        :param nodeScale: scale size for nodes
        :param labelAttr: attribute to be used as label. If None, will use the string representation
                          of the node object.
        :param attrUnit: Unum object with the unit to be removed from label attribute (if applicable)

        """

        return self.tree.draw(particleColor=particleColor,
                              smColor=smColor,
                              pvColor=pvColor, nodeScale=nodeScale,
                              labelAttr=labelAttr, attrUnit=attrUnit)

    def copy(self, emptyTree=False):
        """
        Create a copy of self.
        Faster than deepcopy.

        :param emptyTree: If True, creates an element without a
                          tree.

        :returns: copy of element (Element object)
        """

        # Allows for derived classes (like inclusive classes)
        if emptyTree:
            newel = self.__class__()
        else:
            newel = self.__class__(info=self.tree.copyTree())
        newel.weightList = self.weightList.copy()
        newel.motherElements = self.motherElements[:]
        newel.elID = self.elID
        newel.coveredBy = set(list(self.coveredBy)[:])
        newel.testedBy = set(list(self.testedBy)[:])

        return newel

    def getFinalStates(self):
        """
        Get the list of particles which have not decayed (within the element)

        :returns: list of ParticleNode objects
        """

        finalStates = self.tree.getFinalStates()
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

        ancestorsDict = {igen + 1: []}
        for mother in self.motherElements:
            if mother is self:
                continue
            ancestorsDict[igen + 1].append(mother)
            for jgen, elList in mother._getAncestorsDict(igen + 1).items():
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

    @property
    def canonName(self):
        """
        Returns the canonincal name for the tree. If not defined,
        compute it.
        """

        canonName = self.tree.canonName
        if not canonName:
            self.tree.setGlobalProperties()
            canonName = self.tree.canonName

        return canonName

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
                    # (elements sharing the same parent are removed during clustering)
                    if (newel is not None) and (newel not in newElements):
                        newElements.append(newel)
                        added = True

            # Check for invisible compressed topologies (look for effective
            # LSP, such as LSP + neutrino = LSP')
            if doInvisible:
                for element in newElements:
                    newel = element.invisibleCompress()
                    # Avoids double counting
                    # (elements sharing the same parent are removed during clustering)
                    if (newel is not None) and (newel not in newElements):
                        newElements.append(newel)
                        added = True

        newElements.pop(0)  # Remove original element

        return newElements

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

        tree = newelement.tree
        root = tree.root.node
        # Loop over nodes from root to leaves:
        for mom, daughters in list(tree.bfs_successors(root)):
            if mom == root:  # Skip primary vertex
                continue
            if mom not in tree.successors:  # In case the mother has been removed by compression
                continue
            # Convert node index to node object
            mom = tree.nodesMapping[mom]
            if not mom.particle.isPrompt():  # Skip long-lived
                continue
            bsmDaughter = []
            smDaughters = []
            for d in daughters:
                # Convert to node object
                d = tree.nodesMapping[d]
                # Split daughters into final states SM and others (BSM)
                if hasattr(d, 'isSM') and d.isSM and tree.out_degree(d.node) == 0:
                    smDaughters.append(d)
                else:
                    bsmDaughter.append(d)

            # Skip decays to multiple BSM particles or to SM particles only
            if len(bsmDaughter) != 1:
                continue
            bsmDaughter = bsmDaughter[0]

            # Check mass difference:
            massDiff = mom.mass - bsmDaughter.mass
            if massDiff > minmassgap:
                continue

            # Get grandmother:
            gMomIndex = tree.predecessors[mom.node]
            # Remove mother and all SM daughters and mom:
            removeIndices = [d.node for d in smDaughters] + [mom.node]
            tree.remove_nodes_from(removeIndices)

            # Attach BSM daughter to grandmother:
            tree.add_edge(tree.nodesMapping[gMomIndex], bsmDaughter)

        # Recompute the canonical name and
        newelement.tree.setGlobalProperties()
        newelement.sort()

        # If element was not compressed, return None
        if newelement.canonName == self.canonName:
            return None
        else:
            return newelement

    def invisibleCompress(self):
        """
        Perform invisible compression.

        :returns: compressed copy of the element, if element ends with invisible
                  particles; None, if compression is not possible
        """

        newelement = self.copy()
        newelement.motherElements = [self]
        keepCompressing = True
        # Check for compression until tree can no longer be compressed
        previousName = self.canonName

        while keepCompressing:
            tree = newelement.tree
            root = tree.root.node
            # Loop over nodes:
            for mom, daughters in tree.bfs_successors(root):
                if mom == root:  # Skip primary vertex
                    continue
                # Skip node if its daughters are not stable
                if any(tree.out_degree(d) != 0 for d in daughters):
                    continue
                # Convert indices to node objects:
                mom = tree.nodesMapping[mom]
                daughters = [tree.nodesMapping[d] for d in daughters]
                # Check if all daughters can be considered MET
                neutralDaughters = all(d.particle.isMET() for d in daughters)
                # Check if the mother is MET or prompt:
                neutralMom = (mom.isMET() or mom.isPrompt())

                # If mother and daughters are neutral, remove daughters
                if (not neutralDaughters) or (not neutralMom):
                    continue

                removeIndices = [d.node for d in daughters]
                tree.remove_nodes_from(removeIndices)
                # Replace mom particle by invisible (generic) particle
                # with its width equal to the maximum width amongst the daughters
                maxWidth = max([d.totalwidth for d in daughters])
                invParticle = Particle(label='inv', mass=mom.mass,
                                       eCharge=0, colordim=1,
                                       _isInvisible=True,
                                       totalwidth=maxWidth,
                                       pdg=mom.pdg, isSM=mom.isSM)
                mom.particle = invParticle

                # For safety break loop since tree structure changed
                break

            # Recompute the canonical name and
            newelement.tree.setGlobalProperties()
            # If iteration has not changed element, break loop
            name = newelement.canonName
            if name == previousName:
                keepCompressing = False
            else:
                keepCompressing = True
                previousName = name
        newelement.sort()

        # If element was not compressed, return None
        if newelement.canonName == self.canonName:
            return None
        else:
            return newelement

    def compressToFinalStates(self):
        """
        Compress the elements to its final states. After the compression
        the element will hold a tree with one root (PV). The root's daughters
        are the final state nodes of self.

        :returns: compressed copy of the element.

        """

        newElement = self.copy(emptyTree=True)
        newElement.tree = self.tree.compressToFinalStates()

        return newElement
