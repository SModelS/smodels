"""
.. module:: tree
   :synopsis: This is a class for describing Simplified Model Topologies
              used for decomposing BSM models.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.genericSMS import GenericSMS
from smodels.theory.particle import Particle
from smodels.theory import crossSection
from itertools import product
from smodels.tools.physicsUnits import fb

class TheorySMS(GenericSMS):
    """
    A class for describing Simplified Model Topologies generated by the decompostion
    of full BSM models.
    """

    def __init__(self):
        """
        Initialize basic attributes.
        """

        GenericSMS.__init__(self)

        # Include additional attributes
        self._maxWeight = None
        self._weightList = None
        self._sorted = False
        # List of SMS topologies which could have generated self (such as during compression)
        self.ancestors = [self]
        self._allAncestors = None
        self.smsID = 0  # SMS identifier
        # Type of analyses which have SMS matching self:
        self.coveredBy = set()
        # Type of analyses which have SMS matching self and for
        # which the physical parameters are covered
        self.testedBy = set()

    def __eq__(self, other):
        """
        SMS equality based on the compareTo method.

        :parameter other: TheorySMS object

        :return: True if objects are equivalent.
        """

        cmp = self.compareTo(other)

        return (cmp == 0)

    def __lt__(self, other):
        """
        SMS comparison based on the compareTo method.

        :parameter other: TheorySMS object

        :return: True if self < other, False otherwise.
        """

        cmp = self.compareTo(other)

        return (cmp < 0)

    def __gt__(self, other):
        """
        SMS comparison based on the compareTo method.

        :parameter other: TheorySMS object

        :return: True if self > other, False otherwise.
        """
        return not (self < other)

    def __getattr__(self, attr):
        """
        If the attribute has not been defined for self
        try to fetch it from its nodes
        :param attr: Attribute name
        :return: Attribute value
        """

        # If calling another special method, return default (required for pickling)
        if attr.startswith('__') and attr.endswith('__'):
            return object.__getattr__(attr)

        try:
            val = [getattr(node, attr) if str(node) != 'PV' else None
                   for node in self.nodes]
            return val
        except AttributeError:
            raise AttributeError("Neither SMS nor nodes have attribute ``%s''" % attr)

    def __add__(self, other):
        """
        Combines the equivalent nodes (add particles) and the weightLists.
        Can only be done if self and other have the same topology and ordering.

        :param other: TheorySMS object

        :return: new SMS with the combined particles (TheorySMS object)
        """

        if not isinstance(other,TheorySMS):
            raise SModelSError("Can not compare TheorySMS and %s" %str(type(other)))
        elif self.canonName != other.canonName:
            raise SModelSError("Can not add elements with distinct topologies")

        # Make a copy
        newSMS = self.copy()
        # Add nodes from other
        newSMS.addNodesFrom(other)
        # Add other attributes
        newSMS.ancestors = self.ancestors[:] + other.ancestors[:]
        newSMS._weightList = self._weightList + other._weightList
        newSMS._maxWeight = self._maxWeight + other._maxWeight
        if other._allAncestors is not None:
            newSMS._allAncestors += other._allAncestors[:]

        return newSMS

    def setGlobalProperties(self,sort=True,canonName=True,weight=True):
        """
        Compute and set global properties for the SMS
        (sort, renumber nodes, compute the canonical name and its total weight).
        Should only be called once the SMS will no longer be modified.

        :parameter canonName: If True, compute its canonical name
        :parameter sort: If True, sort the SMS
        :parameter weight: If True, compute its total weight
        """

        if canonName:
            self._canonName = self.getCanonName()
        if sort:
            self.sort(force=True)
        if weight:
            self._weightList = self.computeWeightList()

    def compareTo(self, other):
        """
        Compare self to other.
        If the SMS are not sorted, sort them and then do a direct comparison of each
        node with the same nodeIndex.

        :param other: TheorySMS object to be compared against self

        :return: 0, if objects are equal, -1 if self < other, 1 if other > sekf
        """

        if not isinstance(other,TheorySMS):
            raise SModelSError("Can not compare TheorySMS and %s" %str(type(other)))

        # Make sure the SMS are sorted
        # (if it has already been sorted, it does nothing)
        self.sort()
        other.sort()

        # Compare the nodes
        #(the first node is the dummy PV node, so only the canonical name is compared)
        for nodeIndex in self.nodeIndices:
            nodeA = self.indexToNode(nodeIndex)
            nodeB = other.indexToNode(nodeIndex)
            cmp = nodeA.compareTo(nodeB)
            if cmp != 0:
                return cmp

    def copy(self):
        """
        Returns a shallow copy of self.

        : return: TheorySMS object
        """

        newSMS = TheorySMS()
        newSMS._successors.update({n: daughters[:]
                                       for n, daughters in self._successors.items()})
        newSMS._predecessors = {k: v for k, v in self._predecessors.items()}
        newSMS._nodesMapping = {n: node for n, node in self.nodesMapping.items()}
        newSMS._rootIndex = self._rootIndex
        newSMS._canonName = self._canonName
        newSMS._maxWeight = self._maxWeight
        newSMS._sorted = self._sorted
        newSMS.weightList = self.weightList.copy()
        newSMS.ancestors = self.ancestors[:]
        if self._allAncestors is not None:
            newSMS._allAncestors = self._allAncestors[:]
        newSMS.smsID = self.smsID
        newSMS.coveredBy = set(list(self.coveredBy)[:])
        newSMS.testedBy = set(list(self.testedBy)[:])

        return newSMS

    def addNodesFrom(self, other):
        """
        Combines the nodes (add particles) in equivalent nodes in each trees.
        Can only be done if the trees have the same topology and ordering.

        :param other: TheorySMS object
        """

        if other.canonName != self.canonName:
            raise SModelSError("Can not add trees with distinct topologies")

        nodesB = other.nodes
        newMapping = {}
        for n in self.nodeIndices:
            newMapping[n] = self.indexToNode(n) + other.indexToNode(n)

        # Just change the nodeIndex -> node mapping:
        self.relabelNodes(nodeObjectDict=newMapping)

    def attachDecay(self, motherIndex, decayNodes, copy=True):
        """
        Attaches a decay to self. If copy = True, returns a copy of self
        with the decay attached.

        : param decay: A tuple containing the mother node (carring the BR as the nodeWeight)
                      and the daughter nodes.

        : param copy: if True, return a copy of self, with the decay attached.

        : return: new tree with the other composed with self.
        """

        motherNode = decayNodes[0]
        daughterNodes = decayNodes[1]

        if copy:
            newSMS = self.copy()
        else:
            newSMS = self

        oldMother = self.indexToNode(motherIndex)
        # Update maximum weight
        if newSMS._maxWeight:
            oldMotherWeight = oldMother.nodeWeight
            if oldMotherWeight:
                newSMS._maxWeight = newSMS._maxWeight/oldMotherWeight
            newSMS._maxWeight = newSMS._maxWeight*motherNode.nodeWeight

        # Update mother node:
        self.updateNodeObjects({motherIndex : motherNode})
        # Add daughter nodes and edges:
        for d in daughterNodes:
            dIndex = self.add_node(d)
            newSMS.add_edge(motherIndex, dIndex)

        # The tree is no longer sorted
        newSMS._sorted = False

        return newSMS

    def getFinalStates(self, nodeIndex=None):
        """
        Get the list of particles which have not decayed (appear at the top of the tree).
        If nodeIndex is defined, get the final states generated by the cascade decay of the corresponding node. It also sets the finalState attribute for the node.

        :param nodeIndex: Node index for which to get the final states.

        :returns: list of Particle objects
        """

        if nodeIndex is None:
            nodeIndex = self.rootIndex

        node = self.indexToNode(nodeIndex)
        if node.finalStates is not None:
            return node.finalStates

        # For leaves, the final state is themselves:
        if self.out_degree(nodeIndex) == 0:
            finalStates = [node.particle]
        else:
            finalStates = []
            for dIndex in self.daughterIndices(nodeIndex):
                finalStates += self.getFinalStates(dIndex)

        # Set final state attribute for node:
        node.finalStates = sorted(finalStates)

        return node.finalStates

    @property
    def weightList(self):
        """
        Returns the SMS weight list (production cross-sections*BRs).
        It does not include the weight of final state (undecayed) particles.

        :return: CrossSectionList object
        """

        if self._weightList is None:
            self._weightList = self.computeMaxWeightList()

        return self._weightList

    def computeWeightList(self):
        """
        Computes the SMS weight (production cross-section*BRs).
        It does not include the weight of final state (undecayed) particles.

        :return: CrossSectionList object
        """

        root = self.root
        prodXSec = root.nodeWeight
        maxXSec = prodXSec.getMaxXsec().asNumber(fb)
        maxWeight = self.getMaxWeight()
        brs = maxWeight/maxXSec
        weightList = prodXSec*brs

        return weightList

    @property
    def maxWeight(self):
        """
        Returns the SMS maximum weight (max production cross-section*BRs).
        Does not include the weight of final state (undecayed) particles.
        If it has not yet been computed, compute it.

        :return: SMS weight in fb (float).
        """

        if self._maxWeight is None:
            self._maxWeight = self.computeMaxWeight()

        return self._maxWeight

    def computeMaxWeight(self):
        """
        Computes the SMS maximum weight (max production cross-section*BRs).
        Does not include the weight of final state (undecayed) particles.

        :return: SMS weight in fb (float).
        """

        # Get maximum production cross-section
        root = self.root
        prodXSec = root.nodeWeight
        maxXsec = prodXSec.getMaxXsec().asNumber(fb)

        # Get product of branching ratios:
        brs = 1.0
        for nodeIndex in self.nodeIndices:
            if nodeIndex == self.rootIndex:
                continue
            if self.out_degree(nodeIndex) == 0:
                continue
            node = self.indexToNode(nodeIndex)
            brs *= node.nodeWeight

        weight = maxXsec*brs

        return weight

    def sort(self, nodeIndex=None, force=False):
        """
        Sort subtree of self generated by nodeIndex.
        If nodeIndex is None, sort the tree and re-number the nodes
        according to the bfs order (after sorting).
        If the self is already tagged as sorted and force = False,
        do nothing.

        : param nodeIndex: Node index
        : param force: If True, will sort even if self is tagged as sorted.
        """

        # If tree is already sorted, do nothing
        if self._sorted and not force:
            return

        if nodeIndex is None:
            cName = self.canonName  # Just to make sure canonName is defined
            if cName is None:
                return
            nodeIndex = self.rootIndex

        daughters = self.daughterIndices(nodeIndex)
        if daughters:
            for d in daughters:
                self.sort(d, force=force)
            sorted_daughters = self.sortSubTrees(daughters)


            # Remove nodeIndex -> daughters edges
            self.remove_edges(product([nodeIndex],daughters))
            # Add edges with the correct ordering:
            self.add_edges(product([nodeIndex],sorted_daughters))

        # Finally, after sorting the subtrees,
        # make sure the nodes are sorted according
        # to the generations (breadth-first search)
        if nodeIndex == self.rootIndex:
            self.bfs_sort()
            self.numberNodes()
            # Tag the tree as sorted
            self._sorted = True

    def sortSubTrees(self, subtreeList):
        """
        Sorts a list of subtrees of self generated by the nodes
        in subtreeList.

        : param subtreeList: List of node indices to be considered as roots
                            of the subtrees.

        : return: Sorted list of node indices.
        """

        if len(subtreeList) == 1 or len(subtreeList) == 0:
            return subtreeList

        # First group subtrees by their canonical name:
        nameDict = {}
        for nodeIndex in subtreeList:
            cName = self.indexToNode(nodeIndex).canonName
            if cName not in nameDict:
                nameDict[cName] = [nodeIndex]
            else:
                nameDict[cName].append(nodeIndex)

        # Now sort subtrees with common canonical names:
        sorted_trees = []
        for cName in sorted(nameDict.keys()):
            # Within equal canonincal names sort daughters by the generated subtrees
            sorted_trees += self.sortCommonSubTrees(nameDict[cName])

        return sorted_trees

    def sortCommonSubTrees(self, subtreeList):
        """
        Sorts a list of subtrees of self generated by the nodes
        in subtreeList using a quicksort algorithm.
        All the subtrees should have a common topology
        (same canonical name).

        : param subtreeList: List of node indices

        : return: Sorted list of node indices.
        """

        if len(subtreeList) == 1 or len(subtreeList) == 0:
            return subtreeList

        s = subtreeList[:]
        pivot = s[0]
        i = 0
        for j in range(len(s)-1):
            if self.compareSubTrees(s[j+1], pivot) < 0:
                s[j+1], s[i+1] = s[i+1], s[j+1]
                i += 1
        s[0], s[i] = s[i], s[0]
        first_part = self.sortCommonSubTrees(s[:i])
        second_part = self.sortCommonSubTrees(s[i+1:])
        first_part.append(s[i])

        sortedList = first_part + second_part

        return sortedList

    def compareSubTrees(self, n1, n2):
        """
        Compare the subtrees generated by the nodes n1 and n2.

        : param n1: Node index for the root of subtree1
        : param n2: Node index for the root of subtree2

        : return: 0, if subtrees are equal, -1 if subtree1 < subtree2, 1 if subtree1 > subtree2
        """

        root1 = self.indexToNode(n1)
        root2 = self.indexToNode(n2)
        cmp = root1.compareTo(root2)
        if cmp != 0:
            return cmp

        # For inclusive nodes always return True (once nodes are equal)
        if root1.isInclusive or root2.isInclusive:
            return 0

        daughters1 = self.daughterIndices(n1)
        daughters2 = self.daughterIndices(n2)
        # If nodes are leaves, return 0
        if len(daughters1) == len(daughters2) == 0:
            return 0

        # Check if the daughters from n2 match the ones from n1:
        # (the daughters should be sorted at this point)
        for i1, d1 in enumerate(daughters1):
            d2 = daughters2[i1]
            cmp = self.compareSubTrees(d1, d2)
            if cmp != 0:
                return cmp
        return 0

    def compressToFinalStates(self):
        """
        Compress the SMS to its final states. After the compression
        the tree will have one root(PV). The root's daughters
        are the final state nodes of self.

        : returns: compressed copy of the Tree.
        """

        newTree = self.copy()
        # Get root:
        root = newTree.root
        # Get final state nodes
        fsNodes = [self.indexToNode(n) for n in self.nodeIndices if
                    self.out_degree(n) == 0]

        # Remove all nodes
        newTree.clear()
        # Add root and final state nodes
        rootIndex = newTree.add_node(root)
        fsIndices = newTree.add_nodes_from(fsNodes)
        # Add root > fsNode edges:
        edges = product([rootIndex],fsIndices)
        newTree.add_edges_from(edges)
        newTree._canonName = newTree.getCanonName()
        newTree.sort()

        return newTree

    def _getAncestorsDict(self, igen=0):
        """
        Returns a dictionary with all the ancestors
        of the SMS. The dictionary keys are integers
        labeling the generation (number of generations away from self)
        and the values are a list of SMS objects (ancestors) for that generation.
        igen is used as the counter for the initial generation.
        The output is also stored in self._ancestorsDict for future use.

        :param igen: Auxiliary integer indicating to which generation self belongs.

        :return: Dictionary with generation index as key and ancestors as values
                 (e.g. {igen+1 : [mother1, mother2], igen+2 : [grandmother1,..],...})
        """

        ancestorsDict = {igen + 1: []}
        for ancestor in self.ancestors:
            if ancestor is self:
                continue
            ancestorsDict[igen + 1].append(ancestor)
            for jgen, elList in ancestor._getAncestorsDict(igen + 1).items():
                if jgen not in ancestorsDict:
                    ancestorsDict[jgen] = []
                ancestorsDict[jgen] += elList

        return ancestorsDict

    def getAncestors(self):
        """
        Get a list of all the ancestors of self.
        The list is ordered so the mothers appear first, then the grandmother,
        then the grandgrandmothers,...

        :return: A list of SMS objects containing all the ancestors sorted by generation.
        """

        # Check if the ancestors have already been obtained (performance gain)
        if hasattr(self, '_allAncestors'):
            return self._allAncestors

        ancestorsDict = self._getAncestorsDict()
        allAncestors = []

        for igen in sorted(ancestorsDict.keys()):
            allAncestors += ancestorsDict[igen]

        self._allAncestors = allAncestors

        return self._allAncestors

    def isRelatedTo(self, other):
        """
        Checks if self has any common ancestors with other or they share
        ancestors. Returns True if self and other have at least one ancestor in common
        or are the same object, otherwise returns False.

        :return: True/False
        """

        if self is other:
            return True

        ancestorsA = set([id(sms) for sms in self.getAncestors()])
        ancestorsB = set([id(sms) for sms in other.getAncestors()])

        if ancestorsA.intersection(ancestorsB):
            return True
        else:
            return False

    def setTestedBy(self, resultType):
        """
        Tag the sms as tested by the result type.
        It also recursively tags all of its ancestors.

        :param resultType: String describing the type of result (e.g. 'prompt', 'displaced')
        """

        self.testedBy.add(resultType)
        for ancestor in self.getAncestors():
            ancestor.testedBy.add(resultType)

    def setCoveredBy(self, resultType):
        """
        Tag the sms as covered by the result type
        (it is tested AND its parameters are within the result grid).
        It also recursively tags all of its ancestors.

        :param resultType: String describing the type of result (e.g. 'prompt', 'displaced')
        """

        self.coveredBy.add(resultType)
        for ancestor in self.getAncestors():
            ancestor.coveredBy.add(resultType)

    def compress(self, doCompress, doInvisible, minmassgap):
        """
        Keep compressing the original SMS and the derived ones till they
        can be compressed no more.

        :parameter doCompress: if True, perform mass compression
        :parameter doInvisible: if True, perform invisible compression
        :parameter minmassgap: value (in GeV) of the maximum
                               mass difference for compression
                               (if mass difference < minmassgap, perform mass compression)
        :returns: list with the compressed SMS (TheorySMS objects)
        """

        if not doCompress and not doInvisible:
            return []

        added = True
        newList = [self]
        # Keep compressing the new topologies generated so far until no new
        # compressions can happen:
        while added:
            added = False
            # Check for mass compressed topologies
            if doCompress:
                for sms in newList:
                    newSMS = sms.massCompress(minmassgap)
                    # Avoids double counting
                    # (elements sharing the same parent are removed during clustering)
                    if (newSMS is not None) and (newSMS not in newList):
                        newList.append(newSMS)
                        added = True

            # Check for invisible compressed topologies (look for effective
            # LSP, such as LSP + neutrino = LSP')
            if doInvisible:
                for sms in newList:
                    newSMS = sms.invisibleCompress()
                    # Avoids double counting
                    # (elements sharing the same parent are removed during clustering)
                    if (newSMS is not None) and (newSMS not in newList):
                        newList.append(newSMS)
                        added = True

        newList.pop(0)  # Remove original element

        return newList

    def massCompress(self, minmassgap):
        """
        Perform mass compression.
        It is only done if there is one decay of type BSM_i -> BSM_j + (any number of SM),
        where mass(BSM_i) - mass(BSM_j) < minmassgap AND BSM_i have a prompt decay.

        :parameter minmassgap: value (in GeV) of the maximum
                               mass difference for compression
                               (if mass difference < minmassgap -> perform mass compression)
        :returns: compressed copy of self, if two masses in the
                  SMS can be considered degenerate; None, if compression is not possible;
        """

        newSMS = self.copy()
        newSMS.ancestors = [self]

        # Loop over nodes from root to leaves:
        for momIndex, daughterIndices in self.genIndexIterator():
            if momIndex is newSMS.rootIndex:  # Skip primary vertex
                continue
            if momIndex not in newSMS.nodeIndices:  # In case the mother has been removed by compression
                continue
            # Convert node index to node object
            mom = newSMS.indexToNode(momIndex)
            if mom.isSM:  # Only compress BSM states
                continue
            bsmDaughter = []
            smDaughters = []
            for dIndex in daughterIndices:
                # Convert to node object
                d = newSMS.indexToNode(dIndex)
                # Split daughters into final states SM and others (BSM)
                if hasattr(d, 'isSM') and d.isSM and newSMS.out_degree(dIndex) == 0:
                    smDaughters.append(dIndex)
                else:
                    # Compute mass different between BSM daughter and mom
                    massDiff = mom.mass - d.mass
                    bsmDaughter.append(dIndex)

            # Skip decays to multiple BSM particles or to SM particles only
            if len(bsmDaughter) != 1:
                continue
            bsmDaughter = bsmDaughter[0]

            # Skip if mass difference is above minimum or if the parent is long-lived
            if massDiff > minmassgap or not mom.particle.isPrompt():
                continue

            # Get grandmother:
            gMomIndex = newSMS.parentIndex(momIndex)
            # Remove mother and all SM daughters:
            removeIndices = smDaughters + [momIndex]
            newSMS.remove_nodes_from(removeIndices)

            # Attach BSM daughter to grandmother:
            newSMS.add_edge(gMomIndex, bsmDaughter)

        # If no compression was made, return None
        if self.number_of_nodes() == newSMS.number_of_nodes():
            return None

        # Recompute the global properties (except for the weightList)
        # and sort the new SMS
        newSMS.setGlobalProperties(weight=False)

        return newSMS

    def invisibleCompress(self):
        """
        Perform invisible compression.
        It is done if there is a decay of the type BSM > xxx
        if all daughters are leaves and can be considered MET and the mother can
        be considered MET OR decays promptly.

        :returns: compressed copy of the SMS, if element ends with invisible
                  particles; None, if compression is not possible
        """

        newSMS = self.copy()
        newSMS.ancestors = [self]

        while True:
            # Loop over nodes:
            nNodes = newSMS.number_of_nodes()
            for momIndex, daughterIndices in newSMS.genIndexIterator():
                if momIndex == newSMS.rootIndex:  # Skip primary vertex
                    continue
                # Skip node if its daughters are not stable
                if any(newSMS.out_degree(d) != 0 for d in daughterIndices):
                    continue
                # Convert indices to node objects:
                mom = newSMS.indexToNode(momIndex)
                daughters = newSMS.indexToNode(daughterIndices)
                # Check if all daughters can be considered MET
                neutralDaughters = all(d.particle.isMET() for d in daughters)
                # Check if the mother is MET or prompt:
                neutralMom = (mom.isMET() or mom.isPrompt())

                # If mother and daughters are neutral, remove daughters
                if (not neutralDaughters) or (not neutralMom):
                    continue

                newSMS.remove_nodes_from(daughterIndices)
                # Replace mother particle by invisible (generic) particle
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

            if newSMS.number_of_nodes() == nNodes:
                break  # No compressions could be made, stop


        # Recompute the global properties (except for the weightList)
        # and sort the new SMS
        newSMS.setGlobalProperties(weight=False)

        return newSMS


