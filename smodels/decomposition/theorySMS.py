"""
.. module:: tree
   :synopsis: This is a class for describing Simplified Model Topologies
              used for decomposing BSM models.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


from smodels.decomposition.exceptions import SModelSDecompositionError as SModelSError
from smodels.base.genericSMS import GenericSMS
from smodels.base.particle import Particle
from smodels.base import crossSection
from smodels.base.physicsUnits import fb
from itertools import product
import unum

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

        # Production cross-section:
        self.prodXSec = None
        # Total branching ratio (product of BRs)
        self.decayBRs = 1.0
        # Maximum weight:
        self.maxWeight = None
        # Include additional attributes
        self._sorted = False
        # List of SMS topologies which could have generated self (such as during compression)
        self.ancestors = [self]
        self.smsID = 0  # SMS identifier
        # Type of analyses which have SMS matching self:
        self.coveredBy = set()
        # Type of analyses which have SMS matching self and for
        # which the physical parameters are covered
        self.testedBy = set()

    def __cmp__(self,other):
        """
        Compare self to other based on the compareTo method.

        :parameter other: TheorySMS object

        :return: 0, if objects are equal, -1 if self < other, 1 if other > sekf
        """

        return self.compareTo(other)

    def __lt__(self, other):
        return self.__cmp__(other) == -1

    def __gt__(self, other):
        return self.__cmp__(other) == 1

    def __eq__(self, other):
        return self.__cmp__(other) == 0

    def __ne__(self, other):
        return self.__cmp__(other) != 0

    def __hash__(self):
        return object.__hash__(self)

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
        newSMS.prodXSec = self.prodXSec + other.prodXSec
        newSMS.ancestors = self.ancestors[:] + other.ancestors[:]
        newSMS.weightList = self.weightList + other.weightList
        newSMS.maxWeight = self.maxWeight + other.maxWeight
        # Decay BRs can no longer be properly defined:
        newSMS.decayBRs = None

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
            self._canonName = self.computeCanonName()
        if sort:
            self.sort(force=True)
        if weight:
            self.weightList = self.computeWeightList()

    def compareTo(self, other):
        """
        Compare self to other.
        If the SMS are not sorted, sort them and then do a direct comparison of each
        node with the same nodeIndex.

        :param other: SMS object to be compared against self

        :return: 0, if objects are equal, -1 if self < other, 1 if other > sekf
        """

        if not isinstance(other,TheorySMS):
            raise SModelSError("Can not compare TheorySMS and %s" %str(type(other)))

        # Make sure the SMS are sorted
        # (if it has already been sorted, it does nothing)
        self.sort()
        other.sort()

        cmp = self.compareSubTrees(other,self.rootIndex,other.rootIndex)
        return cmp

    def copy(self, emptyNodes=False):
        """
        Returns a shallow copy of self.

        :param emptyNodes: If True, does not copy any of the nodes from self.

        :return: TheorySMS object
        """

        newSMS = TheorySMS()

        newSMS.maxWeight = self.maxWeight
        newSMS.prodXSec = self.prodXSec
        newSMS.decayBRs = self.decayBRs
        if hasattr(self,'weightList'):
            newSMS.weightList = self.weightList.copy()
        newSMS._sorted = self._sorted
        newSMS.ancestors = self.ancestors[:]
        newSMS.smsID = self.smsID
        newSMS.coveredBy = set(list(self.coveredBy)[:])
        newSMS.testedBy = set(list(self.testedBy)[:])

        # Set nodes:
        if not emptyNodes:
            nodesObjDict = {n : node for n,node in zip(self.nodeIndices,self.nodes)}
            newSMS.copyTreeFrom(self, nodesObjDict)

        return newSMS

    def addNodesFrom(self, other):
        """
        Combines the nodes (add particles) in equivalent nodes in each trees.
        Can only be done if the trees have the same topology and ordering.

        :param other: TheorySMS object
        """

        if other.canonName != self.canonName:
            raise SModelSError("Can not add trees with distinct topologies")

        newMapping = {}
        for n in self.nodeIndices:
            newMapping[n] = self.indexToNode(n) + other.indexToNode(n)

        # Just change the nodeIndex -> node mapping:
        self.updateNodeObjects(nodeObjectDict=newMapping)

    def attachDecay(self, motherIndex, decayNodes, br=1.0, copy=True):
        """
        Attaches a decay to self. If copy = True, returns a copy of self
        with the decay attached.

        :param motherIndex: Node index for the mother to which the decay should be added.
        :param decayNodes: Particle nodes for the daughters.
        :br: Branching ratio value for the decay
        :param copy: if True, return a copy of self, with the decay attached.

        :return: new tree with the other composed with self.
        """

        if self.out_degree(motherIndex) != 0:
            raise SModelSError("Can not attach decay to an intermediate node.")


        motherNode = decayNodes[0]
        daughterNodes = decayNodes[1]

        if copy:
            newSMS = self.copy()
        else:
            newSMS = self

        # Update maximum weight
        newSMS.maxWeight = self.maxWeight*br
        # Update BRs:
        newSMS.decayBRs = self.decayBRs*br

        # Update mother node:
        newSMS.updateNodeObjects({motherIndex : motherNode})
        # Add daughter nodes and edges:
        for d in daughterNodes:
            dIndex = newSMS.add_node(d)
            newSMS.add_edge(motherIndex, dIndex)

        # The tree is no longer sorted
        newSMS._sorted = False

        return newSMS

    def computeWeightList(self):
        """
        Computes the SMS weight (production cross-section*BRs) using maxWeight
        and the production cross-section.

        :return: CrossSectionList object
        """

        prodXSec = self.prodXSec
        brs = self.decayBRs

        return prodXSec*brs

    def getBSMattr(self,attr):
        """
        Returns the list of attribute values for the BSM particles
        following the node ordering in self.

        :param attr: Attribute label (str)

        :return: List of attribute values
        """

        values = []
        for node in self.nodes:
            if node.isSM:
                continue
            values.append(getattr(node,attr))

        return values

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
            self.add_edges_from(product([nodeIndex],sorted_daughters))

        # Finally, after sorting the subtrees,
        # make sure the nodes are sorted and numbered according
        # to the generations (breadth-first search)
        if nodeIndex == self.rootIndex:
            self.bfs_sort(numberNodes=True)
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
            cName = self.nodeCanonName(nodeIndex)
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
            if self.compareSubTrees(self, s[j+1], pivot) < 0:
                s[j+1], s[i+1] = s[i+1], s[j+1]
                i += 1
        s[0], s[i] = s[i], s[0]
        first_part = self.sortCommonSubTrees(s[:i])
        second_part = self.sortCommonSubTrees(s[i+1:])
        first_part.append(s[i])

        sortedList = first_part + second_part

        return sortedList

    def compareSubTrees(self, other, n1, n2):
        """
        Compare the subtrees generated by the nodeIndex n1 in self
        and the nodeIndex n2 in other.

        :param other: TheorySMS object (if other=self compare subtrees of the same SMS).
        :param n1: Node index for the root of subtree1
        :param n2: Node index for the root of subtree2

        :return: 0, if subtrees are equal, -1 if subtree1 < subtree2, 1 if subtree1 > subtree2
        """

        # Compare node canonical names
        cName1 = self.nodeCanonName(n1)
        cName2 = other.nodeCanonName(n2)
        if cName1 != cName2:
            if cName1 < cName2:
                return -1
            else:
                return 1

        # Compare nodes
        cmp = self.compareNodes(n1,n2)
        if cmp != 0:
            return cmp

        daughters1 = self.daughterIndices(n1)
        daughters2 = other.daughterIndices(n2)
        # If nodes are leaves, return 0
        if len(daughters1) == len(daughters2) == 0:
            return 0

        # Check if the daughters from n2 match the ones from n1:
        # (the daughters should be sorted at this point)
        for i1, d1 in enumerate(daughters1):
            d2 = daughters2[i1]
            cmp = self.compareSubTrees(other, d1, d2)
            if cmp != 0:
                return cmp
        return 0

    def compareNodes(self,nodeIndex1,nodeIndex2):
        """
        Convenience function for defining how nodes are compared
        within the SMS.

        :param nodeIndex1: Index of first node
        :param nodeIndex2: Index of second node

        :return: 1 if node1 > node2, -1 if node1 < node2, 0 if node1 == node2.
        """

        # Comparison parameters:
        node1 = self.indexToNode(nodeIndex1)
        node2 = self.indexToNode(nodeIndex2)

        # Directly use node comparison:
        cmp = node1.compareTo(node2)
        return cmp

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
            for jgen, smsList in ancestor._getAncestorsDict(igen + 1).items():
                if jgen not in ancestorsDict:
                    ancestorsDict[jgen] = []
                ancestorsDict[jgen] += smsList

        return ancestorsDict

    def getAncestors(self):
        """
        Get a list of all the ancestors of self.
        The list is ordered so the mothers appear first, then the grandmother,
        then the grandgrandmothers,...

        :return: A list of SMS objects containing all the ancestors sorted by generation.
        """

        ancestorsDict = self._getAncestorsDict()
        allAncestors = []

        for igen in sorted(ancestorsDict.keys()):
            allAncestors += ancestorsDict[igen]

        return allAncestors

    def isRelatedTo(self, other):
        """
        Checks if self has other as an ancestor or they share
        ancestors. Returns True if self and other have at least one ancestor in common,
                   otherwise returns False.

        :return: True/False
        """

        if self is other:
            return True

        # Get ids for all ancestors of self, including self
        ancestorsA = set([id(self)]+[id(sms) for sms in self.getAncestors()])
        # Get ids for all ancestors of other, including other
        ancestorsB = set([id(other)]+[id(sms) for sms in other.getAncestors()])

        # Check if any ids intersect
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

        # Start with self and make a copy later if needed
        newSMSList = [self]  # Dummy list to store newSMS

        # Loop over nodes from root to leaves:
        for momIndex, daughterIndices in self.genIndexIterator():
            newSMS = newSMSList[0]
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

            # If making first compression, copy self:
            if newSMS is self:
                newSMS = self.copy()
                newSMS.ancestors = [self]
                newSMSList[0] = newSMS

            # Get grandmother:
            gMomIndex = newSMS.parentIndex(momIndex)
            # Remove mother and all SM daughters:
            removeIndices = smDaughters + [momIndex]
            newSMS.remove_nodes_from(removeIndices)

            # Attach BSM daughter to grandmother:
            newSMS.add_edge(gMomIndex, bsmDaughter)

        newSMS = newSMSList[0]
        # If no compression was made, return None
        if newSMS is self:
            return None

        # Set the compressed topology weight as the original weight
        # (it can not longer be computed from its nodes)
        newSMS.weightList = self.weightList
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

        # Start with self and make a copy later if needed
        newSMSList = [self]  # Dummy list to store newSMS

        while True:
            newSMS = newSMSList[0]
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

                # If making first compression, copy self:
                if newSMS is self:
                    newSMS = self.copy()
                    newSMS.ancestors = [self]
                    newSMSList[0] = newSMS

                newSMS.remove_nodes_from(daughterIndices)
                # Replace mother particle by invisible (generic) particle
                # with its width equal to the maximum width amongst the daughters
                maxWidth = max([d.totalwidth for d in daughters])
                invParticle = Particle(label='inv', mass=mom.mass,
                                       eCharge=0, colordim=1,
                                       _isInvisible=True,
                                       totalwidth=maxWidth,
                                       pdg=mom.pdg, isSM=mom.isSM)
                newmom = mom.copy()
                newmom.particle = invParticle
                newSMS.updateNodeObjects({momIndex : newmom})

                # For safety break loop since tree structure changed
                break

            if newSMS.number_of_nodes() == nNodes:
                break  # No compressions could be made, stop

        newSMS = newSMSList[0]
        if newSMS is self:
            return None
        # Recompute the global properties (except for the weightList)
        # and sort the new SMS
        newSMS.setGlobalProperties(weight=False)

        return newSMS


