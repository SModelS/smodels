"""
.. module:: genericSMS
   :synopsis: This is a base class for describing Simplified Model Topologies  using a rooted tree syntax.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


from smodels.base.exceptions import SModelSBaseError as SModelSError
from smodels.base.inclusiveObjects import InclusiveValue
from collections import OrderedDict
from itertools import product

class GenericSMS(object):
    """
    A generic class for describing and manipulating 
    Simplified Model Topologies based on graphs concepts.
    """

    def __init__(self):
        """
        Initialize basic attributes.
        """

        self._canonName = None
        self._rootIndex = None
        self._successors = OrderedDict()  # Stores the nodes and their successors (daughters)
        self._predecessors = {}  # Stores the nodes and their predecessors (parents)
        self._nodesMapping = {}  # Stores the nodeIndex->node object mapping
        self._nodeCanonNames = {}  # Stores the canonical names for the nodes
        self._finalStates = {}  # Stores the final states of each node
        self._sorted = False # Tag SMS as sorted or not

    def __hash__(self):
        return object.__hash__(self)

    def __repr__(self):
        """
        Returns the string representation of the tree.
        """

        return str(self)

    def __str__(self):
        """
        Returns a string representing the process
        described by the tree.
        """

        return self.treeToString()

    def __getattr__(self, attr):
        """
        If the attribute has not been defined for self
        try to fetch it from its nodes

        :param attr: Attribute name
        :return: Attribute value
        """

        # If calling another special method, return default (required for pickling)
        if (attr.startswith('__') and attr.endswith('__')) or attr in dir(self):
            return self.__getattribute__(attr)


        try:
            val = [getattr(node, attr) if node is not self.root else None
                   for node in self.nodes]
            return val
        except AttributeError:
            raise AttributeError(f"Neither SMS nor nodes have attribute ``{attr}''")

    def add_node(self, node, nodeIndex=None):
        """
        Adds a node object to the tree. If nodeIndex is None,
        the node index will be automatically assigned.

        :param node: ParticleNode object
        :param nodeIndex: The index for the ParticleNode. It must not
                          match any other indices already in the tree.

        :return: The node index for the newly added node
        """

        if nodeIndex is None:   
            if not self._successors:
                nodeIndex = 0
            else:
                nodeIndex = max(self.nodeIndices)+1
        elif nodeIndex in self._successors:
            raise SModelSError("Trying to add a node with a nodeIndex already in the tree.")
        self._successors[nodeIndex] = []
        self._nodesMapping[nodeIndex] = node

        return nodeIndex

    def add_nodes_from(self, nodes):
        """
        Adds a list of nodes to the Tree.

        :param nodes: List of ParticleNode objects

        :return: A list of node indices for the newly added nodes
        """

        nodeIndices = []
        for node in nodes:
            nodeIndices.append(self.add_node(node))

        return nodeIndices

    def remove_node(self, nodeIndex):
        """
        Removes a node from the tree if the nodeIndex is in the tree.
        The node is removed as well as its appearence in any edges.

        :param nodeIndex: Node index
        """

        if nodeIndex in self._successors:
            self._successors.pop(nodeIndex)
            self._nodesMapping.pop(nodeIndex)

        if nodeIndex in self._nodeCanonNames:
            self._nodeCanonNames.pop(nodeIndex)

        for nodeA, daughtersA in self._successors.items():
            if nodeIndex not in daughtersA:
                continue
            daughtersA = [d for d in daughtersA[:] if d != nodeIndex]
            self._successors[nodeA] = daughtersA[:]

        if nodeIndex in self._predecessors:
            self._predecessors.pop(nodeIndex)

        for nodeA, momA in list(self._predecessors.items()):
            if momA == nodeIndex:
                self._predecessors.pop(nodeA)

        if nodeIndex in self._finalStates:
            self._finalStates.pop(nodeIndex)

    def remove_nodes_from(self, nodeIndices):
        """
        Removes a list of nodes from the Tree.

        :param nodeIndices: List of node indices
        """

        for nodeIndex in nodeIndices:
            self.remove_node(nodeIndex)

    def add_edge(self, nodeIndexA, nodeIndexB):
        """
        Adds a directed edge to existing nodes in the Tree (nodeA -> nodeB).

        :param nodeIndexA: Index for node A
        :param nodeIndexB: Index for node B
        """

        self._successors[nodeIndexA].append(nodeIndexB)
        self._predecessors[nodeIndexB] = nodeIndexA

    def add_edges_from(self, edges):
        """
        Adds a list of directed edges to the Tree.

        :param edges: List of tuples containing node indices
                      (e.g. [(nodeIndexA,nodeIndexB),(nodeIndexA,nodeIndexC),...])
        """

        for edge in edges:
            self.add_edge(edge[0], edge[1])

    def remove_edge(self, nodeIndexA, nodeIndexB):
        """
        Removes an edge from the tree if the edge
        (nodeIndexA -> nodeIndexB) is in the tree.

        :param nodeIndexA: Index for node A
        :param nodeIndexB: Index for node B
        """

        if nodeIndexA in self._successors:
            daughters = self._successors[nodeIndexA]
            daughters = [d for d in daughters if d != nodeIndexB]
            self._successors[nodeIndexA] = daughters

        if nodeIndexB in self._predecessors:
            if self._predecessors[nodeIndexB] == nodeIndexA:
                self._predecessors.pop(nodeIndexB)

    def remove_edges(self, edges):
        """
        Removes edges from the tree if they appear in the tree.

        :param edges: List of tuples containing node indices
                      (e.g. [(nodeIndexA,nodeIndexB),(nodeIndexA,nodeIndexC),...])
        """
        for edge in edges:
            self.remove_edge(edge[0], edge[1])

    def clear(self):
        """
        Remove all nodes and edges from the graph, but
        keep its canonName.
        """

        self._successors = OrderedDict()
        self._predecessors = {}
        self._nodesMapping = {}
        self._nodeCanonNames = {}
        self._finalStates = {}
        self._rootIndex = None

    def indexToNode(self, nodeIndex):
        """
        Returns the node object with index nodeIndex.
        If nodeIndex is a list of indices, return the corresponding
        list of node objects.

        :param nodeIndex: Integer or list of integers of
                          node indices.

        :return: Node object or list of Node objects
        """

        if isinstance(nodeIndex,int):
            return self._nodesMapping[nodeIndex]
        elif isinstance(nodeIndex,list):
            return [self._nodesMapping[n] for n in nodeIndex]
        elif isinstance(nodeIndex,tuple):
            return tuple([self._nodesMapping[n] for n in nodeIndex])
        else:
            raise SModelSError(f"Can not convert object of type {str(type(nodeIndex))} to nodes")

    def daughterIndices(self, nodeIndex, ignoreInclusiveNodes=False):
        """
        Returns the list of node indices corresponding to the
        daughters of nodeIndex.

        :param nodeIndex: Parent node index
        :param ignoreInclusiveNodes: If True, skips inclusive nodes
        """

        daughters = self._successors[nodeIndex]
        if ignoreInclusiveNodes:
            daughters = [d for d in daughters[:]
                         if not self.indexToNode(d).isInclusive]

        return daughters

    def daughters(self, nodeIndex,  ignoreInclusiveNodes=False):
        """
        Returns the list of node objects corresponding to the
        daughters of nodeIndex.

        :param nodeIndex: Parent node index
        :param ignoreInclusiveNodes: If True, it skips inclusive nodes
        """

        daughterIndices = self.daughterIndices(nodeIndex, ignoreInclusiveNodes)
        daughters = self.indexToNode(daughterIndices)

        return daughters

    def parentIndex(self, nodeIndex):
        """
        Returns the node index corresponding to the
        parent of nodeIndex.

        :param nodeIndex: Daughter node index
        """

        return self._predecessors[nodeIndex]

    def parent(self, nodeIndex):
        """
        Returns the node object corresponding to the
        parent of nodeIndex.

        :param nodeIndex: Daughter node index
        """

        parentIndex = self.parentIndex(nodeIndex)
        parent = self.indexToNode(parentIndex)

        return parent

    @property
    def rootIndex(self):
        """
        Returns the index of the root node (primary vertex) of the tree.
        If it has not been defined, compute it.

        :return: root node index
        """

        if self._rootIndex is None:
            root = [nodeIndex for nodeIndex in self.nodeIndices
                    if self.in_degree(nodeIndex) == 0]
            if len(root) != 1:
                raise SModelSError(f"Malformed Tree, {len(root)} root(s) have been found.")
            self._rootIndex = root[0]

        return self._rootIndex

    @property
    def root(self):
        """
        Returns the root node (primary vertex) of the tree.
        If it has not been defined, compute it.

        :return: root node
        """

        rootIndex = self.rootIndex
        root = self.indexToNode(rootIndex)

        return root

    @property
    def nodeIndices(self):
        """
        Returns the tist of node indices in the Tree.

        :return: List of indices (int)
        """


        nodeIndexList = list(self._successors.keys())

        return nodeIndexList

    @property
    def edgeIndices(self):
        """
        Returns the list of edges indices (pairs of integers) in the Tree.

        :return: List of edge indices
        """

        edgesList = []
        for n in self.nodeIndices:
            edgesList += list(product([n],self.daughterIndices(n)))

        return edgesList

    @property
    def nodes(self):
        """
        Returns the tist of ParticleNode objects in the Tree.

        :return: List of ParticleNode objects
        """


        nodeList = self.indexToNode(self.nodeIndices)

        return nodeList

    @property
    def edges(self):
        """
        Returns the list of edges (pairs of node objects) in the Tree.

        :return: List of edges
        """

        edgesList = [tuple(self.indexToNode(edgeTuple)) for edgeTuple in self.edgeIndices]

        return edgesList

    @property
    def canonName(self):
        """
        Returns the canonName. If not defined, it will be computed.

        :return: Canonical name (int)
        """

        if self.rootIndex not in self._nodeCanonNames:
            self.computeCanonName()

        return self._nodeCanonNames[self.rootIndex]

    def computeCanonName(self, nodeIndex=None):
        """
        Recursively sets the canonName for each node.
        Returns the canonical name in integer form.

        :param nodeIndex: Node index to set the name for. If None, it will use the root

        :return: Integer representing the Tree canonical name
        """

        if not self.number_of_nodes():
            return None

        if nodeIndex is None:
            nodeIndex = self.rootIndex

        # Set the final state
        node = self.indexToNode(nodeIndex)
        # If it is inclusive node set its name to an inclusive integer
        # and return its name (no need to check the children)
        if node.isInclusive or node.inclusiveList:
            canonName = InclusiveValue()
            self._nodeCanonNames[nodeIndex] = canonName
            return canonName

        children = self.daughterIndices(nodeIndex)
        if not children:
            canonName = 10
        else:
            tp = [self.computeCanonName(n) for n in children]
            if any(isinstance(name, InclusiveValue) for name in tp):
                canonName = InclusiveValue()
            else:
                tp = sorted(tp)
                tpStr = '1' + "".join(str(c) for c in tp) + '0'
                canonName = int(tpStr)

        self._nodeCanonNames[nodeIndex] = canonName

        return canonName

    def nodeCanonName(self,nodeIndex):
        """
        Returns the canon name for the node.

        :param nodeIndex: Index of the node

        :return: Canonical name (int)
        """

        if nodeIndex not in self._nodeCanonNames:
            self.computeCanonName(nodeIndex)

        return self._nodeCanonNames[nodeIndex]

    def out_degree(self, nodeIndex):
        """
        Computes the number of outgoing edges from the node
        (number of daughters).

        :param nodeIndex: Node index (int)

        :return: Number of outgoing edges (int)
        """

        if nodeIndex not in self.nodeIndices:
            return 0
        else:
            return len(self.daughterIndices(nodeIndex))

    def in_degree(self, nodeIndex):
        """
        Computes the number of incoming edges to the node
        (number of parents).

        :param nodeIndex: Node index (int)

        :return: Number of incoming edges (1 or 0)
        """

        if nodeIndex in self._predecessors:
            if self._predecessors[nodeIndex] is not None:
                return 1

        return 0

    def number_of_nodes(self):
        """
        Returns the total number of nodes in the Tree.


        :return: Number of nodes (int)
        """

        return len(self.nodeIndices)

    def genIndexIterator(self, nodeIndex=None,
                         includeLeaves=False, ignoreInclusiveNodes=False):
        """
        Returns an iterator over the generations (mother and its daughters)
        of node indices starting at nodeIndex using a breadth first search.

        :param nodeIndex: Node index from tree. If None, starts at tree root.
        :param includeLeaves: If True, it will consider the leaves (undecayed nodes)
                              as moms in the iterator (with an empty daughters list)
        :param ignoreInclusiveNodes: If True, it skip inclusive nodes and its descendents.

        :return: Iterator over nodes.
        """

        if nodeIndex is None:
            nodeIndex = self.rootIndex

        mom = nodeIndex
        if ignoreInclusiveNodes:
            if self.indexToNode(mom).isInclusive:
                return []

        daughters = self.daughterIndices(mom,ignoreInclusiveNodes)

        generation = [(mom, daughters)]
        while generation:
            for pair in generation:
                yield pair
            next_generation = []
            for pair in generation:
                mom, daughters = pair
                for new_mom in daughters:
                    if new_mom not in self.nodeIndices:
                        continue
                    new_daughters = self.daughterIndices(new_mom,ignoreInclusiveNodes)

                    if new_daughters or includeLeaves:
                        next_generation.append((new_mom, new_daughters))
            generation = next_generation

    def dfsIndexIterator(self, nodeIndex=None, ignoreInclusiveNodes=False):
        """
        Iterates over the node indices following a depth-first traversal of the tree
        starting at nodeIndex. If nodeIndex is None, include all nodes.

        :param nodeIndex: Node index to which start the iterator
                          (the corresponding node is NOT included in the iterator)
        :param ignoreInclusiveNodes: If True, it skip inclusive nodes and its descendents.

        :return: Iterator over node indices
        """

        if nodeIndex is None:
            nodeIndex = self.rootIndex
            yield self.rootIndex

        nodes = self.daughterIndices(nodeIndex,ignoreInclusiveNodes)

        visited = set()   # Store visited nodes
        depth_limit = len(nodes)
        # Loop over nodes
        for node in nodes:
            if node in visited:
                continue
            yield node
            visited.add(node)
            daughters = [iter(self.daughterIndices(node,
                                                   ignoreInclusiveNodes))]
            while daughters:
                last_children = daughters[-1]  # Get the last added children
                try:
                    child = next(last_children)
                    if child not in visited:
                        yield child
                        visited.add(child)
                        # Add daughters of child
                        daughters.append(iter(self.daughterIndices(child,
                                                                   ignoreInclusiveNodes)))
                except StopIteration:
                    # All the children have been visited, remove from list
                    daughters.pop()

    def bfs_sort(self, numberNodes=False):
        """
        Sort the nodes according to their appearence in a
        breadth first search.

        :param numberNodes: If True, renumber the nodes according to their bfs order

        :return: Dictionary with old node indices as keys
                 and new indices as values.
        """

        orderedList = []
        for nodeIndex, _ in self.genIndexIterator(includeLeaves=True):
            orderedList.append(nodeIndex)

        # Sort according to the bfs order
        self.sortAccordingTo(orderedList)

        # Re-number nodes according to their order
        if numberNodes:
            indexDict = {n : orderedList.index(n) for n in self.nodeIndices}
            self.relabelNodeIndices(nodeIndexDict=indexDict)
        else:
            # Define dummy dict
            indexDict = {n : n for n in self.nodeIndices}

        return indexDict

    def sortAccordingTo(self,indicesList):
        """
        Sort the nodes according to their order in indicesList.

        :param indicesList: List of node indices used to sort the nodes.
        """

        newSuccessors = OrderedDict()
        indices = self.nodeIndices
        sortList = indicesList[:]
        # Indices not present in indicesList are left
        # put at the end of the list
        for nodeIndex in indices:
            if nodeIndex not in sortList:
                sortList.append(nodeIndex)

        sortedIndices = sorted(indices, key = lambda n: sortList.index(n))
        # Go over sorted indices and create new successors dict
        for nodeIndex in sortedIndices:
            daughters = self.daughterIndices(nodeIndex)
            # Sort daughters according to the list
            sortedDaughters = sorted(daughters, key = lambda n: sortList.index(n))
            newSuccessors[nodeIndex] = sortedDaughters[:]
        self._successors = newSuccessors

    def sort(self, nodeIndex=None, force=False):
        """
        Sort subtree of self generated by nodeIndex.
        If nodeIndex is None, sort the tree and re-number the nodes
        according to the bfs order (after sorting).
        If the self is already tagged as sorted and force = False,
        do nothing.

        :param nodeIndex: Node index
        :param force: If True, will sort even if self is tagged as sorted.

        :return: Dictionary with old node indices as keys
                 and new indices as values.
        """

        # If tree is already sorted, do nothing
        if hasattr(self,'_sorted') and self._sorted and not force:
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
            nodeIndexMap = self.bfs_sort(numberNodes=True)
            # Tag the tree as sorted
            self._sorted = True
            return nodeIndexMap

    def sortSubTrees(self, subtreeList):
        """
        Sorts a list of subtrees of self generated by the nodes
        in subtreeList.

        :param subtreeList: List of node indices to be considered as roots
                            of the subtrees.

        :return: Sorted list of node indices.
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

        :param subtreeList: List of node indices

        :return: Sorted list of node indices.
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
        cmp = self.compareNodes(other,n1,n2)
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

    def compareNodes(self,other,nodeIndex1,nodeIndex2):
        """
        Convenience function for defining how nodes are compared
        within the SMS.

        :param other: TheorySMS object (if other=self compare subtrees of the same SMS).
        :param nodeIndex1: Index of first node
        :param nodeIndex2: Index of second node

        :return: 1 if node1 > node2, -1 if node1 < node2, 0 if node1 == node2.
        """

        # Comparison parameters:
        node1 = self.indexToNode(nodeIndex1)
        node2 = other.indexToNode(nodeIndex2)

        # Directly use node comparison:
        cmp = node1.compareTo(node2)
        return cmp

    def copyTreeFrom(self, other, nodesObjDict):
        """
        Replaces the tree structure (nodes, edges, indices,...)
        by the structure in other. Uses the nodesObjDict to set the
        new node mapping (nodeIndex > nodeObj).

        :param other: SMS object
        :param nodesObjDict: Dictionary where keys are node indices (from other)
                             and values are node objects.
        """

        self._successors = {n : dList[:] for n,dList in other._successors.items()}
        self._predecessors = {d : parent for d,parent in other._predecessors.items()}
        self._canonName = other._canonName
        self._rootIndex = other._rootIndex
        self._nodesMapping = {nodeIndex : n for nodeIndex,n in nodesObjDict.items()}
        self._nodeCanonNames = {nodeIndex : cName for nodeIndex,cName
                                in other._nodeCanonNames.items()}
        self._finalStates = {nodeIndex : pList[:] for nodeIndex,pList
                             in self._finalStates.items()}

    def treeToString(self, 
                     removeIndicesFrom='stable'):
        """
        Convert the tree to a process string (e.g. '(PV(0) > gluino(1),squark(2)), (gluino(1) >
                           MET(3),jet(4),jet(5)), (squark(2) > HSCP(6),u(7))')
        
        Node indices can be removed from specific particles using the removeIndicesFrom option. Allowed values are: 
        None -> keep all indices
        'SM' -> remove indices from SM particles
        'stable' -> remove indices from stable (undecayed) particles
        'all' -> remove indices from all particles
        
        The default is to remove indices from stable particles.
                           
        :param removeIndicesFrom: If defined, will remove indices from particles according to their properties.

        :return: String describing the process
        """

        if removeIndicesFrom not in [None, 'SM', 'stable', 'all']:
            raise SModelSError(f"removeIndicesFrom = {removeIndicesFrom} value not accepted for treeToString")
        rmFrom = removeIndicesFrom
        
        smsStr = ""
        rootIndex = self.rootIndex
        for momIndex, daughterIndices in self.genIndexIterator(rootIndex):

            # Convert from indices to node objects
            mom = self.indexToNode(momIndex)
            daughters = self.indexToNode(daughterIndices)
            if momIndex == rootIndex: # Always remove from PV
                smsStr += f'({mom} > '
            elif rmFrom is not None:
                if rmFrom == 'all':
                    smsStr += f'({mom} > '
                elif rmFrom == 'SM' and mom.isSM:
                    smsStr += f'({mom} > '
                elif rmFrom == 'stable' and self.out_degree(momIndex) == 0:
                    smsStr += f'({mom} > '
                else: # Add index
                    smsStr += '(%s(%i) > ' % (mom, momIndex)
            else: # If None, always add index
                smsStr += '(%s(%i) > ' % (mom, momIndex)

            # Add daughters
            for iD, d in enumerate(daughters):
                dIndex = daughterIndices[iD]
                if rmFrom is not None:
                    if rmFrom == 'all':
                        smsStr += f'{d},'
                    elif rmFrom == 'SM' and d.isSM:
                        smsStr += f'{d},'
                    elif rmFrom == 'stable' and self.out_degree(dIndex) == 0:
                        smsStr += f'{d},'
                    else: # Add index
                        smsStr += '%s(%i),' % (d, dIndex)
                else: # If None, always add index
                    smsStr += '%s(%i),' % (d, dIndex)
            smsStr = smsStr[:-1] + '), '
        smsStr = smsStr[:-2]

        return smsStr

    def treeToBrackets(self):
        """
        Convert the Tree to a nested list with the Z2 even
        final states. The Z2 odd final states (e.g. 'MET', 'HSCP') are
        not included. The Tree must be Z2-preserving and represent the
        pair production cascade decay of two Z2-odd particles.

        :return: Nested list with the strings for the Z2-even final states
                 (e.g. [[['e-','mu'],['L']],[['jet']]])
        """

        branches = self.daughters(self.rootIndex)
        branchIndices = self.daughterIndices(self.rootIndex)
        finalState = []
        intermediateState = []
        branchList = []
        if len(branches) != 2:
            raise SModelSError(f"Can not convert tree to bracket with {len(branches)} branches")
        for ib, b in enumerate(branches):
            bIndex = branchIndices[ib]
            intermediateState.append([])
            branchList.append([])
            # Deal separately with the case where the primary mother is stable:
            if self.out_degree(bIndex) == 0:
                if not b.isSM:
                    finalState.append(str(b))
                    continue
                else:
                    raise SModelSError("Can not convert tree with Z2-violating decays to bracket")
            for momIndex, daughterIndices in self.genIndexIterator(bIndex):
                # Convert from indices to node objects
                mom = self.indexToNode(momIndex)
                daughters = self.indexToNode(daughterIndices)

                vertexList = [str(d) for d in daughters if d.isSM]
                fstates = []
                for idaughter, daughter in enumerate(daughters):
                    if daughter.isSM:
                        continue
                    if self.out_degree(daughterIndices[idaughter]) != 0:
                        continue
                    fstates.append(str(daughter))

                if daughters:
                    if len(vertexList) != len(daughters) - 1 or len(fstates) > 1:
                        raise SModelSError(f"Can not convert tree with Z2-violating decays to bracket: \n  {self.treeToString()}")
                intermediateState[ib].append(str(mom))
                branchList[ib].append(vertexList)
                finalState += fstates

        return branchList, finalState, intermediateState

    def getFinalStates(self, nodeIndex=None):
        """
        Get the list of nodes which have not decayed (appear at the top of the tree).
        If source is defined, get the final states generated by the cascade decay of the source
        node. It also caches the finalState to self._finalStates.

        :param nodeIndex: Node index for which to get the final states.

        :returns: list of node indices
        """

        if nodeIndex is None:
            nodeIndex = self.rootIndex

        if nodeIndex in self._finalStates:
            return self._finalStates[nodeIndex]

        # For leaves, the final state is themselves:
        if self.out_degree(nodeIndex) == 0:
            self._finalStates[nodeIndex] = [nodeIndex]
        else:
            finalStates = []
            for d in self.daughterIndices(nodeIndex):
                finalStates += self.getFinalStates(d)
            self._finalStates[nodeIndex] = finalStates

        return self._finalStates[nodeIndex]

    def relabelNodeIndices(self,nodeIndexDict):
        """
        Relabel node indices according to nodeIndexDict.
        For node indices not appearing in nodeIndexDict nothing is done.

        :param nodeIndexDict: Dictionary with current node indices as keys
                              and new indices as values
        """

        if any(nodeIndex not in nodeIndexDict for nodeIndex in self.nodeIndices):
            raise SModelSError("Dictionary for relabelling nodes must contain all node indices")

        newMapping = {}
        newSuccessors = OrderedDict()
        newPredecessors = {}
        newCanonNames = {}
        newFinalStates = {}
        for oldIndex, newIndex in nodeIndexDict.items():
            # Update daughter indices
            newDaughters = []
            for d in self.daughterIndices(oldIndex):
                if d in nodeIndexDict:
                    newDaughters.append(nodeIndexDict[d])
                else:
                    newDaughters.append(d)

            newPredecessors.update({d : newIndex
                                    for d in newDaughters})
            # Add entry to newSuccessors dict:
            newSuccessors[newIndex] = newDaughters
            # Add entry to newMapping dict
            newMapping[newIndex] = self.indexToNode(oldIndex)
            # Add entry to canonNames dict
            if oldIndex in self._nodeCanonNames:
                newCanonNames[newIndex] = self._nodeCanonNames[oldIndex]
            if oldIndex in self._finalStates:
                newFinalStates[newIndex] = self._finalStates[oldIndex]

        # Update dicts:
        self._successors = newSuccessors
        self._predecessors = newPredecessors
        self._nodesMapping = newMapping
        self._nodeCanonNames = newCanonNames
        self._finalStates = newFinalStates
        self._rootIndex = nodeIndexDict[self.rootIndex]

    def updateNodeObjects(self, nodeObjectDict):
        """
        Update the node index -> node object mapping.
        Only affects the indices appearing in nodeObjectDict.

        :param nodeObjectDict: Dictionary with current node indices as keys
                              and new node objects as values
        """

        for nodeIndex, newObj in nodeObjectDict.items():
            self._nodesMapping[nodeIndex] = newObj

    def checkConsistency(self):
        """
        Make sure the tree has the correct topology(directed rooted tree).
        Raises an error otherwise.
        """

        malformedTree = False
        if len(self.nodeIndices) > 1:
            rootIndex = self.rootIndex

            # Check if root has no parents and at least one daughter
            if self.in_degree(rootIndex) != 0 or self.out_degree(rootIndex) == 0:
                malformedTree = True

            nNodes = len(self.nodeIndices)
            # Check if all nodes (except root) have a unique parent
            if any(self.in_degree(nodeIndex) != 1 for nodeIndex in self.nodeIndices
                   if nodeIndex != rootIndex):
                malformedTree = True

            # Check if all nodes can be reached from the root node
            if len(list(self.dfsIndexIterator())) != nNodes:
                malformedTree = True

        if malformedTree:
            raise SModelSError("Graph created with malformed structure (not  a tree).")

    def switchBranches(self):
        """
        If the SMS has a two branch structure (PV > X,Y), return
        a new SMS with its branches switched (PV > Y,X).
        Otherwise return None.

        :return: A new SMS object with the branches switched or None.
        """
    
        if len(self.daughterIndices(self.rootIndex)) != 2:
            return None

        smsNew = self.copy()
        branchIndices = smsNew.daughterIndices(smsNew.rootIndex)
        nodes = {}
        edges = {}
        for bIndex in branchIndices:
            nodes[bIndex] = smsNew.indexToNode(bIndex)
            edges[bIndex] = [(bIndex,d) for d in smsNew.daughterIndices(bIndex)]
            edges[bIndex].append((smsNew.rootIndex,bIndex))
            smsNew.remove_node(bIndex)
        
        for bIndex in branchIndices[::-1]:
            smsNew.add_node(nodes[bIndex],bIndex)
            smsNew.add_edges_from(edges[bIndex])
        smsNew.bfs_sort(numberNodes=True)

        return smsNew

    def draw(self, particleColor='steelblue2',
                smColor='lightpink2',
                pvColor='darkgray',
                labelAttr='label',
                attrUnit=None, filename=None, view=True,
                maxLabelSize=10,
                usePVimage=False,
                graph_kwargs={'layout' : 'dot', 'ranksep' : '0.3', 'rankdir' : "LR"},
                nodes_kwargs={'style' : 'filled', 'fontsize' : '10', 'color' : 'black','shape' : 'circle','margin' : '0'},
                edges_kwargs={'arrowhead' : 'vee', 'arrowsize' : '0.7',  'color' : 'grey53'}):
        """
        Draws Tree using matplotlib.

        :param particleColor: color for particle nodes
        :param smColor: color used for particles which have the isSM attribute set to True
        :param pvColor: color for primary vertex
        :param fontsize: Font size for labels
        :param labelAttr: attribute to be used as label. If None, will use the string representation of the node object. It can also be a dictionary with node indices as keys and the label strings as values.
        :param attrUnit: Unum object with the unit to be removed from label attribute(if applicable)
        :param filename: Filename to save drawing to.
        :param view: open a viewer after plotting
        :param maxLabelSize: Maximum size for the label string for the node. If the label is larger, it will be truncated.
                             If None/False/0, it will keep the full label.
        :param usePVimage: Path to a image file (png, bmp or jpeg) to be used instead of the primary vertex (PV) node.
        :param graph_kwargs: Dictionary with graph attributes to be used.
        :param nodes_kwargs: Dictionary with nodes attributes to be used.
        :param edges_kwargs: Dictionary with nodes attributes to be used.

        :return: Display a GraphViz Digraph object, if view is true (and save it to file if filename is defined)
        """

        try:
            import graphviz
        except ImportError:
            raise SModelSError("For drawing SMS objects, please install graphviz")

        nodesAndIndices = zip(self.nodes,self.nodeIndices)

        if labelAttr is None:
            labels = {nodeIndex: "" for _,nodeIndex in nodesAndIndices}
        elif isinstance(labelAttr,dict):
            labels = {k : v for k,v in labelAttr.items()}
        elif labelAttr == 'label':
            labels = {nodeIndex: str(n) for n,nodeIndex in nodesAndIndices}
        elif attrUnit is not None:
            labels = {nodeIndex: str(getattr(n, labelAttr).asNumber(attrUnit))
                      if (hasattr(n, labelAttr) and getattr(n, labelAttr) is not None)
                      else str(n) for n,nodeIndex in nodesAndIndices}
        elif labelAttr == 'node':
            labels = {nodeIndex: str(nodeIndex) for _,nodeIndex in nodesAndIndices}
        elif labelAttr == 'canonName':
            labels = {nodeIndex: str(self.nodeCanonName(nodeIndex))
                        for _,nodeIndex in nodesAndIndices}
        else:
            labels = {nodeIndex: str(getattr(n, labelAttr)) if hasattr(n, labelAttr)
                      else str(n) for n,nodeIndex in nodesAndIndices}

        for key in labels:
            if labels[key] == 'anyOdd':
                labels[key] = 'BSM'

        node_color = {}
        for n in self.nodeIndices:
            node = self.indexToNode(n)
            if n == self.rootIndex:
                node_color[n] = pvColor
            elif hasattr(node, 'isSM') and node.isSM:
                node_color[n] = smColor
            else:
                node_color[n] = particleColor

        # Truncate labels if needed:
        if maxLabelSize:
            for key,val in labels.items():
                if len(val) > maxLabelSize:
                    labels[key] = val[:maxLabelSize]+'...'

        dot = graphviz.Digraph()
        for key,val in graph_kwargs.items():
            dot.attr(**{key : str(val)})
        for nodeIndex in self.nodeIndices:            
            if labels[nodeIndex] == 'PV' and usePVimage:
                dot.node(str(nodeIndex),shape='plaintext', 
                         fontsize='16', label="",color='white',
                         width='0.6', height='1.2', fixedsize='true',
                         shapefile=usePVimage)
            else:
                nodeAttrs = {k : v for k,v in nodes_kwargs.items()}
                if 'label' not in nodeAttrs:
                    nodeAttrs['label'] = labels[nodeIndex]
                if 'fillcolor' not in nodeAttrs:
                    nodeAttrs['fillcolor'] = node_color[nodeIndex]
                dot.node(str(nodeIndex), **nodeAttrs)
        for edgeIndex in self.edgeIndices:
            dot.edge(str(edgeIndex[0]),str(edgeIndex[1]),
                     **edges_kwargs)

        # If filename is defined, save image
        if filename is not None:
            import os
            filename = os.path.abspath(filename)
            # dot.format = extension[1:]
            dot.render(outfile=filename, view=view, cleanup=True)

        # Try to display (for various circumstances)
        if view:
            try:
                display(dot) # for notebooks
            except NameError:
                try:
                    import os
                    fname = filename
                    if fname != None:
                        fname, _ = os.path.splitext(filename)
                    dot.view(filename=fname) # for terminals
                except (RuntimeError, graphviz.ExecutableNotFound,\
                        graphviz.CalledProcessError) as e:
                    pass
