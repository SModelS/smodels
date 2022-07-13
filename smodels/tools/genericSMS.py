"""
.. module:: tree
   :synopsis: This is a base class for describing Simplified Model Topologies
              using a rooted tree syntax.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.inclusiveObjects import InclusiveValue
from collections import OrderedDict
from itertools import product


class GenericSMS(object):
    """
    A generic class for describing and manipulating Simplified Model Topologies based
    on graphs concepts.
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
            raise AttributeError("Neither SMS nor nodes have attribute ``%s''" % attr)

    def add_node(self, node):
        """
        Adds a node object to the tree. The node index will be
        automatically assigned.

        :param node: ParticleNode object

        :return: The node index for the newly added node
        """

        if not self._successors:
            nodeIndex = 0
        else:
            nodeIndex = max(self.nodeIndices)+1
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
        elif isinstance(nodeIndex,(list,tuple)):
            return [self._nodesMapping[n] for n in nodeIndex]
        else:
            raise SModelSError("Can not convert object of type %s to nodes" %str(type(nodeIndex)))

    def daughterIndices(self, nodeIndex):
        """
        Returns the list of node indices corresponding to the
        daughters of nodeIndex.

        :param nodeIndex: Parent node index
        """

        return self._successors[nodeIndex]

    def daughters(self, nodeIndex):
        """
        Returns the list of node objects corresponding to the
        daughters of nodeIndex.

        :param nodeIndex: Parent node index
        """

        daughterIndices = self.daughterIndices(nodeIndex)
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
                raise SModelSError("Malformed Tree, %i root(s) have been found." % len(root))
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
        if node.isInclusive:
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

    def genIndexIterator(self, nodeIndex=None, includeLeaves=False):
        """
        Returns an iterator over the generations (mother and its daughters)
        of node indices starting at nodeIndex using a breadth first search.

        :param nodeIndex: Node index from tree. If None, starts at tree root.
        :param includeLeaves: If True, it will consider the leaves (undecayed nodes)
                              as moms in the iterator (with an empty daughters list)

        :return: Iterator over nodes.
        """

        if nodeIndex is None:
            nodeIndex = self.rootIndex

        mom = nodeIndex
        daughters = self.daughterIndices(mom)
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
                    new_daughters = self.daughterIndices(new_mom)
                    if new_daughters or includeLeaves:
                        next_generation.append((new_mom, new_daughters))
            generation = next_generation

    def dfsIndexIterator(self):
        """
        Iterates over the nodes following a depth-first traversal of the tree.

        :return: Iterator over nodes
        """

        nodes = self.nodeIndices
        yield nodes[0]
        visited = set()   # Store visited nodes
        depth_limit = len(nodes)
        # Loop over nodes
        for node in nodes:
            if node in visited:
                continue
            visited.add(node)
            daughters = [iter(self.daughterIndices(node))]
            while daughters:
                last_children = daughters[-1]  # Get the last added children
                try:
                    child = next(last_children)
                    if child not in visited:
                        yield child
                        visited.add(child)
                        # Add daughters of child
                        daughters.append(iter(self.daughterIndices(child)))
                except StopIteration:
                    # All the children have been visited, remove from list
                    daughters.pop()

    def bfs_sort(self, numberNodes=False):
        """
        Sort the nodes according to their appearence in a
        breadth first search.

        :param numberNodes: If True, renumber the nodes according to their bfs order
        """

        new_successors = OrderedDict()
        indexDict = {}
        newIndex = 0
        for mom, daughters in self.genIndexIterator(includeLeaves=True):
            indexDict[mom] = newIndex
            newIndex += 1
            new_successors[mom] = daughters

        self._successors = new_successors
        if numberNodes:
            self.relabelNodeIndices(nodeIndexDict=indexDict)

    def treeToString(self):
        """
        Convert the tree to a process string (e.g. '(PV > gluino(1),squark(2)), (gluino(1) >
                           MET,jet,jet), (squark(2) > HSCP,u)')

        :return: String describing the process
        """

        elStr = ""
        rootIndex = self.rootIndex
        nodesDict = {rootIndex: 0}
        counter = 1
        for momIndex, daughterIndices in self.genIndexIterator(rootIndex):

            # Convert from indices to node objects
            mom = self.indexToNode(momIndex)
            daughters = self.indexToNode(daughterIndices)
            # Add mom (if mom = 0 does not include index)
            if momIndex == 0:
                elStr += '(%s > ' % mom
            else:
                if momIndex not in nodesDict:
                    nodesDict[momIndex] = counter
                    counter += 1
                elStr += '(%s(%i) > ' % (mom, nodesDict[momIndex])

            # Add daughters for the decay, if daughter is unstable, add index
            for iD, d in enumerate(daughters):
                dIndex = daughterIndices[iD]
                if d.isInclusive and d.finalStates:
                    stable = False
                elif self.out_degree(dIndex) == 0:
                    stable = True
                else:
                    stable = False
                if stable:
                    elStr += '%s,' % d  # stable
                else:
                    if dIndex not in nodesDict:
                        nodesDict[dIndex] = counter
                        counter += 1
                    elStr += '%s(%i),' % (d, nodesDict[dIndex])
            elStr = elStr[:-1] + '), '
        elStr = elStr[:-2]

        # Deal separately with daughters for inclusive nodes (if it exists)
        for nodeIndex in self.nodeIndices:
            n = self.indexToNode(nodeIndex)
            if not n.isInclusive:
                continue
            if not n.finalStates:
                continue
            daughters = [str(d) for d in n.finalStates]
            elStr += ', (%s(%i) > %s)' % (n, nodesDict[nodeIndex],
                                          ','.join(daughters))
        return elStr

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
            raise SModelSError("Can not convert tree to bracket with %i branches" % len(branches))
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
                        raise SModelSError("Can not convert tree with Z2-violating decays to bracket: \n  %s" % self.treeToString())
                intermediateState[ib].append(str(mom))
                branchList[ib].append(vertexList)
                finalState += fstates

        return branchList, finalState, intermediateState

    def getFinalStates(self, nodeIndex=None):
        """
        Get the list of particles which have not decayed (appear at the top of the tree).
        If source is defined, get the final states generated by the cascade decay of the source
        node. It also sets the finalState attribute for the source node.

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
            for d in self.daughterIndices(nodeIndex):
                finalStates += self.getFinalStates(d)

        # Set final state attribute for node:
        node.finalStates = sorted(finalStates)

        return node.finalStates

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

        # Update dicts:
        self._successors = newSuccessors
        self._predecessors = newPredecessors
        self._nodesMapping = newMapping
        self._nodeCanonNames = newCanonNames
        self._rootIndex = nodeIndexDict[self._rootIndex]

    def updateNodeObjects(self, nodeObjectDict):
        """
        Update the node index -> node object mapping.
        Only affects the indices appearing in nodeObjectDict.

        :param nodeObjectDict: Dictionary with current node indices as keys
                              and new node objects as values
        """

        for nodeIndex, newObj in nodeObjectDict.items():
            self._nodesMapping[nodeIndex] = newObj

    def copyTreeFrom(self,other,nodesObjDict):
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

    def draw(self, particleColor='lightcoral',
             smColor='skyblue',
             pvColor='darkgray',
             nodeScale=4, labelAttr=None, attrUnit=None):
        """
        Draws Tree using matplotlib.

        : param particleColor: color for particle nodes
        : param smColor: color used for particles which have the isSM attribute set to True
        : param pvColor: color for primary vertex
        : param nodeScale: scale size for nodes
        : param labelAttr: attribute to be used as label. If None, will use the string representation
                          of the node object.
        : param attrUnit: Unum object with the unit to be removed from label attribute(if applicable)

        """
        import matplotlib.pyplot as plt
        import networkx as nx

        if labelAttr is None:
            labels = {n: str(n) if not n.isInclusive else n.longStr() for n in self.nodes}
        elif attrUnit is not None:
            labels = {n: str(getattr(n, labelAttr).asNumber(attrUnit))
                      if (hasattr(n, labelAttr) and getattr(n, labelAttr) is not None)
                      else str(n) for n in self.nodes}
        elif labelAttr == 'node':
            labels = {n: str(nodeIndex) for n,nodeIndex in zip(self.nodes,self.nodeIndices)}
        elif labelAttr == 'canonName':
            labels = {n: str(self.nodeCanonName(nodeIndex))
                        for n,nodeIndex in zip(self.nodes,self.nodeIndices)}
        else:
            labels = {n: str(getattr(n, labelAttr)) if hasattr(n, labelAttr)
                      else str(n) for n in self.nodes}

        for key in labels:
            if labels[key] == 'anyOdd':
                labels[key] = 'BSM'
        node_size = []
        node_color = []
        for n in self.nodes:
            node_size.append(nodeScale * 100 * len(labels[n]))
            if 'pv' == labels[n].lower():
                node_color.append(pvColor)
            elif hasattr(n, 'isSM') and n.isSM:
                node_color.append(smColor)
            else:
                node_color.append(particleColor)

        # Compute position of nodes (in case the nodes have the same string representation, first
        # convert the nodes to integers for computing the position)
        G = nx.DiGraph()
        G.add_nodes_from(self.nodes)
        G.add_edges_from(self.edges)
        H = nx.convert_node_labels_to_integers(G, label_attribute='node_label')
        H_layout = nx.drawing.nx_agraph.graphviz_layout(H, prog='dot')
        pos = {H.nodes[n]['node_label']: p for n, p in H_layout.items()}

        nx.draw(G, pos,
                with_labels=True,
                arrows=True,
                labels=labels,
                node_size=node_size,
                node_color=node_color)

        plt.show()
