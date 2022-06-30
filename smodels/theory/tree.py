"""
.. module:: tree
   :synopsis: Classes used to construct trees (root directed graphs) describing
              simplified model topologies.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.theory.auxiliaryFunctions import bracketToProcessStr
from smodels.tools.inclusiveObjects import InclusiveValue
from smodels.theory.treeMatcher import TreeMatcher, sortTreeList
from smodels.theory.particleNode import ParticleNode, InclusiveParticleNode
from smodels.tools.physicsUnits import fb
from collections import OrderedDict
from itertools import product
import unum


class Tree(object):
    """
    A wrapper for the tree graphs used to describe elements.
    """

    def __init__(self, info=None, finalState=None,
                 intermediateState=None, model=None):

        self._canonName = None
        self._root = None
        self._weight = None
        self.successors = OrderedDict()  # Stores the nodes and their successors (daughters)
        self.predecessors = {}  # Stores the nodes and their predecessors (parents)
        self.nodesMapping = {}  # Stores the mapping nodeIndex->node object

        # Convert from string:
        if isinstance(info, str):
            try:
                self.stringToTree(info, finalState, intermediateState, model)
            except (SModelSError, TypeError):
                raise SModelSError("Can not create element from input %s" % info)
        # Convert from dictionary with edges ({node : [node1,node2,..], ...}):
        elif isinstance(info, dict):
            self.add_nodes_from(info.keys())
            for node1, nodeList in info.items():
                self.add_edges_from(product([node1], nodeList))
        # Convert from Tree or DiGraph objec:
        elif isinstance(info, Tree):
            self.add_nodes_from(info.nodes)
            self.add_edges_from(info.edges)
        elif info is not None:
            raise SModelSError("Can not create element from input type %s" % type(info))

    def __add__(self, other):
        """
        Combines the nodes (add particles) in equivalent nodes in each trees.
        Can only be done if the trees have the same topology and ordering.

        :param other: tree (Tree object)

        :return: new tree with the combined particles (Tree object)
        """

        return self.addTrees(other)

    def __repr__(self):
        return str(self)

    def __str__(self):
        return self.treeToString()

    def add_node(self, node):
        """
        Adds a node to the Tree if the nodeIndex (node.node)
        is not in the Tree.

        :param node: ParticleNode object
        """

        nodeIndex = node.node   # Use node number for fast indexing
        if nodeIndex not in self.successors:
            self.successors[nodeIndex] = []
            self.nodesMapping[nodeIndex] = node

    def add_nodes_from(self, nodes):
        """
        Adds a list of nodes to the Tree.

        :param nodes: List of ParticleNode objects
        """

        for node in nodes:
            self.add_node(node)

    def remove_node(self, nodeIndex):
        """
        Removes a node from the Tree if the nodeIndex (node.node).
        The node is removed as well as its appearing in any edges.

        :param nodeIndex: Node index
        """

        if nodeIndex in self.successors:
            self.successors.pop(nodeIndex)
            self.nodesMapping.pop(nodeIndex)

        for nodeA, daughtersA in self.successors.items():
            if nodeIndex not in daughtersA:
                continue
            daughtersA = [d for d in daughtersA[:] if d != nodeIndex]
            self.successors[nodeA] = daughtersA[:]

        if nodeIndex in self.predecessors:
            self.predecessors.pop(nodeIndex)

        for nodeA, momA in list(self.predecessors.items()):
            if momA == nodeIndex:
                self.predecessors.pop(nodeA)

    def remove_nodes_from(self, nodeIndices):
        """
        Removes a list of nodes from the Tree.

        :param nodeIndices: List of node indices
        """

        for nodeIndex in nodeIndices:
            self.remove_node(nodeIndex)

    def add_edge(self, nodeA, nodeB):
        """
        Adds a directed edge to the Tree (nodeA -> nodeB).
        If the nodes were not present in the Tree, they are also added.

        :param nodeA: ParticleNode object
        :param nodeB: ParticleNode object
        """

        self.add_node(nodeA)
        self.add_node(nodeB)
        nodeIndexA = nodeA.node
        nodeIndexB = nodeB.node
        self.successors[nodeIndexA].append(nodeIndexB)
        self.predecessors[nodeIndexB] = nodeIndexA

    def add_edges_from(self, edges):
        """
        Adds a list of directed edges to the Tree.

        :param edges: List of tuples containing ParticleNode objects
                      (e.g. [(nodeA,nodeB),(nodeA,nodeC),...])
        """

        for edge in edges:
            self.add_edge(edge[0], edge[1])

    def clear(self):
        """
        Remove all nodes and edges from the graph, but
        keep its canonName.
        """

        self.successors = OrderedDict()
        self.predecessors = {}
        self.nodesMapping = {}
        self._root = None

    @property
    def nodes(self, indices=False):
        """
        Returns the tist of ParticleNode objects in the Tree.

        :param indices: If True, return the node indices, instead of
                        the node objects.

        :return: List of ParticleNode objects
        """

        if not indices:
            nodeList = [self.nodesMapping[n] for n in self.successors]
        else:
            nodeList = list(self.successors)

        return nodeList

    @property
    def edges(self, indices=False):
        """
        Returns the list of edges (pairs of ParticleNode objects) in the Tree.

        :param indices: If True, return the node indices, instead of
                        the node objects.

        :return: List of edges
        """

        edgesList = []
        for node, nodeList in self.successors.items():
            for nodeB in nodeList:
                if not indices:
                    edgesList.append((self.nodesMapping[node], self.nodesMapping[nodeB]))
                else:
                    edgesList.append((node, nodeB))
        return edgesList

    def out_degree(self, nodeIndex):
        """
        Computes the number of outgoing edges from the node
        (number of daughters).

        :param nodeIndex: Node index (int)

        :return: Number of outgoing edges (int)
        """

        if nodeIndex not in self.successors:
            return 0
        else:
            return len(self.successors[nodeIndex])

    def in_degree(self, nodeIndex):
        """
        Computes the number of incoming edges to the node
        (number of parents).

        :param nodeIndex: Node index (int)

        :return: Number of incoming edges (1 or 0)
        """

        if nodeIndex in self.predecessors:
            if self.predecessors[nodeIndex] is not None:
                return 1

        return 0

    def number_of_nodes(self):
        """
        Returns the total number of nodes in the Tree.


        :return: Number of nodes (int)
        """

        return len(self.successors)

    def bfs_successors(self, node=None,
                       includeLeaves=False, indices=False):
        """
        Returns an iterator over the mother and daughter
        nodes starting at node using a breadth first search.

        :param node: Node from tree. If None, starts at tree root.
        :param includeLeaves: If True, it will consider the leaves (undecayed nodes)
                              as moms in the iterator (with an empty daughters list)
        :param indices: If True, return the node indices, instead of
                        the node objects.


        :return: Iterator over nodes.
        """

        if node is None:
            node = self.root

        mom = node.node
        daughters = self.successors[mom]
        generation = [(mom, daughters)]
        while generation:
            for pair in generation:
                if not indices:
                    yield (self.nodesMapping[pair[0]], [self.nodesMapping[d] for d in pair[1]])
                else:
                    yield pair
            next_generation = []
            for pair in generation:
                mom, daughters = pair
                for new_mom in daughters:
                    if new_mom not in self.successors:
                        continue
                    new_daughters = self.successors[new_mom]
                    if new_daughters or includeLeaves:
                        next_generation.append((new_mom, new_daughters))
            generation = next_generation

    def dfs_successors(self, node=None,
                       includeLeaves=False, indices=False):
        """
        Returns a dictionary with the mother as keys
        and the daughters as values using a depth-first search.
        nodes starting at node using a breadth first search.

        :param node: Node from tree. If None, starts at tree root.
        :param includeLeaves: If True, it will consider the leaves (undecayed nodes)
                              as moms in the iterator (with an empty daughters list)
        :param indices: If True, return the node indices, instead of
                        the node objects.

        :return: Dictionary with mothers and daughters
        """

        if node is None:
            node = self.root

        mom = node.node
        daughters = self.successors[mom][:]
        if not indices:
            daughtersDict = {self.nodesMapping[mom]: [self.nodesMapping[d] for d in daughters]}
        else:
            daughtersDict = {mom: daughters}

        # Go down in generation until it has no more daughters
        while daughters:
            new_mom = daughters.pop(0)
            new_daughters = self.successors[new_mom][:]
            if new_daughters or includeLeaves:
                if not indices:
                    daughtersDict.update({self.nodesMapping[new_mom]:
                                          [self.nodesMapping[d] for d in new_daughters]})
                else:
                    daughtersDict.update({new_mom: daughters[:]})

            # Attach new daughters to the beginning of the list
            daughters = new_daughters[:] + daughters[:]

        return daughtersDict

    def bfs_tree(self, node=None):
        """
        Returns an oriented tree constructed from of a breadth-first-search
        starting at node.
        """

        if node is None:
            node = self.root

        T = Tree()
        T.nodesMapping = {nodeIndex: n for nodeIndex, n in self.nodesMapping}
        for mom, daughters in self.bfs_successors(node, includeLeaves=True, indices=True):
            T.successors[mom] = daughters
            for d in daughters:
                T.predecessors[d] = mom
        return T

    def stringToTree(self, stringElement, finalState=None,
                     intermediateState=None, model=None):
        """
        Converts a string describing an element to a Tree (Tree object). It accepts
        the (old) bracket notation or the process notation. For the old notation the
        optional arguments finalState and intermediateState can also be defined.
        If the argument model is defined, the particle labels will be converted to
        Particle objects from the model. Otherwise the nodes will hold the particle strings.

        :param stringElement: The process in string format
                              (e.g. '(PV > gluino(1),squark(2)), (gluino(1) >
                               MET,jet,jet), (squark(2) > HSCP,u)' or [[['jet','jet']],[['u']]]).
                               The particle labels should match the particles in the Model
                               (if Model != None).

        :parameter model: The model (Model object) to be used when converting particle labels to
                          particle objects. If None, the nodes will only store the particle labels.

        :parameter finalState: (optional) list containing the final state labels for each branch
                               (e.g. ['MET', 'HSCP'] or ['MET','MET'])
        :parameter intermediateState: (optional) nested list containing intermediate state labels
                                         for each branch  (e.g. [['gluino'], ['gluino']])

        """

        # First check if string is in old format:
        if '[' in stringElement and ']' in stringElement:
            procString = bracketToProcessStr(stringElement, finalState=finalState,
                                             intermediateState=intermediateState)
        elif '>' in stringElement and 'PV' in stringElement:
            procString = stringElement
        else:
            raise SModelSError("Could not recognize string format for element (%s)" % stringElement)

        decays = procString.replace(" ", "").split("),(")
        decays[0] = decays[0][1:]  # Remove trailing parenthesis
        decays[-1] = decays[-1][:-1]  # Remove remaining parenthesis

        # Build a dictionary with all unstable particles:
        nodesDict = {}
        for dec in decays:
            ptcs = [dec.split('>')[0].strip()] + [p.strip() for p in dec.split('>')[1].split(',')]
            for ptc in ptcs:
                if ptc == 'PV':
                    n = 0  # PV is always node zero
                elif '(' in ptc and ')' in ptc:  # Unstable particles should always have a unique numbering
                    n = eval(ptc.split('(')[1].split(')')[0])
                else:
                    continue  # Stable particles will have their nodes defined later
                nodesDict[ptc] = n

        # Make sure particles have unique nodes:
        if len(set(list(nodesDict.values()))) != len(list(nodesDict.values())):
            raise SModelSError("Input string has non unique nodes: %s" % nodesDict)

        # Find maximum node defined so far:
        maxNode = max([n for n in nodesDict.values() if n is not None])

        # First add all nodes:
        edges = []
        for dec in decays:
            mom = dec.split('>')[0].strip()
            daughters = [p.strip() for p in dec.split('>')[1].split(',')]
            ptcs = [mom] + daughters
            for iptc, ptc in enumerate(ptcs):
                if ptc in nodesDict:
                    n = nodesDict[ptc]
                else:
                    n = maxNode + 1
                    maxNode = n
                particleLabel = ptc.replace('(%i)' % n, '')
                # Check for inclusive node:
                if particleLabel.lower() == 'inclusivenode':
                    node = InclusiveParticleNode(nodeNumber=n)
                else:
                    if model is None:
                        particle = particleLabel
                    else:
                        particle = model.getParticlesWith(label=particleLabel)
                        if not particle:
                            raise SModelSError("Final state %s has not been defined in model %s"
                                               % (particleLabel, model))
                        elif len(particle) != 1:
                            raise SModelSError("Ambiguos defintion of label %s in model %s"
                                               % (particleLabel, model))
                        else:
                            particle = particle[0]
                    node = ParticleNode(particle=particle, nodeNumber=n)
                self.add_node(node)  # only added if node is not in the tree
                if iptc == 0:
                    # Make sure to use the correct node object
                    momNode = [n for n in self.nodes if n == node][0]
                else:
                    # For inclusive nodes, add decays as possible final states
                    # and do not include nodes in the tree
                    if momNode.isInclusive:
                        self.remove_node(node.node)
                        momNode.finalStates.append(node.particle)
                    else:
                        edges.append((momNode, node))
        # Now add all edges
        self.add_edges_from(edges)

    def treeToString(self):
        """
        Convert the tree to a process string (e.g. '(PV > gluino(1),squark(2)), (gluino(1) >
                           MET,jet,jet), (squark(2) > HSCP,u)')

        :return: String describing the process
        """

        elStr = ""
        root = self.root
        nodesDict = {0: 0}
        counter = 1
        for mom, daughters in self.bfs_successors(root):

            # Add mom (if mom = 0 does not include index)
            if mom.node == 0:
                elStr += '(%s > ' % mom
            else:
                if mom.node not in nodesDict:
                    nodesDict[mom.node] = counter
                    counter += 1
                elStr += '(%s(%i) > ' % (mom, nodesDict[mom.node])

            # Add daughters for the decay, if daughter is unstable, add index
            for n in daughters:
                if n.isInclusive and n.finalStates:
                    stable = False
                elif self.out_degree(n.node) == 0:
                    stable = True
                else:
                    stable = False
                if stable:
                    elStr += '%s,' % n  # stable
                else:
                    if n.node not in nodesDict:
                        nodesDict[n.node] = counter
                        counter += 1
                    elStr += '%s(%i),' % (n, nodesDict[n.node])
            elStr = elStr[:-1] + '), '
        elStr = elStr[:-2]

        # Deal separately with daughters for inclusive nodes (if it exists)
        for n in self.nodes:
            if not n.isInclusive:
                continue
            if not n.finalStates:
                continue
            daughters = [str(d) for d in n.finalStates]
            elStr += ', (%s(%i) > %s)' % (n, nodesDict[n.node],
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

        branches = [self.nodesMapping[b] for b in self.successors[self.root.node]]
        finalState = []
        intermediateState = []
        branchList = []
        if len(branches) != 2:
            raise SModelSError("Can not convert tree to bracket with %i branches" % len(branches))
        for ib, b in enumerate(branches):
            intermediateState.append([])
            branchList.append([])
            # Deal separately with the case where the primary mother is stable:
            if self.out_degree(b.node) == 0:
                if not b.isSM:
                    finalState.append(str(b))
                    continue
                else:
                    raise SModelSError("Can not convert tree with Z2-violating decays to bracket")
            for mom, daughters in self.bfs_successors(b):
                vertexList = [str(d) for d in daughters if d.isSM]
                fstates = [str(d) for d in daughters
                           if not d.isSM and self.out_degree(d.node) == 0]
                if daughters:
                    if len(vertexList) != len(daughters) - 1 or len(fstates) > 1:
                        raise SModelSError("Can not convert tree with Z2-violating decays to bracket: \n  %s" % self.treeToString())
                intermediateState[ib].append(str(mom))
                branchList[ib].append(vertexList)
                finalState += fstates

        return branchList, finalState, intermediateState

    @property
    def canonName(self):
        """
        Returns the canonName. If not defined, it will be computed.

        :return: Canonical name (int)
        """

        if self._canonName is None:
            self._canonName = self.setGlobalProperties()

        return self._canonName

    def setGlobalProperties(self, nodeIndex=None):
        """
        Define useful global properties for the trees and the nodes.
        Recursively sets the canonName and finalStates for each node.
        Returns the canonical name in integer form.

        :param nodeIndex: Node index to set the name for. If None, it will use the root

        :return: Integer representing the Tree canonical name
        """

        if not self.number_of_nodes():
            return None

        if nodeIndex is None:
            nodeIndex = self.root.node

        # Set the final state
        node = self.nodesMapping[nodeIndex]
        node.finalStates = self.getFinalStates(node)
        # If it is inclusive node set its name to an inclusive integer
        # and return its name (no need to check the children)
        if node.isInclusive:
            node.canonName = InclusiveValue()
            return node.canonName

        children = self.successors[nodeIndex]
        if not children:
            node.canonName = 10
        else:
            tp = [self.setGlobalProperties(n) for n in children]
            if any(isinstance(name, InclusiveValue) for name in tp):
                node.canonName = InclusiveValue()
            else:
                tp = sorted(tp)
                tpStr = '1' + "".join(str(c) for c in tp) + '0'
                node.canonName = int(tpStr)

        # Use the root canon name to assign the tree canon name
        if not self.in_degree(nodeIndex):
            self._canonName = node.canonName

        return node.canonName

    def getFinalStates(self, source=None):
        """
        Get the list of particles which have not decayed (appear at the top of the tree).
        If source is defined, get the final states generated by the cascade decay of the source
        node.

        :param source: ParticleNode belonging to the tree.

        :returns: list of Particle objects
        """

        if source is None:
            source = self.root

        # For inclusive nodes return directly it final state attribute
        if source.isInclusive:
            return sorted(source.finalStates)

        # If source is root, it is quicker to get all final states directly
        if source is self.root:
            finalStates = [n.particle for n in self.nodes if self.out_degree(n.node) == 0]
            return sorted(finalStates)

        # Return source, if it is already a final state
        if not self.successors[source.node]:
            return [source.particle]

        finalStates = []
        for mom, daughters in self.bfs_successors(source):
            finalStates += [d.particle for d in daughters
                            if not self.successors[d.node]]

        return sorted(finalStates)

    @property
    def root(self):
        """
        Returns the root node (primary vertex) of the tree.
        If it has not been defined, compute it.

        :return: root node
        """

        if self._root is None:
            root = [n for n in self.nodes if self.in_degree(n.node) == 0]
            if len(root) != 1:
                print(self.successors)
                print(self.predecessors)
                raise SModelSError("Malformed Tree, %i root(s) have been found." % len(root))
            self._root = root[0]

        return self._root

    def getTreeWeightList(self):
        """
        Computes the tree weight (production cross-section*BRs). If it can not be computed,
        return None. Does not include the weight of final state (undecayed) particles.

        :return: CrossSectionList object or Unum object, if available, None otherwise.
        """

        weight = 1.0
        for n in self.nodes:
            # Skip final state (stable) particles
            if self.out_degree(n.node) == 0:
                continue
            weight *= n.nodeWeight

        return weight

    def getTreeWeight(self):
        """
        Computes the tree maximum weight (max production cross-section*BRs). If it can not be computed,
        return None. Does not include the weight of final state (undecayed) particles.

        :return: tree weight in fb (float) if available, None otherwise.
        """

        if self._weight is not None:
            return self._weight

        weight = self.getTreeWeightList()

        if isinstance(weight, (int, float)):
            self._weight = weight
        elif isinstance(weight, unum.Unum):
            self._weight = weight.asNumber(fb)
        else:
            self._weight = weight.getMaxXsec().asNumber(fb)
        return self._weight

    def getSubTree(self, source):
        """
        Returns the subtree with source as its root.

        :param source: Node
        :return: Tree object with source as its root node
        """

        subTree = Tree()
        subTree.nodesMapping = {nodeIndex: node for nodeIndex, node in self.nodesMapping.items()}
        subTree._root = source
        for mom, daughters in self.bfs_successors(source, indices=True, includeLeaves=True):
            subTree.successors[mom] = daughters[:]
            for d in daughters:
                subTree.predecessors[d] = mom

        return subTree

    def compareTreeTo(self, other):
        """
        Check self equals other.
        The comparison is based on the node comparison.
        If trees match, return a copy of other with its nodes
        sorted according to self.
        If self or other have inclusive nodes, the matching
        subtree is returned.

        :param other: a tree (Tree object)

        :return: (True, new tree) if nodes match, (False, None) otherwise.
        """

        # Compare canon names and make sure they are defined
        canonA = self.canonName
        canonB = other.canonName
        if canonA != canonB:
            if canonA > canonB:
                return 1, None
            else:
                return -1, None

        # Try to find an isormorphism:
        matcher = TreeMatcher(self, other)
        cmp, matchedTree = matcher.compareTrees()
        return cmp, matchedTree

    def sort(self):
        """
        Sort tree. All the subtrees are sorted according to compareTrees.
        The node numbering is not modified.
        """

        cName = self.canonName  # Just to make sure canonName is defined
        if cName is None:
            return

        # Get all subtrees formed by the daughters of the
        # current root node and sort each subtree
        subTrees = []
        for daughter in self.successors[self.root.node][:]:
            subTree = self.getSubTree(source=self.nodesMapping[daughter])
            subTree.sort()
            subTrees.append(subTree)
            self.remove_nodes_from(subTree.successors.keys())

        # Sort the subtrees
        subTrees = sortTreeList(subTrees)
        # Add sorted nodes:
        self.add_nodes_from([subTree.root for subTree in subTrees])
        # Add sorted edges:
        for subTree in subTrees:
            self.add_edge(self.root, subTree.root)
            self.add_edges_from(subTree.edges)

    def numberNodes(self):
        """
        Renumber the nodes from 0,N following a breadth-first search
        for the tree.
        """

        # Get
        nodeOrder = []
        for mom, daughters in self.bfs_successors(includeLeaves=True):
            nodeOrder.append((mom, daughters))

        # Re-number nodes:
        inode = 0
        for node, daughters in nodeOrder:
            node.node = inode
            inode += 1

        new_successors = OrderedDict()
        new_predecessors = {}
        newMapping = {}
        for node, daughters in nodeOrder:
            daughters = sorted(daughters, key=lambda d: d.node)
            new_successors[node.node] = [d.node for d in daughters]
            for d in daughters:
                new_predecessors[d.node] = node.node
            newMapping[node.node] = node

        self.successors = new_successors
        self.predecessors = new_predecessors
        self.nodesMapping = newMapping

    def copyTree(self, emptyNodes=False):
        """
        Returns a copy of self. The copy contains the same canonical name
        and copies of nodes.

        :param emptyNodes: If True, will not define nodes to the tree

        :return: Tree object
        """

        newTree = Tree()
        if not emptyNodes:
            newTree.successors.update({n: daughters[:]
                                       for n, daughters in self.successors.items()})
            newTree.predecessors = {k: v for k, v in self.predecessors.items()}
            newTree.nodesMapping = {n.node: n.copy() for n in self.nodes}
            if self._root is not None:
                newTree._root = newTree.nodesMapping[self._root.node]
        if self._canonName:
            newTree._canonName = self._canonName

        newTree._weight = self._weight

        return newTree

    def compressToFinalStates(self):
        """
        Compress the Tree to its final states. After the compression
        the tree will have one root (PV). The root's daughters
        are the final state nodes of self.

        :returns: compressed copy of the Tree.
        """

        newTree = self.copyTree()
        # Get root:
        root = newTree.root
        # Get final state nodes:
        fsNodes = [node for node in newTree.nodes if newTree.out_degree(node.node) == 0]
        fsNodes = sorted(fsNodes, key=lambda node: node.particle)
        # Remove all nodes
        newTree.clear()

        # Add root > fsNode edges:
        edges = [(root, fsNode) for fsNode in fsNodes]
        newTree.add_edges_from(edges)
        newTree.setGlobalProperties()

        return newTree

    def checkConsistency(self):
        """
        Make sure the tree has the correct topology (directed rooted tree).
        Raises an error otherwise.
        """

        malformedTree = False
        if len(self.successors) > 1:
            rootIndex = self.root.node

            # Check if root has no parents and at least one daughter
            if self.in_degree(rootIndex) != 0 or self.out_degree(rootIndex) == 0:
                malformedTree = True

            # Check if all nodes (except root) have a unique parent
            if any(self.in_degree(nodeIndex) != 1 for nodeIndex in self.successors
                   if nodeIndex != rootIndex):
                malformedTree = True

            # Check if all nodes can be reached from the root node
            if len(self.dfs_successors(includeLeaves=True)) != len(self.successors):
                malformedTree = True

        if malformedTree:
            raise SModelSError("Graph created with malformed structure (not  a tree).")

    def addTrees(self, other):
        """
        Combines the nodes(add particles) in equivalent nodes in each trees.
        Can only be done if the trees have the same topology and ordering.

        : param other: tree(Tree object)

        : return: new tree with the combined particles(Tree object)
        """

        if other.canonName != self.canonName:
            raise SModelSError("Can not add trees with distinct topologies")

        nodesB = other.nodes
        newMapping = {n.node: n + nodesB[inode]
                      for inode, n in enumerate(self.nodes)}

        newTree = Tree()
        newTree._canonName = self._canonName
        newTree._root = newMapping[self._root.node]
        newTree.successors.update({k: v[:] for k, v in self.successors.items()})
        newTree.predecessors.update({k: v for k, v in self.predecessors.items()})
        # Just change the nodeIndex -> node mapping:
        newTree.nodesMapping = newMapping

        return newTree

    def attachTo(self, other, atNode=None):
        """
        Returns a copy of other with self attached to it at node atNode.
        If atNode is not defined, the trees are joined where the nodes overlap.
        The overlapping node from self is kept.

        :param other: tree (Tree object)
        :param atNode: node (ParticleNode object) where the merge has to take place
                        If not defined, will find the common node between the trees.

        : return: new tree with the other composed with self.
        """

        if atNode is None:
            # Find intersection:
            self_nodes = set(list(self.nodes))
            other_nodes = set(list(other.nodes))
            common_nodes = self_nodes.intersection(other_nodes)
            if len(common_nodes) != 1:
                raise SModelSError("Can not merge trees. %i common nodes found" % len(common_nodes))

            atNode = list(common_nodes)[0]  # merge node of self
            # Make sure the node from self replaces the equivalent node from other
            atNode = self.nodesMapping[atNode.node]

        if not any(n is atNode for n in self.nodes):
            raise SModelSError("Attach node must belong to self")

        if other.out_degree(atNode.node) != 0:
            raise SModelSError("Can not attach tree. Common node can not have daughters in the other tree.")
        if self.in_degree(atNode.node) != 0:
            raise SModelSError("Can not attach tree. Common node can not have parents in self.")

        newTree = other.copyTree()
        # Replace node from base tree by atNode:
        newTree.nodesMapping[atNode.node] = atNode

        for mom, daughters in self.bfs_successors(atNode, includeLeaves=True):
            newTree.add_node(mom)
            newTree.add_edges_from(product([mom], daughters))

        # Update weight
        newTree._weight = self.getTreeWeight()*other.getTreeWeight()

        return newTree

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
