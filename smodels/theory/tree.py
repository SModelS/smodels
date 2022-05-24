"""
.. module:: graphTools
   :synopsis: A collection of functions used to create and maninuplate graphs.
              Based on networkx.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


import networkx as nx
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.theory.particle import Particle
from smodels.theory.auxiliaryFunctions import bracketToProcessStr
from smodels.tools.inclusiveObjects import InclusiveValue

# Define a common inclusive particle object
# to be used with the inclusive node
IncluviseParticle = Particle(label='Inclusive')


def equalNodes(nodeA, nodeB):
    """
    Check if nodes canonical names and particles are equal.
    (Auxiliary method for networkx DiGraphMatcher.)

    :param nodeA: ParticleNode object or a dictionary with node attributes
    :param nodeB: ParticleNode object or a dictionary with node attributes

    :return: True if nodes are equal, false otherwise
    """

    if nodeA['canonName'] != nodeB['canonName']:
        return False
    if nodeA['particle'] != nodeB['particle']:
        return False

    return True


def sortTreeList(treeList):
    """
    Sort a list of Tree objects according to the compareTrees method.

    :param treeList: List of Tree objects

    :return: Sorted list of Tree objects
    """

    # Create a listed of sorted nodes with proper canonical names according to the node comparison:
    # Always put the inclusive trees at the beginning
    inclusiveTrees = [t for t in treeList if isinstance(t.root, InclusiveParticleNode)]
    sortedTrees = []
    for tree in treeList:
        if isinstance(tree.root, InclusiveParticleNode):
            continue
        itree = 0
        while itree < len(sortedTrees) and compareTrees(tree, sortedTrees[itree]) > 0:
            itree += 1
        sortedTrees.insert(itree, tree)

    sortedTrees = inclusiveTrees + sortedTrees
    return sortedTrees


def compareTrees(treeA, treeB):
    """
    Recursively compare trees. The comparison is based on the canon names
    and particles of the roots. If these are equal, the subtrees are compared.
    The comparison does not take into account permutations of equivalent subtrees.
    If one of the trees has an InclusiveParticleNode as root, it is always considered
    as the smaller tree.

    :param treeA: Tree object
    :param treeB: Tree object

    :return: 1, if treeA > treeB, -1, if treeA < treeB, 0 if treeA == treeB
    """

    canonA = treeA.canonName
    canonB = treeB.canonName
    if canonA != canonB:
        return (-1)**int(canonA < canonB)

    rootA = treeA.root
    rootB = treeB.root
    if isinstance(rootA, InclusiveParticleNode):
        return 0
    elif isinstance(rootB, InclusiveParticleNode):
        return 0

    cmp = rootA.compareNode(rootB)
    # Return comparison if roots differ or root has no daughters
    if cmp != 0:
        return cmp

    # If root has no daughters and roots agree, then trees are equal
    if treeA.out_degree(rootA) == 0:
        return 0

    # Get daughters and sort them, so the comparison is unique:
    subTreesA = [Tree(nx.bfs_tree(treeA, source=daughter))
                 for daughter in treeA.successors(rootA)]
    subTreesA = sortTreeList(subTreesA)
    subTreesB = [Tree(nx.bfs_tree(treeB, source=daughter))
                 for daughter in treeB.successors(rootB)]
    subTreesB = sortTreeList(subTreesB)

    # Loop over sorted daughters and stop when the first differ:
    for itree, subtreeA in enumerate(subTreesA):
        subtreeB = subTreesB[itree]
        cmp = compareTrees(subtreeA, subtreeB)
        if cmp != 0:
            return cmp

    # If all subtrees were identical, return zero
    return 0


class ParticleNode(object):
    """
    Simple wrapper for creating graphs with Particle objects.
    It is necessary because the same particle can appear multiple
    times within a tree, so
    Particle objects can not be directly used as nodes
    (since the same particle can not appear as distinct nodes)

    :ivar particle: Stores the Particle object
    :ivar nodeNumber: Node identifier
    :ivar nodeWeight: Stores the node weight
                      (1 for stable particles, BR for unstable
                      and production xsec for primary vertex)
    """

    _lastNodeNumber = 0

    def __init__(self, particle, nodeNumber=None, nodeWeight=1.0):
        self.particle = particle

        # Since ParticleNodes are identified by their numbering,
        # if it is not specifically assigned, automatically assign
        # a new number which does not overlap with any previous class instances
        if nodeNumber is None:
            self.node = ParticleNode._lastNodeNumber + 1
            ParticleNode._lastNodeNumber += 1
        else:
            self.node = nodeNumber
        ParticleNode._lastNodeNumber = max(self.node,
                                           ParticleNode._lastNodeNumber)
        self.canonName = None
        # Node weight:
        # (1 for stable particles, BR for unstable and
        # production xsec for primary vertex)
        self.nodeWeight = nodeWeight

        # Flag to tag nodes which should not be decayed
        self.finalState = False

    def __hash__(self):
        return self.node

    def __cmp__(self, other):
        if self.node == other.node:
            return 0
        elif self.node > other.node:
            return 1
        else:
            return -1

    def __lt__(self, other):
        return self.__cmp__(other) == -1

    def __gt__(self, other):
        return self.__cmp__(other) == 1

    def __eq__(self, other):
        return self.__cmp__(other) == 0

    def __ne__(self, other):
        return self.__cmp__(other) != 0

    def __str__(self):

        return str(self.particle)

    def __repr__(self):

        return self.__str__()

    def __getitem__(self, attr):
        if attr in self.__dict__:
            return self.__dict__[attr]
        else:
            return getattr(self, attr)

    def __getattr__(self, attr):
        """
        Returns the attribute from particle.

        :parameter attr: Attribute string

        :return: Attribute from self.particle
        """

        return getattr(self.particle, attr)

    def __getstate__(self):
        """
        Since getattr is defined, we must defined getstate
        for pickling/unpickling.
        """

        attrDict = self.__dict__

        return attrDict

    def __setstate__(self, state):
        """
        Since getattr is defined, we must defined getstate
        for pickling/unpickling.
        """

        self.__dict__.update(state)

    def __add__(self, other):
        """
        Adds two nodes. The properties of self are kept, except
        for the particle and nodeWeight, which are added with other.

        :param other: ParticleNode object

        :return: a copy of self with the particle combined with other.particle
                 and nodeWeight added.
        """

        newNode = self.copy()
        newNode.particle = self.particle + other.particle
        newNode.nodeWeight = self.nodeWeight + other.nodeWeight

        return newNode

    def compareNode(self, other):
        """
        Compare nodes accoring to their canonical name
        and particle.

        :param other: ParticleNode or InclusiveParticleNode object

        :return: 1 if self > other, -1 if self < other and 0 if self == other
        """

        if isinstance(other, InclusiveParticleNode):
            return -other.compareNode(self)

        if not isinstance(other, ParticleNode):
            raise SModelSError("Can not compare node to %s" % type(other))

        if self.canonName != other.canonName:
            return (-1)**int(self.canonName < other.canonName)

        if self.particle > other.particle:
            return 1
        elif self.particle < other.particle:
            return -1

        return 0

    def equalTo(self, other):
        """
        Compare nodes accoring to their canonical name
        and particle.

        :param other: ParticleNode or InclusiveParticleNode object

        :return: True if nodes are equal, false otherwise
        """

        return (self.compareNode(other) == 0)

    def copy(self):
        """
        Makes a shallow copy of itself. The particle attribute
        shares the same object with the original copy.
        :return: ParticleNode object
        """

        newNode = ParticleNode(particle=self.particle,
                               nodeNumber=self.node)
        newNode.canonName = self.canonName
        newNode.nodeWeight = self.nodeWeight

        return newNode


class InclusiveParticleNode(ParticleNode):
    """
    An inclusive ParticleNode class. It will return True when compared to any other
    ParticleNode object or InclusiveParticleNode object.

    :ivar particle: IncluviseParticle (dummy)
    :ivar nodeNumber: Node identifier
    :ivar nodeWeight: 1 (dummy value)
    :ivar finalStates: Allowed final states (final state nodes)
    """

    def __init__(self, nodeNumber=None,
                 nodeWeight=1.0,
                 particle=IncluviseParticle,
                 finalStates=[]):
        ParticleNode.__init__(self, particle=particle,
                              nodeNumber=nodeNumber,
                              nodeWeight=nodeWeight)
        self.canonName = InclusiveValue()
        self.finalStates = finalStates[:]

    def compareNode(self, other):
        """
        Return 0 when compared to any ParticleNode.
        If compared to InclusiveParticleNode, compare the
        allowed BSM and SM states.

        :param other: ParticleNode or InclusiveParticleNode object

        :return: 0 (self == other), 1 (self > other), -1 (self < other)
        """

        if isinstance(other, ParticleNode):
            return 0
        elif isinstance(other, InclusiveParticleNode):
            fsA = sorted(self.finalStates)
            fsB = sorted(other.finalStates)
            if fsA != fsB:
                if fsA > fsB:
                    return 1
                else:
                    return -1
        return 0

    def copy(self):
        """
        Makes a shallow copy of itself. The particle attribute
        shares the same object with the original copy.
        :return: ParticleNode object
        """

        newNode = InclusiveParticleNode(nodeNumber=self.node,
                                        particle=self.particle)
        newNode.canonName = self.canonName
        newNode.nodeWeight = self.nodeWeight
        newNode.finalStates = self.finalStates[:]

        return newNode

    def __getattr__(self, attr):
        """
        Returns None if does not contain the attribute.

        :parameter attr: Attribute string

        :return: Attribute value or None
        """
        if attr not in self.__dict__:
            return None


class Tree(nx.DiGraph):
    """
    A wrapper for the tree graphs used to describe elements.
    """

    def __init__(self, info=None, finalState=None,
                 intermediateState=None, model=None):

        self._canonName = None
        self._root = None
        if not info:  # Initialize empty tree
            nx.DiGraph.__init__(self)
        elif isinstance(info, str):
            nx.DiGraph.__init__(self)
            try:
                self.stringToTree(info, finalState, intermediateState, model)
            except (SModelSError, TypeError):
                raise SModelSError("Can not create element from input %s" % info)
        elif isinstance(info, dict):
            nx.DiGraph.__init__(self, info)
        elif isinstance(info, nx.DiGraph):
            nx.DiGraph.__init__(self, info)
        else:
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
                    if isinstance(momNode, InclusiveParticleNode):
                        self.remove_node(node)
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
        T = self
        root = self.root
        nodesDict = {0: 0}
        counter = 1
        for mom, daughters in nx.bfs_successors(T, root):

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
                if not list(T.successors(n)):
                    elStr += '%s,' % n  # stable
                else:
                    if n.node not in nodesDict:
                        nodesDict[n.node] = counter
                        counter += 1
                    elStr += '%s(%i),' % (n, nodesDict[n.node])
            elStr = elStr[:-1] + '), '
        elStr = elStr[:-2]

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

        tree = self
        branches = [b for b in tree.successors(self.root)]
        finalState = []
        intermediateState = []
        branchList = []
        if len(branches) != 2:
            raise SModelSError("Can not convert tree to bracket with %i branches" % len(branches))
        for ib, b in enumerate(branches):
            intermediateState.append([])
            branchList.append([])
            # Deal separately with the case where the primary mother is stable:
            if tree.out_degree(b) == 0:
                if not b.isSM:
                    finalState.append(str(b))
                    continue
                else:
                    raise SModelSError("Can not convert tree with Z2-violating decays to bracket")
            for mom, daughters in nx.bfs_successors(tree, b):
                vertexList = [str(d) for d in daughters if d.isSM]
                fstates = [str(d) for d in daughters if not d.isSM and tree.out_degree(d) == 0]
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
            self._canonName = self.setCanonName()

        return self._canonName

    def setCanonName(self, node=None):
        """
        Recursively compute the canonical name for the Tree T.
        Returns the name in integer form.

        :param T: Tree (Tree object)
        :param node: Node to get the name for. If None, it will use the root

        :return: Integer representing the Tree
        """

        if not self.number_of_nodes():
            return None

        if node is None:
            node = self.root

        # If it is inclusive node set its name to an inclusive integer
        # and return its name (no need to check the children)
        if isinstance(node, InclusiveParticleNode):
            node.canonName = InclusiveValue()
            return node.canonName

        children = list(self.successors(node))
        if not children:
            node.canonName = 10
        else:
            tp = [self.setCanonName(n) for n in children]
            if any(isinstance(name, InclusiveValue) for name in tp):
                node.canonName = InclusiveValue()
            else:
                tp = sorted(tp)
                tpStr = '1' + "".join(str(c) for c in tp) + '0'
                node.canonName = int(tpStr)

        # Use the root canon name to assign the tree canon name
        if self.in_degree[node] == 0:
            self._canonName = node.canonName

        return node.canonName

    def getFinalStates(self):
        """
        Get the list of particles which have not decayed (appear at the top of the tree)

        :returns: list of ParticleNode objects
        """

        finalStates = [n for n in self.nodes if self.out_degree(n) == 0]
        return finalStates

    @property
    def root(self):
        """
        Returns the root node (primary vertex) of the tree.
        If it has not been defined, compute it.

        :return: root node
        """

        if self._root is None:
            root = [n for n in self.nodes if self.in_degree[n] == 0]
            if len(root) != 1:
                raise SModelSError("Malformed Tree, %i root(s) have been found." % len(root))
            self._root = root[0]

        return self._root

    def getTreeWeight(self):
        """
        Computes the tree weight (production cross-section*BRs). If it can not be computed,
        return None. Does not include the weight of final state (undecayed) particles.

        :return: tree weight (Unum object) if available, None otherwise.
        """

        weight = 1.0
        for n in self.nodes:
            # Skip final state (stable) particles
            if self.out_degree(n) == 0:
                continue
            if not hasattr(n, 'nodeWeight') or n.nodeWeight is None:
                return None
            weight *= n.nodeWeight

        return weight

    def bfs_successors(self, node=None):
        """
        Returns an iterator over the mother and daughter
        nodes starting at node using a breadth first search.

        :param node: Node from tree. If None, starts at tree root.

        :return: Iterator over nodes.
        """

        if node is None:
            node = self.root

        return nx.bfs_successors(self, node)

    def compareTreeTo(self, other):
        """
        Check self equals other.
        The comparison is based on the node canonical name (node structure)
        and the node particles. It is done recursively, so two nodes are equal
        only if all their daughters are also equal.
        For daughters with the same canonical name, all the permutations are tried
        when comparing nodes.
        If nodes match, return a tree (Tree obj) with nodeB as root
        and its daughters as nodes, but with the daughters sorted according
        to the ordering to which they matched nodeA.

        :param other: a tree (Tree object)

        :return: (True, new tree) if nodes match, (False, None) otherwise.
        """

        # Compare canon names and make sure they are defined
        canonA = self.canonName
        canonB = other.canonName
        if canonA != canonB:
            return (-1)**int(canonA < canonB), None

        # Compare root nodes:
        rootA = self.root
        rootB = other.root
        cmpNode = rootA.compareNode(rootB)
        if cmpNode != 0:
            return cmpNode, None

        # Make sure node attributes are set and updated:
        attrs = {n: n.__dict__ for n in self.nodes}
        nx.set_node_attributes(self, attrs)
        attrs = {n: n.__dict__ for n in other.nodes}
        nx.set_node_attributes(other, attrs)

        # Check for inclusive tree:
        if any(isinstance(d, InclusiveParticleNode) for d in self.successors(rootA)):
            # Check if self is a subgraph of other
            matcher = nx.algorithms.isomorphism.DiGraphMatcher(other, self,
                                                               node_match=equalNodes)
            matches = matcher.subgraph_isomorphisms_iter()
        elif any(isinstance(d, InclusiveParticleNode) for d in other.successors(rootB)):
            # Check if other is a subgraph of self
            matcher = nx.algorithms.isomorphism.DiGraphMatcher(self, other,
                                                               node_match=equalNodes)
            matches = matcher.subgraph_isomorphisms_iter()
        else:
            matcher = nx.algorithms.isomorphism.DiGraphMatcher(self, other,
                                                               node_match=equalNodes)
            # Check if there is any isomorphism which matches the trees
            # using the equalNodes function to compare node attributes
            matches = matcher.match()

        # If there is a match, return the first match
        # with a copy of the other and its nodes ordered according
        # to self.
        for match in matches:
            # Make a copy of other
            matchedTree = other.copyTree()
            # Remove all nodes and edges
            matchedTree.clear()
            # Store the nodes in the correct order:
            newNodes = []
            for n in self.nodes:
                newNode = match[n].copy()  # match = {self.node : other.node}
                newNode.node = n.node
                newNodes.append(newNode)
                match[n] = newNode  # Update the dictionary with the copied node
            # Store the edges in the correct order:
            newEdges = [(match[nA], match[nB]) for nA, nB in self.edges]
            # Add the nodes and the edges to the new tree:
            matchedTree.add_nodes_from(newNodes)
            matchedTree.add_edges_from(newEdges)
            return 0, matchedTree

        # If no matches were found, compare the trees:
        cmp = compareTrees(self, other)
        # We should always have cmp != 0, since no isomorphism was found.
        # Sanity check:
        if cmp == 0:
            raise SModelSError("Trees are equal, but no isomorphism was found. Something went wrong.")
        return cmp, None

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
        for daughter in list(self.successors(self.root)):
            subTree = Tree(nx.bfs_tree(self, source=daughter))
            subTree.sort()
            subTrees.append(subTree)
            self.remove_nodes_from(subTree.nodes())

        # Sort the subtrees
        subTrees = sortTreeList(subTrees)
        # Add sorted nodes:
        for subTree in subTrees:
            self.add_node(subTree.root)

        # Add sorted edges:
        for subTree in subTrees:
            self.add_edge(self.root, subTree.root)
            self.add_edges_from(subTree.edges)

    def numberNodes(self):
        """
        Renumber the nodes from 0,N following a breadth-first search
        for the tree.
        """

        newTree = nx.bfs_tree(self, source=self.root)
        nodeList = list(newTree.nodes)
        edgesList = list(newTree.edges)

        # Re-number nodes:
        inode = 0
        for node in nodeList:
            node.node = inode
            inode += 1

        # Finally update tree with the node dictionary:
        nx.to_networkx_graph(edgesList, create_using=self)

    def copyTree(self):
        """
        Returns a copy of self. The copy contains the same canonical name
        and copies of nodes.

        : return: Tree object
        """

        newNodesDict = {n: n.copy() for n in self.nodes}
        newTree = nx.relabel_nodes(self, newNodesDict, copy=True)
        newTree._canonName = self._canonName

        return newTree

    def checkConsistency(self):
        """
        Make sure the tree has the correct topology(directed rooted tree).
        Raises an error otherwise.
        """
        if self.number_of_nodes() and not nx.is_arborescence(self):
            raise SModelSError("Graph created with malformed structure (not  a tree).")

    def addTrees(self, other):
        """
        Combines the nodes(add particles) in equivalent nodes in each trees.
        Can only be done if the trees have the same topology and ordering.

        : param other: tree(Tree object)

        : return: new tree with the combined particles(Tree object)
        """

        nodesB = list(other.nodes)
        nodesDict = {n: n + nodesB[inode]
                     for inode, n in enumerate(self.nodes)}
        newTree = nx.relabel_nodes(self, nodesDict, copy=True)

        return newTree

    def compose(self, other):
        """
        Returns a new graph of self composed with other.
        The trees are joined where the nodes overlap and the nodes
        from self are kept.

        : param other: tree(Tree object)

        : return: new tree with the other composed with self.
        """

        return nx.compose(self, other)

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

        if labelAttr is None:
            labels = {n: str(n) for n in self.nodes}
        elif attrUnit is not None:
            labels = {n: str(getattr(n, labelAttr).asNumber(attrUnit)) if hasattr(n, labelAttr)
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
        H = nx.convert_node_labels_to_integers(self, label_attribute='node_label')
        H_layout = nx.drawing.nx_agraph.graphviz_layout(H, prog='dot')
        pos = {H.nodes[n]['node_label']: p for n, p in H_layout.items()}

        nx.draw(self, pos,
                with_labels=True,
                arrows=True,
                labels=labels,
                node_size=node_size,
                node_color=node_color)

        plt.show()
