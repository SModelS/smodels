"""
.. module:: graphTools
   :synopsis: A collection of functions used to create and maninuplate graphs. Based on networkx.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


import networkx as nx
import itertools
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.theory.auxiliaryFunctions import bracketToProcessStr


class ParticleNode(object):
    """
    Simple wrapper for creating graphs with Particle objects.
    It is necessary because the same particle can appear multiple times within a tree, so
    Particle objects can not be directly used as nodes
    (since the same particle can not appear as distinct nodes)

    :ivar particle: Stores the Particle object
    :ivar nodeNumber: Node identifier
    :ivar nodeWeight: Stores the node weight
                      (1 for stable particles, BR for unstable and production xsec for primary vertex)
    """

    _lastNodeNumber = 0

    def __init__(self, particle, nodeNumber=None, nodeWeight=1.0):
        self.particle = particle

        # Since ParticleNodes are identified by their numbering,
        # if it is not specifically assigned, automatically assign
        # a new number which does not overlap with any previous class instances
        if nodeNumber is None:
            self.node = ParticleNode._lastNodeNumber+1
            ParticleNode._lastNodeNumber += 1
        else:
            self.node = nodeNumber
        ParticleNode._lastNodeNumber = max(self.node, ParticleNode._lastNodeNumber)
        self.canonName = None
        # Node weight:
        # (1 for stable particles, BR for unstable and production xsec for primary vertex)
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

    def __getattr__(self, attr):
        """
        Returns the attribute from particle.

        :parameter attr: Attribute string

        :return: Attribute from self.particle
        """

        return getattr(self.particle, attr)

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

    def cmpNode(self, other):
        """
        Compare nodes. Nodes are equal if their canonNames are equal
        and the particles match.

        :param other: node (ParticleNode object)
        """

        if not isinstance(other, ParticleNode):
            return -1
        if self.canonName == other.canonName:
            return self.cmpParticle(self, other)
        elif self.canonName > other.canonName:
            return 1
        else:
            return -1

    def eqNode(self, other):
        return (self.cmpNode(other) == 0)

    def cmpParticle(self, other):
        return self.particle.__cmp__(other.particle)

    def eqParticle(self, other):
        return self.particle == other.particle

    def isFinalState(self):
        """
        Checks if the node should be considered as a final state
        (should not be decayed). Returns True if .finalState has been
        defined and set to True
        """


def compareNodes(treeA, treeB, nodeA, nodeB):
    """
    Check if nodeA from treeA equals nodeB from treeB.
    The comparison is based on the node canonical name (node structure)
    and the node particles. It is done recursively, so two nodes are equal
    only if all their daughters are also equal.
    For daughters with the same canonical name, all the permutations are tried
    when comparing nodes.
    If nodes match return a tree (Tree obj) with nodeB as root
    and its daughters as nodes, but with the daughters sorted according to the ordering
    to which they matched nodeA.

    :param treeA: a tree (Tree object)
    :param treeB: a tree (Tree object)
    :param nodeA: a node belonging to treeA (ParticleNode object)
    :param nodeB: a node belonging to treeB (ParticleNode object)

    :return: (True, new tree) if nodes match, (False, None) otherwise.
    """

    # logger.debug('Comparing nodes', nodeA, nodeB)

    if not isinstance(nodeA, ParticleNode):
        return -1, None
    if not isinstance(nodeB, ParticleNode):
        return -1, None

    if nodeA.canonName != nodeB.canonName:
        # logger.debug('\t node names differ', nodeA, nodeB)
        if nodeA.canonName > nodeB.canonName:
            return 1, None
        else:
            return -1, None

    if nodeA.particle != nodeB.particle:
        # print('\t node particles differ', nodeA, nodeB)
        if nodeA.particle > nodeB.particle:
            return 1, None
        else:
            return -1, None

    daughtersA = list(treeA.successors(nodeA))
    daughtersB = list(treeB.successors(nodeB))
    if not daughtersA and not daughtersB:
        newNodeB = Tree()
        newNodeB.add_node(nodeB)
        return 0, newNodeB

    # Compute all permutations within each subgroup with a common canon name:
    allPerms = []
    for key, group in itertools.groupby(daughtersB, lambda d: d.canonName):
        allPerms.append(list(itertools.permutations(group)))

    # Construct all permutations and check against daughtersA:
    cmp = 0
    cmp_orig = None
    for dB_perm in itertools.product(*allPerms):
        dB_perm = list(itertools.chain.from_iterable(dB_perm))
        cmp = 0
        newDaughtersB = []
        for iB, dB in enumerate(dB_perm):
            cmp, newNodeB = compareNodes(treeA, treeB, daughtersA[iB], dB)
            newDaughtersB.append(newNodeB)
            if cmp != 0:  # If one node differs, try next permutation
                break

        # The first permutation is the original ordering.
        # Store the comparison for the original ordering in case it does not match,
        # since this will be the value returned in this case.
        if cmp_orig is None:
            cmp_orig = cmp

        if cmp == 0:  # If permutation matches, stop
            break

    # If matches also return a tree with the nodeB as root and its daughters
    # ordered according to how they matched nodeA
    if cmp == 0:
        newNodeB = Tree()
        newNodeB.add_node(nodeB)
        for dB in newDaughtersB:
            dB.add_node(nodeB)
            dB.add_edge(nodeB, list(dB.nodes())[0])
            newNodeB = nx.compose(newNodeB, dB)
        # logger.debug('\t returning cmp = ', 0)
        return 0, newNodeB
    else:
        # logger.debug('\t returning cmp = ', cmp_orig)
        return cmp_orig, None


class Tree(nx.DiGraph):
    """
    A wrapper for the tree graphs used to describe elements.
    """

    def __init__(self, info=None, finalState=None,
                 intermediateState=None, model=None):

        self.canonName = None
        if not info:  # Initialize empty tree
            nx.DiGraph.__init__(self)
        elif info and isinstance(info, str):
            nx.DiGraph.__init__(self)
            try:
                self.stringToTree(info, finalState, intermediateState, model)
            except (SModelSError, TypeError):
                raise SModelSError("Can not create element from input %s" % info)
        elif isinstance(info, dict):
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

        decays = procString.replace(" ,", ",").replace(", ", ",").split("),(")
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
            ptcs = [dec.split('>')[0].strip()] + [p.strip() for p in dec.split('>')[1].split(',')]

            for iptc, ptc in enumerate(ptcs):
                if ptc in nodesDict:
                    n = nodesDict[ptc]
                else:
                    n = maxNode+1
                    maxNode = n
                particleLabel = ptc.replace('(%i)' % n, '')
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
                if iptc == 0:
                    momNode = node
                self.add_node(node)
                if node != momNode:
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
        root = self.getTreeRoot()
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
            elStr = elStr[:-1]+'), '
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
        branches = [b for b in tree.successors(self.getTreeRoot())]
        finalState = []
        intermediateState = []
        branchList = []
        if len(branches) != 2:
            raise SModelSError("Can not convert tree to bracket with %i branches" % len(branches))
        for ib, b in enumerate(branches):
            intermediateState.append([])
            branchList.append([])
            for mom, daughters in nx.bfs_successors(tree, b):
                vertexList = [str(d) for d in daughters if d.Z2parity == 1]
                fstates = [str(d) for d in daughters if d.Z2parity == -1 and not list(tree.successors(d))]
                if len(vertexList) != len(daughters)-1 or len(fstates) > 1:
                    raise SModelSError("Can not convert tree with Z2-violating decays to bracket")
                branchList[ib].append(vertexList)
                intermediateState[ib].append(str(mom))
                finalState += fstates

        return branchList, finalState, intermediateState

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
            node = self.getTreeRoot()

        children = list(self.successors(node))
        if not children:
            node.canonName = 10
        else:
            tp = sorted([self.setCanonName(n) for n in children])
            tpStr = '1'+"".join(str(c) for c in tp)+'0'
            node.canonName = int(tpStr)
            # Use the root canon name to assign the tree canon name
            if self.in_degree[node] == 0:
                self.canonName = node.canonName

        return node.canonName

    def getFinalStates(self):
        """
        Get the list of particles which have not decayed (appear at the top of the tree)

        :returns: list of ParticleNode objects
        """

        finalStates = [n for n in self.nodes() if self.out_degree(n) == 0]
        return finalStates

    def getTreeRoot(self):
        """
        Get the root node (primary vertex) of the tree.

        :return: node
        """

        root = [n for n in self.nodes() if self.in_degree[n] == 0]
        if len(root) != 1:
            raise SModelSError("Malformed Tree, %i root(s) have been found." % len(root))

        return root[0]

    def getTreeWeight(self):
        """
        Computes the tree weight (production cross-section*BRs). If it can not be computed,
        return None. Does not include the weight of final state (undecayed) particles.

        :return: tree weight (Unum object) if available, None otherwise.
        """

        weight = 1.0
        for n in self.nodes():
            # Skip final state (stable) particles
            if self.out_degree(n) == 0:
                continue
            if not hasattr(n, 'nodeWeight') or n.nodeWeight is None:
                return None
            weight *= n.nodeWeight

        return weight

    def bfs_successors(self, node):
        """
        Returns an iterator over the mother and daughter
        nodes starting at node using a breadth first search.

        :param node: Node from tree

        :return: Iterator over nodes.
        """

        return nx.bfs_successors(self, node)

    def sort(self):
        """
        Sort the tree according to the canonName and then particle. For each node,
        all the daughters are sorted according to their canonName.

        :return: sorted tree (Tree object)
        """

        # If canon names have not been defined, compute them:
        if self.canonName is None:
            self.setCanonName()

        sortedTreeDict = {}
        for mom, daughters in nx.bfs_successors(self, self.getTreeRoot()):
            sortedDaughters = sorted(daughters, key=lambda d: (d.canonName, d.particle))
            sortedTreeDict[mom] = sortedDaughters

        # Remove all nodes and edges:
        self.remove_nodes_from(list(self.nodes()))
        # Replace with sorted nodes:
        for mom, daughters in sortedTreeDict.items():
            self.add_node(mom)
            for d in daughters:
                self.add_edge(mom, d)

    def copyTree(self):
        """
        Returns a shallow copy of self.

        :return: Tree object
        """

        newTree = self.copy()
        newTree.canonName = self.canonName

        return newTree

    def checkConsistency(self):
        """
        Make sure the tree has the correct topology (directed rooted tree).
        Raises an error otherwise.
        """
        if self.number_of_nodes() and not nx.is_arborescence(self):
            raise SModelSError("Graph created with malformed structure (not  a tree).")

    def addTrees(self, other):
        """
        Combines the nodes (add particles) in equivalent nodes in each trees.
        Can only be done if the trees have the same topology and ordering.

        :param other: tree (Tree object)

        :return: new tree with the combined particles (Tree object)
        """

        nodesDict = {n: n+list(other.nodes())[inode]
                     for inode, n in enumerate(self.nodes())}
        newTree = nx.relabel_nodes(self, nodesDict, copy=True)

        return newTree

    def compose(self, other):
        """
        Returns a new graph of self composed with other.
        The trees are joined where the nodes overlap and the nodes
        from self are kept.

        :param other: tree (Tree object)

        :return: new tree with the other composed with self.
        """

        return nx.compose(self, other)

    def draw(self, particleColor='lightcoral',
             smColor='skyblue',
             pvColor='darkgray',
             nodeScale=4, labelAttr=None, attrUnit=None):
        """
        Draws Tree using matplotlib.

        :param particleColor: color for particle nodes
        :param smColor: color used for particles which have the _isSM attribute set to True
        :param pvColor: color for primary vertex
        :param nodeScale: scale size for nodes
        :param labelAttr: attribute to be used as label. If None, will use the string representation
                          of the node object.
        :param attrUnit: Unum object with the unit to be removed from label attribute (if applicable)

        """
        import matplotlib.pyplot as plt

        if labelAttr is None:
            labels = {n: str(n) for n in self.nodes()}
        elif attrUnit is not None:
            labels = {n: str(getattr(n, labelAttr).asNumber(attrUnit)) if hasattr(n, labelAttr)
                      else str(n) for n in self.nodes()}
        else:
            labels = {n: str(getattr(n, labelAttr)) if hasattr(n, labelAttr)
                      else str(n) for n in self.nodes()}

        for key in labels:
            if labels[key] == 'anyOdd':
                labels[key] = 'BSM'
        node_size = []
        node_color = []
        for n in self.nodes():
            node_size.append(nodeScale*100*len(labels[n]))
            if 'pv' == labels[n].lower():
                node_color.append(pvColor)
            elif hasattr(n, '_isSM') and n._isSM:
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
