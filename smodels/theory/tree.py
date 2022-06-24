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
from collections import OrderedDict


class Tree(object):
    """
    A wrapper for the tree graphs used to describe elements.
    """

    def __init__(self, info=None, finalState=None,
                 intermediateState=None, model=None):

        self._canonName = None
        self._root = None
        self.nodes_and_edges = OrderedDict()

        # Convert from string:
        if isinstance(info, str):
            try:
                self.stringToTree(info, finalState, intermediateState, model)
            except (SModelSError, TypeError):
                raise SModelSError("Can not create element from input %s" % info)
        # Convert from dictionary with edges ({node : [node1,node2,..], ...}):
        elif isinstance(info, dict):
            for node1, nodeList in info.items():
                self.add_edges_from(zip([node1]*len(nodeList), nodeList))
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

        if node not in self.nodes_and_edges:
            self.nodes_and_edges[node] = []

    def add_nodes_from(self, nodes):

        for node in nodes:
            self.add_node(node)

    def add_edge(self, edge):

        nodeA = edge[0]
        nodeB = edge[1]
        if nodeA not in self.nodes_and_edges:
            self.nodes_and_edges[nodeA] = [nodeB]
        elif nodeB not in self.nodes_and_edges[nodeA]:
            self.nodes_and_edges[nodeA].append(nodeB)

    def add_edges_from(self, edges):

        for edge in edges:
            self.add_edge(edge)

    @property
    def nodes(self):

        for node in self.nodes_and_edges.keys():
            yield node

    @property
    def edges(self):

        for node, nodeList in self.nodes_and_edges.items():
            for nodeB in nodeList:
                yield (node, nodeB)

    def out_degree(self, node):

        if node not in self.nodes_and_edges:
            return 0
        else:
            return len(self.nodes_and_edges[node])

    def in_dregree(self, node):

        degree = 0
        # Count how many times node appears as a daughter
        for nodeList in self.nodes_and_edges.values():
            degree += nodeList.count(node)

        return degree

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
                elif self.out_degree(n) == 0:
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

        branches = [b for b in self.successors(self.root)]
        finalState = []
        intermediateState = []
        branchList = []
        if len(branches) != 2:
            raise SModelSError("Can not convert tree to bracket with %i branches" % len(branches))
        for ib, b in enumerate(branches):
            intermediateState.append([])
            branchList.append([])
            # Deal separately with the case where the primary mother is stable:
            if self.out_degree(b) == 0:
                if not b.isSM:
                    finalState.append(str(b))
                    continue
                else:
                    raise SModelSError("Can not convert tree with Z2-violating decays to bracket")
            for mom, daughters in self.bfs_successors(b):
                vertexList = [str(d) for d in daughters if d.isSM]
                fstates = [str(d) for d in daughters if not d.isSM and self.out_degree(d) == 0]
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

    def setGlobalProperties(self, node=None):
        """
        Define useful global properties for the trees and the nodes.
        Recursively sets the canonName and finalStates for each node.
        Returns the canonical name in integer form.

        :param T: Tree (Tree object)
        :param node: Node to get the name for. If None, it will use the root

        :return: Integer representing the Tree canonical name
        """

        if not self.number_of_nodes():
            return None

        if node is None:
            node = self.root

        # Set the final state
        node.finalStates = self.getFinalStates(node)
        # If it is inclusive node set its name to an inclusive integer
        # and return its name (no need to check the children)
        if node.isInclusive:
            node.canonName = InclusiveValue()
            return node.canonName

        children = list(self.successors(node))
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
        if self.in_degree[node] == 0:
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
            finalStates = [n.particle for n in self.nodes if self.out_degree(n) == 0]
            return sorted(finalStates)

        # Return source, if it is already a final state
        if self.out_degree(source) == 0:
            return [source.particle]

        finalStates = []
        for mom, daughters in self.bfs_successors(source):
            finalStates += [d.particle for d in daughters if self.out_degree(d) == 0]

        return sorted(finalStates)

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

    def getNode(self, node):
        """
        Return the node object in the tree corresponding to node.
        Node can be an integer or a node object. The node.node number
        is used to retrive the object. If the node is not found, return None.

        :param node: ParticleNode, InclusiveParticleNode or int

        :return: Node object if found, otherwise None.
        """

        # Allows comparison with integers (by node.node) or node objs
        if isinstance(node, int):
            nodeNumber = node
        elif isinstance(node, (ParticleNode, InclusiveParticleNode)):
            nodeNumber = node.node
        else:
            raise SModelSError("Can not fetch node object from type %s" % (type(node)))

        for n in self.nodes:
            if n.node == nodeNumber:
                return n

        return None

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

    def bfs_successors(self, node=None, includeLeaves=False):
        """
        Returns an iterator over the mother and daughter
        nodes starting at node using a breadth first search.

        :param node: Node from tree. If None, starts at tree root.
        :param includeLeaves: If True, it will consider the leaves (undecayed nodes)
                              as moms in the iterator (with an empty daughters list)

        :return: Iterator over nodes.
        """

        if node is None:
            node = self.root

        mom = node
        daughters = self.nodes_and_edges[mom]
        generation = [(mom, daughters)]
        while generation:
            for pair in generation:
                yield pair
            next_generation = []
            for pair in generation:
                mom, daughters = pair
                for new_mom in daughters:
                    new_daughters = self.nodes_and_edges[new_mom]
                    if new_daughters or includeLeaves:
                        next_generation.append((new_mom, new_daughters))
            generation = next_generation

    def dfs_successors(self, node=None, includeLeaves=False):
        """
        Returns a dictionary with the mother as keys
        and the daughters as values using a depth-first search.
        nodes starting at node using a breadth first search.

        :param node: Node from tree. If None, starts at tree root.
        :param includeLeaves: If True, it will consider the leaves (undecayed nodes)
                              as moms in the iterator (with an empty daughters list)

        :return: Dictionary with mothers and daughters
        """

        if node is None:
            node = self.root

        mom = node
        daughters = self.nodes_and_edges[node][:]
        daughtersDict = {mom: daughters[:]}

        # Go down in generation until it has no more daughters
        while daughters:
            new_mom = daughters.pop(0)
            new_daughters = self.nodes_and_edges[new_mom][:]
            if new_daughters or includeLeaves:
                daughtersDict.update({new_mom: new_daughters[:]})
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
        T.add_node(node)
        for mom, daughters in self.bfs_successors(node):
            T.add_node(mom)
            for daughter in daughters:
                T.add_edge((mom, daughter))
        return T

    def getSubTree(self, source):
        """
        Returns the subtree with source as its root.

        :param source: Node
        :return: Tree object with source as its root node
        """

        subTree = Tree()
        subTree.add_node(source)
        for mom, daughters in self.bfs_successors(source):
            for d in daughters:
                subTree.add_edge(mom, d)

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
        for daughter in list(self.successors(self.root))[:]:
            subTree = self.getSubTree(source=daughter)
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

        # Get
        nodeOrder = []
        for mom, daughters in self.bfs_successors(includeLeaves=True):
            nodeOrder.append((mom, daughters))

        # Re-number nodes:
        inode = 0
        for node, daughters in nodeOrder:
            node.node = inode
            inode += 1

        new_nodes_and_edges = OrderedDict()
        for node, daughters in nodeOrder:
            daughters = sorted(daughters, key=lambda d: d.node)
            new_nodes_and_edges[node] = daughters

        self.nodes_and_edges = new_nodes_and_edges

    def copyTree(self):
        """
        Returns a copy of self. The copy contains the same canonical name
        and copies of nodes.

        : return: Tree object
        """

        newNodesDict = {n: n.copy() for n in self.nodes}
        newTree = self.relabel_nodes(newNodesDict, copy=True)
        newTree._canonName = self._canonName

        return newTree

    def relabel_nodes(self, newNodesDict, copy=True):
        """
        Replace the nodes in the Tree according to newNodesDict.
        If copy = True, create a new Tree object, otherwise replace inplace.

        :param newNodesDict: A dictionary with current nodes as keys as new nodes as values.
        :param copy: If True, create a new Tree object, otherwise replace inplace.

        :return: A new tree if copy = True, otherwise returns self
        """

        nodes_and_edges = list(self.nodes_and_edges.list())
        new_nodes_and_edges = OrderedDict()
        for entry in nodes_and_edges:
            node = entry[0]
            daughters = entry[1][:]
            new_node = newNodesDict[node]
            new_daughters = [newNodesDict[d] for d in daughters]
            new_nodes_and_edges[new_node] = new_daughters

        if copy:
            newTree = self.copyTree()
        else:
            newTree = self

        newTree.nodes_and_edges = new_nodes_and_edges

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
        fsNodes = [node for node in newTree.nodes if newTree.out_degree(node) == 0]
        fsNodes = sorted(fsNodes, key=lambda node: node.particle)
        # Remove all nodes
        newTree.remove_nodes_from(list(newTree.nodes))

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
        if len(self.nodes_and_edges) > 1:
            root = self.root

            # Check if root has no parents and at least one daughter
            if self.in_dregree(root) != 0 or self.out_degree(root) == 0:
                malformedTree = True

            # Check if all nodes (except root) have a unique parent
            if any(self.in_dregree(node) != 1 for node in self.nodes if node is not root):
                malformedTree = True

            # Check if all nodes can be reached from the root node
            if len(self.dfs_successors(includeLeaves=True)) != len(self.nodes_and_edges):
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

        nodesB = list(other.nodes)
        nodesDict = {n: n + nodesB[inode]
                   for inode, n in enumerate(self.nodes)}
        newTree = self.relabel_nodes(nodesDict, copy=True)

        return newTree

    def compose(self, other, atNode=None, copy=True):
        """
        Returns a new graph of self composed with other at node atNode.
        If atNode is not defined, the trees are joined where the nodes overlap and the nodes
        from self are kept.

        :param other: tree (Tree object)
        :param atNode: node (ParticleNode object) where the merge has to take place
                        If not defined, will find the common node between the trees.
        :param copy: If True, return a copy of self, merged with other. Otherwise,
                     extend self.

        : return: new tree with the other composed with self.
        """

        if atNode is None:
            # Find intersection:
            self_nodes = set(list(self.nodes))
            other_nodes = set(list(other.nodes))
            common_nodes = self_nodes.intersection(other_nodes)
            if len(common_nodes) != 1:
                raise SModelSError("Can not merge trees. %i common nodes found" % len(common_nodes))

            atNode = common_nodes[0]  # merge node of self

        if self.out_degree(atNode) != 0:
            raise SModelSError("Can not merge trees. Common node in Tree A can not have daughters.")
        if other.in_degree(atNode) != 0:
            raise SModelSError("Can not merge trees. Common node in Tree B can not have parents.")

        if copy:
            newTree = self.copyTree()
        else:
            newTree = self

        for mom, daughters in other.bfs_successors(atNode):
            if mom == atNode:
                mom = atNode  # Make sure the merge node of self is kept
            newTree.nodes_and_edges[mom] = daughters[:]

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
