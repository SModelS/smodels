"""
.. module:: genericSMS
   :synopsis: This is a base class for describing Simplified Model Topologies  using a rooted tree syntax.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


from smodels.base.exceptions import SModelSBaseError as SModelSError
from collections import OrderedDict
from itertools import product


class GenericGraph(object):

    """
    A generic class for describing and manipulating 
    simple graphs.
    """

    def __init__(self):
        """
        Initialize basic attributes.
        """

        self._successors = OrderedDict()  # Stores the nodes and their successors (daughters)
        self._predecessors = {}  # Stores the nodes and their predecessors (parents)
        self._nodesMapping = {}  # Stores the nodeIndex->node object mapping

    def __hash__(self):
        return object.__hash__(self)

    def __repr__(self):
        """
        Returns the string representation of the graph.
        """

        return str(self)
    
    def __str__(self):
        """
        Returns a string listing the graph nodes.
        """

        gStr = ",".join([f'({self.indexToNode(nI)} > {str([d for d in self.daughters(nI)])})' 
                         for nI in self.nodeIndices])

        return gStr


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
            val = [getattr(node, attr) for node in self.nodes]
            return val
        except AttributeError:
            raise AttributeError("Neither Graph nor its nodes have attribute ``%s''" % attr)
        
    def add_node(self, node, nodeIndex=None):
        """
        Adds a node object to the graph. If nodeIndex is None,
        the node index will be automatically assigned.

        :param node: node object
        :param nodeIndex: The index for the node. It must not
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
        Adds a list of nodes to the Graph.

        :param nodes: List of node objects

        :return: A list of node indices for the newly added nodes
        """

        nodeIndices = []
        for node in nodes:
            nodeIndices.append(self.add_node(node))

        return nodeIndices
    
    def remove_node(self, nodeIndex):
        """
        Removes a node from the Graph if the nodeIndex is in the Graph.
        The node is removed as well as its appearence in any edges.

        :param nodeIndex: Node index
        """

        if nodeIndex in self._successors:
            self._successors.pop(nodeIndex)
            self._nodesMapping.pop(nodeIndex)

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
        Removes a list of nodes from the Graph.

        :param nodeIndices: List of node indices
        """

        for nodeIndex in nodeIndices:
            self.remove_node(nodeIndex)

    def add_edge(self, nodeIndexA, nodeIndexB):
        """
        Adds a directed edge to existing nodes in the Graph (nodeA -> nodeB).

        :param nodeIndexA: Index for node A
        :param nodeIndexB: Index for node B
        """

        self._successors[nodeIndexA].append(nodeIndexB)
        self._predecessors[nodeIndexB] = nodeIndexA

    def add_edges_from(self, edges):
        """
        Adds a list of directed edges to the Graph.

        :param edges: List of tuples containing node indices
                      (e.g. [(nodeIndexA,nodeIndexB),(nodeIndexA,nodeIndexC),...])
        """

        for edge in edges:
            self.add_edge(edge[0], edge[1])

    def remove_edge(self, nodeIndexA, nodeIndexB):
        """
        Removes an edge from the tree if the edge
        (nodeIndexA -> nodeIndexB) is in the Graph.

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
        Removes edges from the tree if they appear in the Graph.

        :param edges: List of tuples containing node indices
                      (e.g. [(nodeIndexA,nodeIndexB),(nodeIndexA,nodeIndexC),...])
        """
        for edge in edges:
            self.remove_edge(edge[0], edge[1])

    def clear(self):
        """
        Remove all nodes and edges from the graph.
        """

        self._successors = OrderedDict()
        self._predecessors = {}
        self._nodesMapping = {}

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
            raise SModelSError("Can not convert object of type %s to nodes" %str(type(nodeIndex)))

    def daughterIndices(self, nodeIndex):
        """
        Returns the list of node indices corresponding to the
        successors of nodeIndex.

        :param nodeIndex: Parent node index
        """

        successors = self._successors[nodeIndex]

        return successors

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
            
        # Update dicts:
        self._successors = newSuccessors
        self._predecessors = newPredecessors
        self._nodesMapping = newMapping

    def updateNodeObjects(self, nodeObjectDict):
        """
        Update the node index -> node object mapping.
        Only affects the indices appearing in nodeObjectDict.

        :param nodeObjectDict: Dictionary with current node indices as keys
                              and new node objects as values
        """

        for nodeIndex, newObj in nodeObjectDict.items():
            self._nodesMapping[nodeIndex] = newObj

    def draw(self, particleColor='steelblue2',
                smColor='lightpink2',
                pvColor='darkgray',
                labelAttr='label',
                attrUnit=None, filename=None, view=True,
                maxLabelSize=10,
                graph_kwargs={'layout' : 'dot', 'ranksep' : '0.3', 'rankdir' : "LR"},
                nodes_kwargs={'style' : 'filled', 'fontsize' : '10', 'color' : 'black','shape' : 'circle','margin' : '0'},
                edges_kwargs={'arrowhead' : 'vee', 'arrowsize' : '0.7',  'color' : 'grey53'}):
        """
        Draws Tree using matplotlib.

        :param particleColor: color for particle nodes
        :param smColor: color used for particles which have the isSM attribute set to True
        :param fontsize: Font size for labels
        :param labelAttr: attribute to be used as label. If None, will use the string representation of the node object. It can also be a dictionary with node indices as keys and the label strings as values.
        :param attrUnit: Unum object with the unit to be removed from label attribute(if applicable)
        :param filename: Filename to save drawing to.
        :param view: open a viewer after plotting
        :param maxLabelSize: Maximum size for the label string for the node. If the label is larger, it will be truncated.
                             If None/False/0, it will keep the full label.
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
        else:
            labels = {nodeIndex: str(getattr(n, labelAttr)) if hasattr(n, labelAttr)
                      else str(n) for n,nodeIndex in nodesAndIndices}

        for key in labels:
            if labels[key] == 'anyOdd':
                labels[key] = 'BSM'

        node_color = {}
        for n in self.nodeIndices:
            node = self.indexToNode(n)
            if hasattr(node, 'isSM') and node.isSM:
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
