"""
.. module:: graphTools
   :synopsis: A collection of functions used to create and maninuplate graphs. Based on networkx.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


import networkx as nx
from smodels.theory.exceptions import SModelSTheoryError as SModelSError


class ParticleNode(object):
    """
    Simple wrapper for creating graphs with Particle objects.
    It is necessary because the same particle can appear multiple times within a tree, so
    Particle objects can not be directly used as nodes
    (since the same particle can not appear as distinct nodes)
    """

    def __init__(self, particle, nodeNumber):
        self.particle = particle
        self.node = nodeNumber
        self.canonName = None

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


def stringToTree(stringElement, model=None, finalState=None, intermediateState=None):
    """
    Converts a string describing an element to a Tree (DiGraph object). It can
    be using the (old) bracket notation or the process notation.

    :param procString: The process in string format  (e.g. '(PV > gluino(1),squark(2)), (gluino(1) >
                       MET,jet,jet), (squark(2) > HSCP,u)' or [[['jet','jet']],[['u']]]). The particle
                       labels should match the particles in the Model (if Model != None).

    :parameter model: The model (Model object) to be used when converting particle labels to
                          particle objects. If None, the nodes will only store the particle labels.

    :parameter finalState: (optional) list containing the final state labels for each branch
                           (e.g. ['MET', 'HSCP'] or ['MET','MET'])
    :parameter intermediateState: (optional) nested list containing intermediate state labels
                                     for each branch  (e.g. [['gluino'], ['gluino']])

    :return: tree describing the process (DiGraph object)
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

    # Find maximum node defined so far:
    maxNode = max([n for n in nodesDict.values() if n is not None])

    # Now build the Graph:
    T = nx.DiGraph()
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
            T.add_node(node)
            if node != momNode:
                edges.append((momNode, node))
    # Now add all edges
    T.add_edges_from(edges)

    return T


def treeToString(tree):
    """
    Convert a tree (DiGraph object) to a process string (e.g. '(PV > gluino(1),squark(2)), (gluino(1) >
                       MET,jet,jet), (squark(2) > HSCP,u)')

    :param T: Tree (DiGraph object)

    :return: String describing the process
    """

    elStr = ""
    T = tree
    root = getTreeRoot(T)
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


def bracketToProcessStr(stringEl, finalState=None, intermediateState=None):
    """
    :parameter info: string describing the element in bracket notation
                         (e.g. [[[e+],[jet]],[[e-],[jet]]])

    :parameter finalState: list containing the final state labels for each branch
                         (e.g. ['MET', 'HSCP'] or ['MET','MET'])
    :parameter intermediateState: nested list containing intermediate state labels
                                     for each branch  (e.g. [['gluino'], ['gluino']])
    """

    branches = eval(stringEl.replace('[*]', "['InclusiveBranch']"))

    if finalState is None:
        fStates = ['MET']*len(branches)
    else:
        fStates = finalState[:]

    if intermediateState is None:
        intStates = [['anyOdd' for dec in br] for br in branches]
    else:
        intStates = intermediateState[:]

    stringProc = '(PV > ' + ','.join(['%s(%i)' % (intStates[i][0], i+1) for i, br in enumerate(branches)]) + ')'

    j = len(branches)
    for ibr, br in enumerate(branches):
        fs = fStates[ibr]
        brStr = ''
        i = ibr+1
        for idec, dec in enumerate(br):
            mom = intStates[ibr][idec]
            if idec < len(br)-1:
                daughter = intStates[ibr][idec+1]
                j += 1
                brStr += ', (%s(%i) > %s(%i),' % (mom, i, daughter, j) + ','.join(dec) + ')'
            else:
                brStr += ', (%s(%i) > %s,' % (mom, i, fs) + ','.join(dec) + ')'
            i = j

        stringProc += brStr

    return stringProc


def fromTreeToList(T, node=None):
    """
    Convert a Tree to a nested list with the Z2 even
    final states. The Z2 odd final states (e.g. 'MET', 'HSCP') are
    not included.

    :param T: DiGraph object
    :param node: Name of node. If None it will start with the primary node

    :return: Nested list with Z2-even particle objects (e.g. [[[e-,mu],[L]],[[jet]]])
    """

    if not isinstance(T, nx.DiGraph):
        raise SModelSError("Input must be a DiGraph object.")

    # Get the primary vertex:
    if node is None:
        node = getTreeRoot(T)

    dList = []
    children = list(T[node])
    if children:
        for n in children:
            if list(T[n]):
                dList.append(fromTreeToList(T, n))
            else:
                if T.nodes[n].Z2parity == 'even':
                    dList.append(T.nodes[n])

    return dList


def getCanonName(T, node=None):
    """
    Recursively compute the canonical name for the Tree T.
    Returns the name in integer form.

    :param T: Tree (networkX DiGraph object)
    :param node: Node to get the name for. If None, it will use the root

    :return: Integer representing the Tree
    """

    if not isinstance(T, nx.DiGraph):
        raise SModelSError("Input must be a DiGraph object.")

    if node is None:
        node = getTreeRoot(T)

    children = list(T.successors(node))
    if not children:
        node.canonName = 10
    else:
        tp = sorted([getCanonName(T, n) for n in children])
        tpStr = '1'+"".join(str(c) for c in tp)+'0'
        node.canonName = int(tpStr)

    return node.canonName


def getNodeLevelDict(T):
    """
    Return a dictionary with the nodes in each level of the tree T
    (e.g. {0 : [PV], 1 : [A,B],...})

    :param T: DiGraph object
    :return: Dictionary with the levels as keys and a list of nodes in each level as values.
    """

    if not isinstance(T, nx.DiGraph):
        raise SModelSError("Input must be a DiGraph object.")

    root = getTreeRoot(T)
    levelNodes = {}
    for n in T.nodes():
        d = nx.shortest_path_length(T, root, n)  # Distance to root
        if d not in levelNodes:
            levelNodes[d] = [n]
        else:
            levelNodes[d].append(n)

    return levelNodes


def getTreeRoot(T):
    """
    Get the root node (primary vertex) of a tree T.

    :param T: DiGraph object
    :return: node
    """

    if not isinstance(T, nx.DiGraph):
        raise SModelSError("Input must be a DiGraph object.")

    root = [n for n in T.nodes() if T.in_degree[n] == 0]
    if len(root) != 1:
        raise SModelSError("Malformed Tree, %i root(s) have been found." % len(root))

    return root[0]


def drawTree(tree, oddColor='lightcoral', evenColor='skyblue',
             pvColor='darkgray', genericColor='lightgray',
             nodeScale=4, labelAttr=None):
    """
    Draws Tree using matplotlib.

    :param tree: tree to be drawn
    :param oddColor: color for Z2-odd particles
    :param evenColor: color for Z2-even particles
    :param genericColor: color for particles without a defined Z2-parity
    :param pvColor: color for primary vertex
    :param nodeScale: scale size for nodes
    :param labelAttr: attribute to be used as label. If None, will use the string representation
                      of the node object.

    """
    import matplotlib.pyplot as plt

    if labelAttr is None:
        labels = {n: str(n) for n in tree.nodes()}
    else:
        labels = {n: str(getattr(n, labelAttr)) if hasattr(n, labelAttr) else str(n) for n in tree.nodes()}

    for key in labels:
        if labels[key] == 'anyOdd':
            labels[key] = 'BSM'
    node_size = []
    node_color = []
    for n in tree.nodes():
        node_size.append(nodeScale*100*len(labels[n]))
        if 'pv' == labels[n].lower():
            node_color.append(pvColor)
        else:
            node_color.append(genericColor)

    # Compute position of nodes (in case the nodes have the same string representation, first
    # convert the nodes to integers for computing the position)
    H = nx.convert_node_labels_to_integers(tree, label_attribute='node_label')
    H_layout = nx.drawing.nx_agraph.graphviz_layout(H, prog='dot')
    pos = {H.nodes[n]['node_label']: p for n, p in H_layout.items()}

    nx.draw(tree, pos,
            with_labels=True,
            arrows=True,
            labels=labels,
            node_size=node_size,
            node_color=node_color)

    plt.show()
