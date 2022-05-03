"""
.. module:: graphTools
   :synopsis: A collection of functions used to create and maninuplate graphs. Based on networkx.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


import networkx as nx
from smodels.theory.exceptions import SModelSTheoryError as SModelSError


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

    # Build a dictionary with all parent nodes:
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
        momNode = nodesDict[ptcs[0]]
        for ptc in ptcs:
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
            T.add_node(n, particle=particle)
            if n != momNode:
                edges.append((momNode, n))
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
        if mom == 0:
            elStr += '(%s > ' % T.nodes[mom]['particle']
        else:
            if mom not in nodesDict:
                nodesDict[mom] = counter
                counter += 1
            elStr += '(%s(%i) > ' % (T.nodes[mom]['particle'], nodesDict[mom])

        # Add daughters for the decay, if daughter is unstable, add index
        for n in daughters:
            if not list(T.successors(n)):
                elStr += '%s,' % (T.nodes[n]['particle'])  # stable
            else:
                if n not in nodesDict:
                    nodesDict[n] = counter
                    counter += 1
                elStr += '%s(%i),' % (T.nodes[n]['particle'], nodesDict[n])
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
                if T.nodes[n]['particle'].Z2parity == 'even':
                    dList.append(T.nodes[n]['particle'])

    return dList


def getTopologyName(T, node=None):
    """
    Recursively compute the topology (canonical) name for the Tree T.
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
        T.nodes[node]['canonName'] = 10
    else:
        tp = sorted([getTopologyName(T, n) for n in children])
        tpStr = '1'+"".join(str(c) for c in tp)+'0'
        T.nodes[node]['canonName'] = int(tpStr)

    return T.nodes[node]['canonName']


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


def drawTree(T, oddColor='lightcoral', evenColor='skyblue',
            pvColor='darkgray', genericColor='lightgray',
            nodeScale=4):
    """
    Draws Tree using matplotlib.
    If outputFile is defined, it will save plot to this file.
    """
    import matplotlib.pyplot as plt

    labels = dict([[n, str(T.nodes[n]['particle'])]
                 for n in T.nodes()])
    for key in labels:
        if labels[key] == 'anyOdd':
            labels[key] = 'BSM'
    node_size = []
    node_color = []
    for n in T.nodes():
        node_size.append(nodeScale*100*len(labels[n]))
        if 'pv' == labels[n].lower():
            node_color.append(pvColor)
        else:
            node_color.append(genericColor)
    pos = nx.drawing.nx_agraph.graphviz_layout(T, prog='dot')
    nx.draw(T, pos,
            with_labels=True,
            arrows=True,
            labels=labels,
            node_size=node_size,
            node_color=node_color)

    plt.show()
