"""
.. module:: graphTools
   :synopsis: A collection of functions used to create and maninuplate graphs. Based on networkx.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


import networkx as nx
import itertools
from smodels.theory.exceptions import SModelSTheoryError as SModelSError


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
    Check if nodeA from treeA equlas nodeB from treeB.
    The comparison is based on the node canonical name (node structure)
    and the node particles. It is done recursively, so two nodes are equal
    only if all their daughters are also equal.
    For daughters with the same canonical name, all the permutations are tried
    when comparing nodes.
    If nodes match return a tree (DiGraph obj) with nodeB as root
    and its daughters as nodes, but with the daughters sorted according to the ordering
    to which they matched nodeA.

    :param treeA: a tree (DiGraph object)
    :param treeB: a tree (DiGraph object)
    :param nodeA: a node belonging to treeA (ParticleNode object)
    :param nodeB: a node belonging to treeB (ParticleNode object)

    :return: (True, new tree) if nodes match, (False, None) otherwise.
    """

    if not isinstance(nodeA, ParticleNode):
        return -1, None
    if not isinstance(nodeB, ParticleNode):
        return -1, None

    if nodeA.canonName != nodeB.canonName:
        if nodeA.canonName > nodeB.canonName:
            return 1, None
        else:
            return -1, None
    # print('\n\nComparing %s with %s' %(nodeA.particle,nodeB.particle))
    # print('\t %s = %s (%s)' %(nodeA.particle,nodeB.particle,nodeA.particle == nodeB.particle))
    if nodeA.particle != nodeB.particle:
        if nodeA.particle > nodeB.particle:
            return 1, None
        else:
            return -1, None

    daughtersA = list(treeA.successors(nodeA))
    daughtersB = list(treeB.successors(nodeB))
    if not daughtersA and not daughtersB:
        newNodeB = nx.DiGraph()
        newNodeB.add_node(nodeB)
        return 0, newNodeB

    #  print('daughters A =',daughtersA)
    # print('daughters B =',daughtersB)

    # Compute all permutations within each subgroup with a common canon name:
    allPerms = []
    for key, group in itertools.groupby(daughtersB, lambda d: d.canonName):
        allPerms.append(list(itertools.permutations(group)))

    # Construct all permutations and check against daughtersA:
    cmp = 0
    cmp_orig = None
    for dB_perm in itertools.product(*allPerms):
        dB_perm = list(itertools.chain.from_iterable(dB_perm))
        # print('Checking permutation:',dB_perm)
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
        newNodeB = nx.DiGraph()
        newNodeB.add_node(nodeB)
        for dB in newDaughtersB:
            dB.add_node(nodeB)
            dB.add_edge(nodeB, list(dB.nodes())[0])
            newNodeB = nx.compose(newNodeB, dB)
        return 0, newNodeB
    else:
        return cmp_orig, None


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

    # Make sure particles have unique nodes:
    if len(set(list(nodesDict.values()))) != len(list(nodesDict.values())):
        raise SModelSError("Input string has non unique nodes: %s" % nodesDict)

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


def treeToBrackets(tree):
    """
    Convert a Tree to a nested list with the Z2 even
    final states. The Z2 odd final states (e.g. 'MET', 'HSCP') are
    not included. The Tree must be Z2-preserving and represent the
    pair production cascade decay of two Z2-odd particles.

    :param tree: DiGraph object

    :return: Nested list with the strings for the Z2-even final states
             (e.g. [[['e-','mu'],['L']],[['jet']]])
    """

    if not isinstance(tree, nx.DiGraph):
        raise SModelSError("Input must be a DiGraph object.")

    branches = [b for b in tree.successors(getTreeRoot(tree))]
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


def getMaxNode(tree):
    """
    Returns the maximum node number (ParticleNode.node) in the nodes from the tree.
    It can be relevant when adding new ParticleNodes, since nodes with same number will
    be considered equal within the tree and would not be added.

    :param tree: Tree (DiGraph object)

    :return: Maximum node number (int)
    """

    maxNode = max([n.node for n in tree.nodes()])
    return maxNode


def getFinalStates(tree):
    """
    Get the list of particles which have not decayed (appear at the top of the tree)

    :returns: list of ParticleNode objects
    """

    finalStates = [n for n in tree.nodes() if tree.out_degree(n) == 0]
    return finalStates


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


def getTreeWeight(tree):
    """
    Computes the tree weight (production cross-section*BRs). If it can not be computed,
    return None. Does not include the weight of final state (undecayed) particles.

    :param: tree (DiGraph object)

    :return: tree weight (Unum object) if available, None otherwise.
    """

    weight = 1.0
    for n in tree.nodes():
        # Skip final state (stable) particles
        if tree.out_degree(n) == 0:
            continue
        if not hasattr(n, 'nodeWeight') or n.nodeWeight is None:
            return None
        weight *= n.nodeWeight

    return weight


def sortTree(tree):
    """
    Sort the tree according to the canonName and then particle. For each node,
    all the daughters are sorted according to their canonName.

    :param tree: tree (DiGraph object)

    :return: sorted tree (DiGraph object)
    """

    sortedTree = nx.DiGraph()
    for mom, daughters in nx.bfs_successors(tree, getTreeRoot(tree)):
        sortedTree.add_node(mom)
        sortedDaughters = sorted(daughters, key=lambda d: (d.canonName, d.particle))
        for d in sortedDaughters:
            sortedTree.add_node(d)
            sortedTree.add_edge(mom, d)

    return sortedTree


def addTrees(treeA, treeB):
    """
    Combines the particles in equivalent nodes in each trees.
    Can only be done if the trees have the same topology and ordering.

    :param treeA: tree (DiGraph object)
    :param treeB: tree (DiGraph object)

    :return: new tree with the combined particles (DiGraph object)
    """

    nodesDict = {n: n+list(treeB.nodes())[inode]
                 for inode, n in enumerate(treeA.nodes())}
    newTree = nx.relabel_nodes(treeA, nodesDict, copy=True)

    return newTree


def getDecayTrees(mother):
    """
    Generates a simple list of trees with all the decay channels
    for the mother. In each tree the mother appears as the root
    and each of its decays as daughters.
    The  mother node weight is set to the respective decay branching ratio.
    (The node numbering for the root/mother node is kept equal,
    while the numbering of the daughters is automatically assigned to
    avoid overlap with any previously created nodes, so the
    decay tree can be directly merged to any other tree.)


    :param mother: Mother for which the decay trees will be generated (ParticleNode object)

    :return: List of simple trees representing the decay channels of daughter sorted in reverse
             order according to the BR (largest BR first)
    """

    decayTrees = []

    # Sort decays:
    decays = []
    for decay in mother.decays:
        if decay is not None:
            decays.append(decay)
        else:
            # Include possibility of mother appearing as a final state
            mom = mother.copy()
            mom.finalState = True  # Forbids further node decays
            decayTrees.append(nx.DiGraph({mom: []}))

    decays = sorted(decays, key=lambda dec: dec.br, reverse=True)

    # Loop over decays of the daughter
    for decay in decays:
        if not decay.br:
            continue  # Skip decays with zero BRs
        daughters = []
        mom = mother.copy()
        mom.nodeWeight = decay.br
        for ptc in decay.daughters:
            ptcNode = ParticleNode(particle=ptc)
            daughters.append(ptcNode)

        decayTrees.append(nx.DiGraph({mom: daughters}))

    return decayTrees


def addOneStepDecays(tree, sigmacut=None):
    """
    Given a tree, generates a list of new trees (DiGraph objects),
    where all the (unstable) nodes appearing at the end of the original tree
    have been decayed. Each entry in the list corresponds to a different combination
    of decays. If no decays were possible, return an empty list.

    :param tree: Tree (DiGraph object) for which to add the decays
    :param sigmacut: Cut on the tree weight (xsec*BR). Any tree with weights
                     smaller than sigmacut will be ignored.


    :return: List of trees with all possible 1-step decays added.
    """

    treeList = [tree]
    # Get all (current) final states which are the mothers
    # of the decays to be added:
    mothers = [n for n in tree.nodes() if tree.out_degree(n) == 0]
    for mom in mothers:
        # Check if mom should decay:
        if mom.finalState:
            continue
        if mom.particle.isStable():
            mom.finalState = True
            continue  # Skip if particle is stable
        # Skip if particle has no decays
        if not hasattr(mom, 'decays'):
            mom.finalState = True
            continue
        if not mom.decays:
            mom.finalState = True
            continue

        # Get all decay trees for final state:
        decayTrees = getDecayTrees(mom)
        if not decayTrees:
            mom.finalState = True
            continue

        # Add all decay channels to all the trees
        newTrees = []
        for T in treeList:
            for decay in decayTrees:
                # The order below matters,
                # since we want to keep the mother from the decay tree (which holds the BR value)
                newTree = nx.compose(decay, T)
                if sigmacut is not None:
                    treeWeight = getTreeWeight(newTree)
                    if treeWeight is not None and treeWeight < sigmacut:
                        # Do not consider this decay or the next ones,
                        # since they are sorted accodring to BR and the subsequent
                        # decays will only contain smaller weights
                        break
                newTrees.append(newTree)

        if not newTrees:
            continue
        treeList = newTrees

    if len(treeList) == 1 and treeList[0] == tree:
        return []
    else:
        return treeList


def cascadeDecay(tree, sigmacut=None):
    """
    Given a tree, generates a list of new trees (DiGraph objects),
    where all the particles have cascade decayed to stable final states.

    :param tree: Tree (DiGraph object) for which to add the decays
    :param sigmacut: Cut on the tree weight (xsec*BR). Any tree with weights
                     smaller than sigmacut will be ignored.

    :return: List of trees with all possible decays added.
    """

    treeList = [tree]
    newTrees = True
    finalTrees = []
    while newTrees:
        newTrees = []
        for T in treeList:
            newT = addOneStepDecays(T, sigmacut)
            if not newT:
                finalStates = getFinalStates(T)
                # Make sure all the final states have decayed
                # (newT can be empty if there is no allowed decay above sigmacut)
                if any(fs.finalState is False for fs in finalStates):
                    continue
                finalTrees.append(T)  # It was not possible to add any new decay to the tree
            else:
                newTrees += newT  # Add decayed trees to the next iteration

        if not newTrees:  # It was not possible to add any new decay
            break
        treeList = newTrees[:]

    return finalTrees


def drawTree(tree, particleColor='lightcoral',
             smColor='skyblue',
             pvColor='darkgray',
             nodeScale=4, labelAttr=None, attrUnit=None):
    """
    Draws Tree using matplotlib.

    :param tree: tree to be drawn
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
        labels = {n: str(n) for n in tree.nodes()}
    elif attrUnit is not None:
        labels = {n: str(getattr(n, labelAttr).asNumber(attrUnit)) if hasattr(n, labelAttr)
                  else str(n) for n in tree.nodes()}
    else:
        labels = {n: str(getattr(n, labelAttr)) if hasattr(n, labelAttr)
                  else str(n) for n in tree.nodes()}

    for key in labels:
        if labels[key] == 'anyOdd':
            labels[key] = 'BSM'
    node_size = []
    node_color = []
    for n in tree.nodes():
        node_size.append(nodeScale*100*len(labels[n]))
        if 'pv' == labels[n].lower():
            node_color.append(pvColor)
        elif hasattr(n, '_isSM') and n._isSM:
            node_color.append(smColor)
        else:
            node_color.append(particleColor)

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
