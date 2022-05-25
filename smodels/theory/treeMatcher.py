"""
.. module:: treeMatcher
   :synopsis: A collection of classes and functions used to compare
              trees, compute their isomorphism and semantic equivalence.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import networkx.algorithms.isomorphism.isomorphvf2 as vf2
from smodels.tools.smodelsLogging import logger

class TreeMatcher(vf2.DiGraphMatcher):
    """
    VF2 isomorphism checker for directed graphs.
    It uses the node comparison method to define equivalent nodes.
    """

    def __init__(self, G1, G2):
        """
        Initialize graph matcher.

        Parameters
        ----------
        G1, G2 : graph
            The graphs to be tested.

        """
        vf2.DiGraphMatcher.__init__(self, G1, G2)

        # These will be modified during checks to minimize code repeat.
        self.G1_adj = self.G1.adj
        self.G2_adj = self.G2.adj

    def semantic_feasibility(self, G1_node, G2_node):
        """
        Returns True if G1_node to G2_node are semantically feasible.
        It uses the node comparison to define semantically equivalent nodes.
        """

        # Test node_match and also test edge_match on successors
        feasible = G1_node.equalTo(G2_node)
        if not feasible:
            return False

        # Test edge_match on predecessors
        self.G1_adj = self.G1.pred
        self.G2_adj = self.G2.pred
        feasible = G1_node.equalTo(G2_node)
        self.G1_adj = self.G1.adj
        self.G2_adj = self.G2.adj

        return feasible

    def getMatches(self, invert=False):
        """
        Compare self.G1 to self.G2, which contains an inclusive node.
        There is a match if G2 is a subtree of G1
        and the inclusive node matches the equivalent node of the subtree.

        The matched tree has its tree reduced to the subtree,
        with an InclusiveNode.

        :param other: Tree containing an InclusiveParticleNode
        :param invert: If True, will return the copy of the
                       G1 subtree correspoding to G2.

        :return: Iterator over the matching trees. If no matches were found, return None.
        """

        treeA = self.G1
        treeB = self.G2
        inclusiveA = any(n.isInclusive for n in treeA.nodes)
        inclusiveB = any(n.isInclusive for n in treeB.nodes)
        if inclusiveA and inclusiveB:
            msg = 'Comparing two trees with inclusive nodes.'
            msg += 'The comparison will likely produce wrong results!'
            logger.warning(msg)
        if inclusiveA and not inclusiveB:
            matcher = TreeMatcher(treeB, treeA)
            treeA = self.G2
            treeB = self.G1
            invert = not invert
        else:
            matcher = self

        inclusive = (inclusiveA or inclusiveB)

        # Get the isomorphism matches (if inclusive, other = subtree of self)
        if inclusive:
            matches = matcher.subgraph_isomorphisms_iter()
        else:
            matches = matcher.match()

        for match in matches:
            if invert:
                match = {v: k for k, v in match.items()}
                sourceTree = treeA
                baseTree = treeB
            else:
                sourceTree = treeB
                baseTree = treeA

            # Make a copy of the source tree
            matchedTree = sourceTree.copyTree()

            # Remove all nodes and edges
            matchedTree.clear()
            # Store the nodes in the correct order:
            newNodes = []
            for n in baseTree.nodes:  # match = {baseTree : sourceTree}
                if n not in match:  # In case of inclusive nodes the match is partial
                    continue
                if n.isInclusive:
                    newNode = n.copy()
                else:
                    # Make a copy of the original node
                    newNode = match[n].copy()
                newNode.node = n.node  # Force the node numbering to be equal
                newNodes.append(newNode)
                match[n] = newNode  # Update the dictionary with the copied node

            # Store the edges in the correct order:
            newEdges = [(match[nA], match[nB]) for nA, nB in baseTree.edges
                        if (nA in match and nB in match)]
            # Add the nodes and the edges to the new tree:
            matchedTree.add_nodes_from(newNodes)
            matchedTree.add_edges_from(newEdges)
            matchedTree.setGlobalProperties()
            yield matchedTree

        # If there was not good match, return None
        return None


def sortTreeList(treeList):
    """
    Sort a list of Tree objects according to the compareTrees method.

    :param treeList: List of Tree objects

    :return: Sorted list of Tree objects
    """

    # Create a listed of sorted nodes with proper canonical names according to the node comparison:
    # Always put the inclusive trees at the beginning
    inclusiveTrees = [t for t in treeList if t.root.isInclusive]
    sortedTrees = []
    for tree in treeList:
        if tree.root.isInclusive:
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
    if rootA.isInclusive:
        return 0
    elif rootB.isInclusive:
        return 0

    cmp = rootA.compareNode(rootB)
    # Return comparison if roots differ or root has no daughters
    if cmp != 0:
        return cmp

    # If root has no daughters and roots agree, then trees are equal
    if treeA.out_degree(rootA) == 0:
        return 0

    # Get daughters and sort them, so the comparison is unique:
    subTreesA = [treeA.getSubTree(source=daughter) for daughter in treeA.successors(rootA)]
    subTreesA = sortTreeList(subTreesA)
    subTreesB = [treeB.getSubTree(source=daughter)  for daughter in treeB.successors(rootB)]
    subTreesB = sortTreeList(subTreesB)

    # Loop over sorted daughters and stop when the first differ:
    for itree, subtreeA in enumerate(subTreesA):
        subtreeB = subTreesB[itree]
        cmp = compareTrees(subtreeA, subtreeB)
        if cmp != 0:
            return cmp

    # If all subtrees were identical, return zero
    return 0
