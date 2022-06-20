"""
.. module:: treeMatcher
   :synopsis: A collection of classes and functions used to compare
              trees, compute their isomorphism and semantic equivalence.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import itertools


class TreeMatcher(object):
    """
    Isomorphism checker for labelled trees. It starts at the root and checks for subtree
    isomosphisms. The tree and subtree topology is encoded in its encoding (canonical name),
    so the comparison is made based on its canonical name. If a set of neighboor subtrees
    have identical names, all the permutations are tested for isomorphism.
    The comparison relies on the method for node comparison, which compares canonical names
    and particle content.
    """

    def __init__(self, T1, T2, inclusive=True):
        """
        Initialize graph matcher.

        :parameter T1: Tree object
        :parameter T2: Tree object
        :param inclusive: If False, particles are required to be identical
                          (the inclusiveness of MultiParticles or InclusiveNodes
                          are not considered when comparing)


        :return: A dictionary mappinn the nodes of treeA and treeB ({nA : nB}).

        """

        # These will be modified during checks to minimize code repeat.
        self.T1 = T1
        self.T2 = T2
        self._comps = {n.node: {} for n in T1.nodes}  # Cache node comparison (for debugging)
        self.mappingDict = {n: None for n in T1.nodes}
        self.inclusive = inclusive  # Whether to allow for inclusive labels or not

    def swapTrees(self):
        """
        Swap the trees T1<->T2, so the comparison will be made as T2 as the base tree.
        """

        T1, T2 = self.T1, self.T2
        self.__init__(T2, T1, self.inclusive)

    def compareSubTrees(self, T1_node, T2_node):
        """
        Compare the subtrees with T1_node and T2_node as roots.
        The comparison is made according to their names, particle content and final states.
        It uses the node comparison to define semantically equivalent nodes.

        :param T1_node: Particle node belonging to self.T1
        :param T2_node: Particle node belonging to self.T2

        :return: 0 if nodes are equal, 1 if T1_node > T2_node, -1 otherwise
        """

        if T2_node.node in self._comps[T1_node.node]:
            return self._comps[T1_node.node][T2_node.node]

        # Compare nodes directly (canon name and particle content)
        cmp = T1_node.compareTo(T2_node, self.inclusive)
        if cmp != 0:
            self._comps[T1_node.node].update({T2_node.node: cmp})
            return cmp

        # For inclusive nodes always return True (once nodes are equal)
        if T1_node.isInclusive or T2_node.isInclusive:
            self._comps[T1_node.node].update({T2_node.node: 0})  # Cache comparison
            self.mappingDict.update({T1_node: T2_node})
            return 0

        # Check for equality of daughters
        successors1 = list(self.T1.successors(T1_node))
        successors2 = list(self.T2.successors(T2_node))

        # Sort the daughters. Inclusive nodes will always come first,
        # and the remaining nodes are sorted by canonName and particle, so if the nodes differ,
        # the comparison is independent of the original daughters position
        sortedDaughters1 = sorted(successors1,
                                  key=lambda n: (not n.isInclusive, n.canonName, n.particle))

        sortedDaughters2 = sorted(successors2,
                                  key=lambda n: (not n.isInclusive, n.canonName, n.particle))

        # Get comparison dictionary for daughters:
        # compDict = {}
        # matchesDict = {}
        # for d1 in sortedDaughters1:
        # compDict[d1.node] = {d2.node: self.compareSubTrees(d1, d2)
        # for d2 in sortedDaughters2}
        # compValues = list(compDict.values())
        # If d2 does not match ANY node in daughters1, stop comparison
        # if 0 not in compValues:
        # cmp = max(compValues)
        # self._comps[T1_node.node].update({T2_node.node: cmp})  # Cache comparison
        # return cmp
        # else:
        # matchesDict[d2.node] = compDict

        # Check all permutations within each set of daughters with the
        # same canonical name
        allPerms = []
        for key, group in itertools.groupby(sortedDaughters2, lambda d: d.canonName):
            allPerms.append(itertools.permutations(group))

        # Construct all permutations and check against daughters1:
        cmp_max = -10
        for daughters2_perm in itertools.product(*allPerms):
            daughters2_perm = itertools.chain.from_iterable(daughters2_perm)
            daughters2_perm = list(daughters2_perm)

            for i2, d2 in enumerate(daughters2_perm):
                d1 = sortedDaughters1[i2]
                # cmp = compDict[d1.node][d2.node]
                cmp = self.compareSubTrees(d1, d2)
                if cmp != 0:
                    cmp_max = max(cmp, cmp_max)
                    break
            else:
                # Found one permutation where all nodes match:
                self._comps[T1_node.node].update({T2_node.node: 0})  # Cache comparison
                self.mappingDict.update({T1_node: T2_node})
                return 0

        self._comps[T1_node.node].update({T2_node.node: cmp_max})  # Cache comparison
        # Return maximum comparison value,
        # so the result is independent of the nodes ordering
        return cmp_max

    def compareTrees(self, invert=False):
        """
        Compare self.T1 to self.T2. Returns an isomorphism of self.T2 which
        matches T1. If T2 contains an inclusive node, T2 can be a subtree of T1 (T2 C T1).
        If T1 contains an inclusive node, T1 can be a subtree of T2 (T1 C T2).
        If both T1 and T2 have inclusive nodes, both possibilites are tried out:
        T1 C T2 or T2 C T1.

        :param invert: If True, will return a copy of the
                       T1 tree matching T2.

        :return: Returns the comparison value and the matching tree.
                 If no matches were found (cmp != 0), return None and the comparison value.
        """

        isInclusiveT1 = any(n.isInclusive for n in self.T1)
        isInclusiveT2 = any(n.isInclusive for n in self.T2)
        if isInclusiveT1:
            self.swapTrees()
            invert = not invert
        # print('First check', invert)
        cmp = self.compareSubTrees(self.T1.root, self.T2.root)
        if invert:
            cmp = -cmp
        if cmp != 0:
            # If both trees are inclusive, try the opposite comparison:
            if isInclusiveT1 and isInclusiveT2:
                self.swapTrees()
                invert = not invert
                # print('Second check', invert)
                cmp_swap = self.compareSubTrees(self.T1.root, self.T2.root)
                if cmp_swap != 0:
                    return cmp, None
            else:
                return cmp, None

        # Remove unmatched nodes (which happens in case of InclusiveNodes)
        match = {n1: n2 for n1, n2 in self.mappingDict.items() if n2 is not None}
        # Define sourceTree (tree to be returned with its ordering
        # and node numbering matching the baseTree)
        if not invert:
            sourceTree = self.T2
            baseTree = self.T1
        else:
            match = {v: k for k, v in match.items()}
            sourceTree = self.T1
            baseTree = self.T2

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
            newNode.node = n.node   # Force the node numbering to be equal
            newNodes.append(newNode)
            match[n] = newNode  # Update the dictionary with the copied node

        # Store the edges in the correct order:
        newEdges = [(match[nA], match[nB]) for nA, nB in baseTree.edges
                    if (nA in match and nB in match)]
        # Add the nodes and the edges to the new tree:
        matchedTree.add_nodes_from(newNodes)
        matchedTree.add_edges_from(newEdges)
        matchedTree.setGlobalProperties()

        return 0, matchedTree


def sortTreeList(treeList):
    """
    Sort a list of Tree objects according to the compareTrees method.

    :param treeList: List of Tree objects

    :return: Sorted list of Tree objects
    """

    # Create a listed of sorted nodes with proper canonical names according to the node comparison.
    # Make sure trees which have an include node as root are always at the beginning
    inclusiveTrees = [t for t in treeList if t.root.isInclusive]
    sortedTrees = []
    for tree in treeList:
        itree = 0
        while itree < len(sortedTrees) and tree.compareTreeTo(sortedTrees[itree])[0] > 0:
            itree += 1
        sortedTrees.insert(itree, tree)

    sortedTrees = inclusiveTrees + sortedTrees
    return sortedTrees
