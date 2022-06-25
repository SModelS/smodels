"""
.. module:: treeMatcher
   :synopsis: A collection of classes and functions used to compare
              trees, compute their isomorphism and semantic equivalence.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import collections


class TreeMatcher(object):
    """
    Isomorphism checker for labelled trees. It starts at the root and checks for subtree
    isomosphisms. The tree and subtree topology is encoded in its encoding (canonical name),
    so the comparison is made based on its canonical name. If a set of neighboor subtrees
    have identical names, all the permutations are tested for isomorphism.
    The comparison relies on the method for node comparison, which compares canonical names
    and particle content.
    """

    def __init__(self, T1, T2):
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

    def swapTrees(self):
        """
        Swap the trees T1<->T2, so the comparison will be made as T2 as the base tree.
        """

        T1, T2 = self.T1, self.T2
        self.__init__(T2, T1)

    def compareSubTrees(self, T1_node, T2_node):
        """
        Compare the subtrees with T1_node and T2_node as roots.
        The comparison is made according to their names, particle content and final states.
        It uses the node comparison to define semantically equivalent nodes.

        :param T1_node: Particle node belonging to self.T1
        :param T2_node: Particle node belonging to self.T2

        :return: (cmp,matcDict), where: cmp = 0 (nodes are equal), 1 (T1_node > T2_node) or -1 (T1_node < T2_node)
                 and matchDict = None (nodes differ) or a dictionary with the mapping of the nodes daughters
                 ({d1 : d2}).
        """

        # Compare nodes directly (canon name and particle content)
        cmp = T1_node.compareTo(T2_node)
        if cmp != 0:
            return (cmp, None)

        # For inclusive nodes always return True (once nodes are equal)
        if T1_node.isInclusive or T2_node.isInclusive:
            return (0, {T1_node: T2_node})

        # Check for equality of daughters
        successors1 = self.T1.successors(T1_node)
        successors2 = self.T2.successors(T2_node)
        if len(successors1) == len(successors2) == 0:
            return (0, {T1_node: T2_node})

        # Sort the daughters. Inclusive nodes will always come first,
        # and the remaining nodes are sorted by canonName and particle, so if the nodes differ,
        # the comparison is independent of the original daughters position
        sortedDaughters1 = sorted(successors1,
                                  key=lambda n: (not n.isInclusive, n.canonName, n.particle))

        sortedDaughters2 = sorted(successors2,
                                  key=lambda n: (not n.isInclusive, n.canonName, n.particle))

        # Compare daughters directly:
        cmpDict = {}
        cmp_max = 0
        for d1, d2 in zip(sortedDaughters1, sortedDaughters2):
            cmp, mapDict = self.compareSubTrees(d1, d2)
            cmp_max = max(abs(cmp), cmp_max)  # Stores 1, if any pair of nodes differs
            cmpDict[d1] = {d2: (cmp, mapDict)}
        # Check for a direct match first
        # (likely since daughters are sorted):
        if cmp_max == 0 and len(sortedDaughters1) == len(sortedDaughters2):
            mapDict = {T1_node: T2_node}
            mapDict.update({d1: d2 for d1, d2 in zip(sortedDaughters1, sortedDaughters2)})
            for node1 in cmpDict:
                daughtersMap = list(cmpDict[node1].values())[0][1]
                mapDict.update(daughtersMap)
            return (0, mapDict)

        # Define left and right nodes in order to compute matching:
        left_nodes = sortedDaughters1[:]
        right_nodes = sortedDaughters2[:]
        edges = {}
        for d1 in left_nodes:
            if d1 not in cmpDict:
                cmpDict[d1] = {}
            for d2 in sortedDaughters2:
                if d2 in cmpDict[d1]:
                    cmp, mapDict = cmpDict[d1][d2]
                else:
                    cmp, mapDict = self.compareSubTrees(d1, d2)
                    cmpDict[d1].update({d2: (cmp, mapDict)})
                if cmp == 0:
                    if d1 not in edges:
                        edges[d1] = {d2: mapDict}
                    else:
                        edges[d1].update({d2: mapDict})

            # If node had no matches (was not added to the graph),
            # we already know T1_node and T2_node differs
            if d1 not in edges:
                break

        # Check if all nodes were included in the graph
        # (they had at least one match)
        if len(edges) != len(left_nodes):
            matched = False
        else:
            # Compute the maximal matching
            # (mapping where each node1 is connected to a single node2)
            mapDict = self.maximal_matching(left_nodes, right_nodes, edges)

            # Check if the match was successful.
            # Consider a successful match if all nodes in daughters1 were matched
            # or if all nodes in daughers2 were matched and the unmatched nodes
            # in daughters1 match an InclusiveNode.
            matched = False
            if len(mapDict) == len(left_nodes):
                matched = True
            elif len(mapDict) == len(right_nodes):
                matched = True
                unMatched = [d1 for d1 in left_nodes if d1 not in mapDict]
                for d1 in unMatched:
                    if any(not d2.isInclusive for d2 in edges[d1]):
                        matched = False
                        break

        if matched:
            # Include matches according to the selected edges:
            for d1, d2 in list(mapDict.items()):
                daughtersMap = edges[d1][d2]
                mapDict.update(daughtersMap)
            mapDict[T1_node] = T2_node

            # If all nodes appear in mapDict, one permutation was found
            # where all nodes match:
            return (0, mapDict)

        else:
            # Otherwise compare subtrees according to their sorted particles:
            for d1, d2 in zip(sortedDaughters1, sortedDaughters2):
                cmp, _ = cmpDict[d1][d2]
                if cmp != 0:
                    return (cmp, None)

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

        isInclusiveT1 = any(n.isInclusive for n in self.T1.nodes)
        isInclusiveT2 = any(n.isInclusive for n in self.T2.nodes)
        if isInclusiveT1:
            self.swapTrees()
            invert = not invert

        cmp, mappingDict = self.compareSubTrees(self.T1.root, self.T2.root)
        if invert:
            cmp = -cmp
        if cmp != 0:
            # If both trees are inclusive, try the opposite comparison:
            if isInclusiveT1 and isInclusiveT2:
                self.swapTrees()
                invert = not invert
                cmp_swap, _ = self.compareSubTrees(self.T1.root, self.T2.root)
                if cmp_swap != 0:
                    return cmp, None
            else:
                return cmp, None

        # Remove unmatched nodes (which happens in case of InclusiveNodes)
        match = {n1: n2 for n1, n2 in mappingDict.items() if n2 is not None}
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

    def maximal_matching(self, left, right, edges):
        """
        Computes the maximal matching from left nodes to right nodes.
        The maximal matching is the maximal number of left nodes which can be
        connected to the right nodes without any node belonging to more than one edge.
        Adpated from networkx.algorithms.bipartite.matching.hopcroft_karp_matching.

        :param left: List of left nodes
        :param right: List of right nodes
        :param edges: Nested dictionary with left nodes as keys and macthing right nodes as values
                      (e.g. {nL1 : {nR2 : {}, nR3 : {}}, nL2 : {nR2 : {}, nR1 : {}},... })
        """

        INFINITY = float("inf")
        # Initialize the "global" variables that maintain state during the search.
        leftmatches = {v: None for v in left}
        rightmatches = {v: None for v in right}
        distances = {}
        queue = collections.deque()

        def breadth_first_search():
            for v in left:
                if leftmatches[v] is None:
                    distances[v] = 0
                    queue.append(v)
                else:
                    distances[v] = INFINITY
            distances[None] = INFINITY
            while queue:
                v = queue.popleft()
                if distances[v] < distances[None]:
                    for u in edges[v]:
                        if distances[rightmatches[u]] is INFINITY:
                            distances[rightmatches[u]] = distances[v] + 1
                            queue.append(rightmatches[u])
            return distances[None] is not INFINITY

        def depth_first_search(v):
            if v is not None:
                for u in edges[v]:
                    if distances[rightmatches[u]] == distances[v] + 1:
                        if depth_first_search(rightmatches[u]):
                            rightmatches[u] = v
                            leftmatches[v] = u
                            return True
                distances[v] = INFINITY
                return False
            return True

        # Implementation note: this counter is incremented as pairs are matched but
        # it is currently not used elsewhere in the computation.
        num_matched_pairs = 0
        while breadth_first_search():
            for v in left:
                if leftmatches[v] is None:
                    if depth_first_search(v):
                        num_matched_pairs += 1

        # Strip the entries matched to `None`.
        leftmatches = {k: v for k, v in leftmatches.items() if v is not None}

        return leftmatches


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
        if tree.root.isInclusive:
            continue

        lo = 0
        hi = len(sortedTrees)
        while lo < hi:
            mid = (lo+hi)//2
            if tree.compareTreeTo(sortedTrees[mid])[0] > 0:
                lo = mid+1
            else:
                hi = mid
        sortedTrees.insert(lo, tree)

    sortedTrees = inclusiveTrees + sortedTrees
    return sortedTrees
