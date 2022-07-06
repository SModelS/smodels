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
    isomosphisms. The tree and subtree topology is encoded in its canonical name,
    so the comparison is made based on it. If a set of neighboor subtrees
    have identical names, all the permutations are tested for isomorphism.
    The comparison relies on the method for node comparison, which compares canonical names
    and particle content.
    """

    def __init__(self, T1, T2):
        """
        Initialize graph matcher.

        :parameter T1: Tree object
        :parameter T2: Tree object
        """

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

        :param T1_node: Node index belonging to self.T1
        :param T2_node: Node index belonging to self.T2

        :return: matcDict is None (subtrees differ) or a dictionary
                 with the mapping of the nodes daughters ({d1 : d2}).
        """

        # Compare nodes directly (canon name and particle content)
        node1 = self.T1.nodesMapping[T1_node]
        node2 = self.T2.nodesMapping[T2_node]
        if node1.isInclusive:
            self.T2.getFinalStates(node2.node)  # Make sure final states are defined
        if node2.isInclusive:
            self.T1.getFinalStates(node1.node)  # Make sure final states are defined

        cmp = node1.compareTo(node2)
        if cmp != 0:
            return None

        # For inclusive nodes always return True (once nodes are equal)
        if node1.isInclusive or node2.isInclusive:
            return {T1_node: T2_node}

        # Check for equality of daughters
        daughters1 = self.T1.successors[T1_node]
        daughters2 = self.T2.successors[T2_node]
        if len(daughters1) == len(daughters2) == 0:
            return {T1_node: T2_node}

        # Compare daughters directly:
        matchesDict = {}
        directMatch = True
        for d1, d2 in zip(daughters1, daughters2):
            mapDict = self.compareSubTrees(d1, d2)
            matchesDict[d1] = {d2: mapDict}  # Store the comparison
            if mapDict is None:
                directMatch = False

        # Check for a direct match first
        if directMatch and len(daughters1) == len(daughters2):
            mapDict = {T1_node: T2_node}
            mapDict.update({d1: d2 for d1, d2 in zip(daughters1, daughters2)})
            for node1 in matchesDict:
                daughtersMap = list(matchesDict[node1].values())[0]
                mapDict.update(daughtersMap)
            return mapDict

        # Define left and right nodes in order to compute matching:
        left_nodes = daughters1[:]
        right_nodes = daughters2[:]
        edges = {}
        for d1 in left_nodes:
            if d1 not in matchesDict:
                matchesDict[d1] = {}
            for d2 in right_nodes:
                if d2 in matchesDict[d1]:
                    mapDict = matchesDict[d1][d2]
                else:
                    mapDict = self.compareSubTrees(d1, d2)
                    matchesDict[d1].update({d2: mapDict})

                # If daughters do not match, do not include in edges
                if mapDict is None:
                    continue
                if d1 not in edges:
                    edges[d1] = {d2: mapDict}
                else:
                    edges[d1].update({d2: mapDict})

            # If node had no matches (was not added to the graph),
            # we already know T1_node and T2_node differs
            if d1 not in edges:
                return None

        # Check if all nodes were included in the graph
        # (they had at least one match)
        if len(edges) != len(left_nodes):
            return None
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
                    if any(not self.T2.nodesMapping[d2].isInclusive
                           for d2 in edges[d1]):
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
            return mapDict

        else:
            return None

    def matchTrees(self, invert=False):
        """
        Compare self.T1 to self.T2. Returns an isomorphism of self.T2 which
        matches T1. If T2 contains an inclusive node, T2 can be a subtree of T1 (T2 C T1).
        If T1 contains an inclusive node, T1 can be a subtree of T2 (T1 C T2).
        If both T1 and T2 have inclusive nodes, both possibilites are tried out:
        T1 C T2 or T2 C T1.

        :param invert: If True, will return a copy of the
                       T1 tree matching T2.

        :return: Returns the matched tree, if trees match. Otherwise returns None.
        """

        isInclusiveT1 = any(n.isInclusive for n in self.T1.nodes)
        isInclusiveT2 = any(n.isInclusive for n in self.T2.nodes)
        if isInclusiveT1:
            self.swapTrees()
            invert = not invert

        mappingDict = self.compareSubTrees(self.T1.root.node, self.T2.root.node)
        if mappingDict is None:
            # If both trees are inclusive, try the opposite comparison:
            if isInclusiveT1 and isInclusiveT2:
                self.swapTrees()
                invert = not invert
                mappingDict = self.compareSubTrees(self.T1.root.node, self.T2.root.node)
                if mappingDict is None:
                    return None
            else:
                return None

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

        # Make a new tree
        matchedTree = sourceTree.copyTree(emptyNodes=True)

        # Store the nodes in the correct order:
        newNodes = []
        for n in baseTree.nodes:  # match = {baseTree : sourceTree}
            if n.node not in match:  # In case of inclusive nodes the match is partial
                continue
            elif n.isInclusive:
                newNode = n.copy()
            else:
                # Make a copy of the original node
                newNode = sourceTree.nodesMapping[match[n.node]].copy()
            newNode.node = n.node   # Force the node numbering to be equal
            newNodes.append(newNode)
            match[n.node] = newNode  # Update the dictionary with the copied node

        # Store the edges in the correct order:
        newEdges = [(match[nA.node], match[nB.node]) for nA, nB in baseTree.edges
                    if (nA.node in match and nB.node in match)]
        # Add the nodes and the edges to the new tree:
        matchedTree.add_nodes_from(newNodes)
        matchedTree.add_edges_from(newEdges)

        return matchedTree

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
