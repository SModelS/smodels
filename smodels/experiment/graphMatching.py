"""
.. module:: graphMatching
   :synopsis: A collection of functions to compute the matching between graphs

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from smodels.base.smodelsLogging import logger

from collections import deque


def getCycle(G):
    """
    Given a directed graph G, return a cycle, if it
    exits.
    
    :param G: Dictionary with the directed edges 
              ({A : [B,C,D], B : [],..})
    
    :return: List of nodes generating the cycle.
             The first and last entries are the
             same node.
             ([A,C,D,A])
    """

    
    visited=set()    
    stack = []
    
    def dfs_visit(u):

        # Add u to the current stack of visited nodes
        stack.append(u)                         
        for v in G[u]:
            if v in stack: # A cycle has been found with a path starting and ending in v                
                stack.append(v) # Include v at the end of the stack
                return True
            if v in visited:
                continue
            cycleFound = dfs_visit(v)
            if cycleFound:
                return True
        visited.add(u)
        stack.remove(u)
        return False    
                                 
    for u in G:
        # If u has been visited, it means no cycle 
        # can have u as the initial node, so we can skip it.
        if u in visited: 
            continue
        # Visit all the descendents of u
        # (if any the "family trees" has a node
        # with a descendent which returns to the node, a cycle
        # has been found)
        cycleFound = dfs_visit(u)
        if cycleFound:
            break
    
    if not cycleFound:
        return None

    # If a cycle has been found the cycle is a subset
    # of the current stack (stack = A -> B -> C -> D -> E -> C)
    # where the last entry is the beginning and end of the cycle
    headNode = stack[-1]
    path = stack[stack.index(headNode):]
    
    return path

def getDirectedEdges(edges,match):    
    """
    From the edges of a bipartite graph (left -> right)
    and a match dictionary, split the edges into left -> right
    (if they appear in the match) and right -> left (otherwise).

    :param edges: Dictionary with all edges between left and right nodes
                  (left nodes as keys and right nodes as values)
    :param match: Dictionary containing the edges for one perfect matching 
                  (left nodes as keys and right nodes as values)
    
    :return: Dictionary of directed edges between left and right nodes
             with edges from matchDict pointing from left to right
             and all the other edges pointing from right to left.
    """

    edgesL2R = {}
    edgesR2L = {}
    
    for nL in edges:
        for nR in edges[nL]:
            if match[nL] == nR:
                if nL not in edgesL2R:
                    edgesL2R[nL] = []
                edgesL2R[nL].append(nR)
            else:
                if nR not in edgesR2L:
                    edgesR2L[nR] = []
                edgesR2L[nR].append(nL)
                
    return edgesL2R,edgesR2L    

def maximal_matching(left, right, edges):
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
    queue = deque()

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

def getNewMatch(left,right,edges,match):
    """
    Given a perfect match, a list of left and right nodes in a bipartite
    graph and their edges (left > right), compute a new match.
    
    :param left: list of left nodes ([nL1,nL2,...])
    :param right: list of right nodes ([nR1,nR2,..])
    :param edges: Dictionary with left->right edges
                  ({nL1 : [nR2,nR3], nL2 : [nR1],...})
    :param match: Dictionary with a perfect matching
                  ({nL1 : nR3, nL2 : nR1,...})
    """
        
    # First define the directed edges
    # (edges in matchDict are left -> right and all the other
    # ones are right -> left)
    edgesL2R,edgesR2L = getDirectedEdges(edges,match)
    
    # For bookkepping create a dictionary with
    # all nodes and directed edges.
    # In order to distinguish between left and right
    # nodes, convert them to tuples:
    # (0,nL) -> left node, (1,nR) -> right node    
    G = {}
    for nL in left:
        if nL in edgesL2R:
            G[(0,nL)] = [(1,nR) for nR in edgesL2R[nL]]
        else:
            G[(0,nL)] = []
    for nR in right:
        if nR in edgesR2L:
            G[(1,nR)] = [(0,nL) for nL in edgesR2L[nR]]
        else:
            G[(1,nR)] = []
            
    # Find one cyle in the directed graph
    # (if a cycle exists, it means a new match can be created)
    path = getCycle(G)
    if not path:
        return None

    # Create a new match removing the left->right
    # edges in the cycle and replacing them by the right->left
    # edges
    new_match = dict(match.items())
    removeEdges = {}
    addEdges = {}
    for n1,n2 in list(zip(path[:-1],path[1:])):
        if n1[0] == 0: # cycle has left -> right (to be removed)            
            l = n1[1]
            r = n2[1]
            removeEdges[l] = r            
        else: # cycle has right -> left (to be added)
            r = n1[1]
            l = n2[1]
            addEdges[l] = r
    
    for l in removeEdges:
        new_match.pop(l)
        
    new_match.update(addEdges)
    if len(new_match) != len(match):
        raise SModelSError("Error constructing new match dict")
    
    return new_match

def getPerfectMatchings(left,right,edges):
    """
    Find all perfect matchings in an undirected bipartite graph.

    :param left: List of left nodes
    :param right: List of right nodes
    :param edges: Dictionary with all edges between left and right nodes
                  (left nodes as keys and right nodes as values)

    :return: list with all matching dictionaries (left nodes as keys and right nodes as values)
    """

    all_matches=[]

    # Find the first match
    match = maximal_matching(left,right,edges)
    if len(match) != len(left):
        return []

    match = dict(sorted(match.items()))
    all_matches.append(match)
    # Enter recursion
    all_matches = perfectMatchingIter(left,right,edges,match,all_matches,None)

    return iter(all_matches)

def perfectMatchingIter(left,right,edges,match,all_matches,add_e=None):
    """
    Iterate over all perfect matchings.

    :param left: List of left nodes
    :param right: List of right nodes
    :param edges: Dictionary with all edges between left and right nodes
                  (left nodes as keys and right nodes as values)
    :param match: Dictionary with a perfect matching 
                  (left nodes as keys and right nodes as values)
    :param all_matches: list with pertfect matching dicts.
                   Newly found matchings will be appendedinto this list.
    :param add_e: List of tuples with the edges used to form subproblems. If not None,
                  will be added to each newly found matchings.

    :return all_matches: updated list of all perfect matchings.
    """

    # Find a new match using the current match:
    new_match = getNewMatch(left,right,edges,match)

    if new_match is None:
        return all_matches
    else:
        
        # Add previously removed edges to the new match:
        if add_e is not None:
            for e in add_e:
                new_match[e[0]] = e[1]

        # Sort new_match for convenience
        new_match = dict(sorted(new_match.items()))
        all_matches.append(new_match)

        # Choose an edge from the original match not present in the new_match
        nL,nR = (set(match.items()).difference(set(new_match.items()))).pop()
        
        # Form subproblems
        # Create G+ with the nodes nL,nR removed and all of
        # its edges
        left_plus = [n for n in left if n != nL]
        right_plus = [n for n in right if n != nR]
        edges_plus = {n : nList[:] for n,nList in edges.items() if n != nL}
        for n,nList in edges_plus.items():
            if nR not in nList:
                continue
            nList.remove(nR)
            edges_plus[n] = nList[:]
        
        # Create G- with the edge nL->nR removed:
        left_minus = left[:]
        right_minus = right[:]
        edges_minus = {n : nList[:] for n,nList in edges.items()}
        edges_minus[nL].remove(nR)
        
        add_e_new = [(nL,nR)]
        if add_e is not None:
            add_e_new.append(add_e)

        all_matches = perfectMatchingIter(left_minus,right_minus,edges_minus,new_match,all_matches,add_e)
        all_matches = perfectMatchingIter(left_plus,right_plus,edges_plus,match,all_matches,add_e_new)

    return all_matches        
