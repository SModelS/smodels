"""
.. module:: vertexGraph
   :synopsis: This is a base class for describing vertices as simple graphs.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


from smodels.base.exceptions import SModelSBaseError as SModelSError
from smodels.base.genericGraph import GenericGraph
import itertools


class VertexGraph(GenericGraph):

    """
    A class for describing and manipulating vertices described as simple graphs.
    """

    def __init__(self,incoming=[],outgoing=[]):
        """
        Initialize the vertex with the given incoming and outgoing particles.
        The canonical representation of the graph is a simple list of particles
        with all particles considered as outgoing, so incoming particles are conjugated.
        In addition if the particle with the largest PDG (in absolute value) is an anti-particle
        (pdg < 0), all the particles in the list are conjugated. 
        With these convention any cycle of the vertex as well as the conjugated vertex
        are represented by the same canonical representation.

        :param incoming: List of Particle objects entering the vertex (i.e. X for X -> a,b)
        :param outgoing: List of Particle objects exiting the vertex (i.e. a,b for X -> a,b)
        """

        GenericGraph.__init__(self)
        self._canonical_rep_dict = {}
        self._canonical_is_conjugated = False
        self.incoming_nodes = []
        self.outgoing_nodes = []

        for p in incoming+outgoing:
            if (not hasattr(p,'pdg')) and (p.getPdg() is None):
                raise SModelSError(f"VertexGraphs should only be created with particles with the pdg attribute (missing in {p})")
            if not hasattr(p,'isSM'):
                raise SModelSError(f"VertexGraphs should only be created with particles with the isSM attribute (missing in {p})")
            if not hasattr(p,'isSelfConjugate'):
                raise SModelSError(f"VertexGraphs should only be created with particles with the isSelfConjugated attribute (missing in {p})")
        
        # Sort incoming and outgoing by pdg and iSM:
        incoming = sorted(incoming, 
                              key = lambda p: (not p.isSM,p.getPdg()), reverse=True)
        outgoing = sorted(outgoing, 
                              key = lambda p: (not p.isSM,p.getPdg()), reverse=True)
        
        # Store all particles in the canonical rep as outgoing 
        for p in incoming:
            nodeIndex = self.add_node(p)    
            self._canonical_rep_dict[nodeIndex]= p.chargeConjugate()
            self.incoming_nodes.append(nodeIndex)
        for p in outgoing:
            nodeIndex = self.add_node(p)
            self._canonical_rep_dict[nodeIndex]= p
            self.outgoing_nodes.append(nodeIndex)

        # Check if the maximum PDG out of the non self-conjufate
        # particles is negative
        allPDGs = [p.getPdg() 
                    for p in self._canonical_rep_dict.values() if not p.isSelfConjugate]
        # If not, conjugate the canonical representation
        if allPDGs and max(allPDGs) != max([abs(pdg) for pdg in allPDGs]):
            self._canonical_is_conjugated = True
            for nodeIndex,p in self._canonical_rep_dict.items():
                self._canonical_rep_dict[nodeIndex] = p.chargeConjugate()

        self.add_edges_from(itertools.product(self.incoming_nodes,self.outgoing_nodes))


    def __repr__(self):
        """
        Returns the string representation of the graph.
        """

        return str(self)
    
    def __str__(self):
        """
        Returns a string listing the graph nodes.
        """

        inStr = ','.join([str(self.indexToNode(inNode)) for inNode in self.incoming_nodes])
        outStr = ','.join([str(self.indexToNode(inNode)) for inNode in self.outgoing_nodes])

        gStr = inStr + " > " + outStr

        return gStr
    
    def __eq__(self, other):
        """
        Define vertex equality based on the vertex particles' PDGs.

        :return: True/False
        """

        if len(self._canonical_rep_dict) != len(other._canonical_rep_dict):
            return False
        
        pdgsA = [p.getPdg() for p in self._canonical_rep_dict.values()]
        pdgsB = [p.getPdg() for p in other._canonical_rep_dict.values()]
        
        return sorted(pdgsA) == sorted(pdgsB)
        
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def chargeConjugate(self):
        """
        Returns a new vertex with the charge conjugation of
        all the incoming and outgoing particles
        """

        incoming = [self.indexToNode(nodeIndex).chargeConjugate()
                    for nodeIndex in self.incoming_nodes]
        
        outgoing = [self.indexToNode(nodeIndex).chargeConjugate()
                    for nodeIndex in self.outgoing_nodes]


        conjVertex = VertexGraph(incoming=incoming, outgoing=outgoing)

        return conjVertex

    
    def matchTo(self,other):
        """
        Check if it matches other, i.e. if there is any cycle permutation of
        the particles in self which matches other.
        If a match is found returns a new vertex with the incoming and outgoing
        particles sorted according to other.
        The comparison is made based on a flat list of particles from self and other
        where all particles are outgoing.
        
        :param other: VertexGraph object.
        
        :return: a new VertexGraph object with the incoming and outgoing particles
                 ordered according to other.
        """

        dictA = self._canonical_rep_dict
        dictB = other._canonical_rep_dict
        if len(dictA) != len(dictB):
            return None

        # Get all permutations of particlesA:
        permutations = itertools.permutations(dictA.items())
        matchDict = {}
        # Go over permutations until there is a full match
        while len(matchDict) != len(dictB):
            matchDict = {}
            # Try to get the next permutation
            # If done, exit the while loop
            try:
                permA = next(permutations)
            except StopIteration:
                break
            
            # Check if particlesA = particlesB
            for i,nodeB in enumerate(dictB):
                partB = dictB[nodeB]
                nodeA,partA = permA[i]
                comp = partA.cmpProperties(partB,
                                            properties=['isSM','pdg'])
                
                # If particles do not match, go to next permutation
                if comp != 0:
                    break
                matchDict[nodeB] = nodeA

        # If after going over all permutations
        # not all particles were matched, return None    
        if len(matchDict) != len(dictB):
            return None
        
        # Create new vertex following the indices in other:
        matchedVertex = VertexGraph()
        matchedVertex._canonical_is_conjugated = other._canonical_is_conjugated
        for nodeB in other.nodeIndices:
            nodeA = matchDict[nodeB]
            partA = dictA[nodeA]
            matchedVertex._canonical_rep_dict[nodeB] = partA
            # If it is an incoming node, add the conjugated particle
            if nodeB in other.incoming_nodes:
                partA = partA.chargeConjugate()
            # If the canonical represation of other includes charge conjugation
            # conjugate the particle
            if other._canonical_is_conjugated:
                partA = partA.chargeConjugate()
            matchedVertex.add_node(partA,nodeB)
        matchedVertex.incoming_nodes = other.incoming_nodes[:]
        matchedVertex.outgoing_nodes = other.outgoing_nodes[:]

        self.add_edges_from(itertools.product(self.incoming_nodes,self.outgoing_nodes))
        
        return matchedVertex
