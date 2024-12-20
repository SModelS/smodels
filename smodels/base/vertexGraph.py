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

        :param incoming: List of Particle objects entering the vertex (i.e. X for X -> a,b)
        :param outgoing: List of Particle objects exiting the vertex (i.e. a,b for X -> a,b)
        """

        GenericGraph.__init__(self)
        self.all_particles_dict = {}
        self.incoming_nodes = []
        self.outgoing_nodes = []

        for p in incoming+outgoing:
            if (not hasattr(p,'pdg')) and p.getPdg() is None:
                raise SModelSError(f"VertexGraphs should only be created with particles with the pdg attribute (missing in {p})")
            if not hasattr(p,'isSM'):
                raise SModelSError(f"VertexGraphs should only be created with particles with the isSM attribute (missing in {p})")
        
        # Sort incoming and outgoing by pdg and iSM:
        incoming = sorted(incoming, 
                              key = lambda p: (not p.isSM,p.getPdg()), reverse=True)
        outgoing = sorted(outgoing, 
                              key = lambda p: (not p.isSM,p.getPdg()), reverse=True)
        
        # Store all particles as final states (outgoing) for fast comparison
        n = 0
        for p in incoming:
            nodeIndex = self.add_node(n)    
            self.all_particles_dict[nodeIndex]= p.chargeConjugate()
            self.incoming_nodes.append(nodeIndex)
            n += 1
        for p in outgoing:
            nodeIndex = self.add_node(n)
            self.all_particles_dict[nodeIndex]= p
            self.outgoing_nodes.append(nodeIndex)
            n += 1

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

        if sorted(self.inPDGs()) != sorted(other.inPDGs()):
            return False

        if sorted(self.outPDGs()) != sorted(other.outPDGs()):
            return False
        
        return True
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def indexToNode(self, nodeIndex):
        """
        Returns the particle with index nodeIndex.
        If nodeIndex is in outgoing_nodes, return particle
        stored in self.all_particles_dict, otherwise
        returns the conjugate of the particle in 
        self.all_particles_dict.

        :param nodeIndex: Integer

        :return: Particle object
        """
        
        if isinstance(nodeIndex,int):
            particle = self.all_particles_dict[nodeIndex]
            if nodeIndex in self.outgoing_nodes:
                return particle
            else:
                return particle.chargeConjugate()
        elif isinstance(nodeIndex,(list,tuple)):
            particles = []
            for n in nodeIndex:
                particle = self.all_particles_dict[n]
                if nodeIndex in self.outgoing_nodes:
                    particles.append(particle)
                else:
                    particles.append(particle.chargeConjugate())
            
            return particles
        else:
            raise SModelSError("Can not convert object of type %s to nodes" %str(type(nodeIndex)))
            
        

    def inPDGs(self):

        pdgList = []
        for nodeIndex in self.incoming_nodes:
            pdgList.append(self.indexToNode(nodeIndex).getPdg())
        
        return pdgList
    
    def outPDGs(self):

        pdgList = []
        for nodeIndex in self.outgoing_nodes:
            pdgList.append(self.indexToNode(nodeIndex).getPdg())
        
        return pdgList
    
    
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

        dictA = self.all_particles_dict
        dictB = other.all_particles_dict
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
        for nodeB in other.nodeIndices:
            matchedA = matchDict[nodeB]
            partA = self.indexToNode(matchedA)
            matchedVertex.add_node(partA,nodeB)
        matchedVertex.incoming_nodes = other.incoming_nodes[:]
        matchedVertex.outgoing_nodes = other.outgoing_nodes[:]

        self.add_edges_from(itertools.product(self.incoming_nodes,self.outgoing_nodes))
        
        return matchedVertex
