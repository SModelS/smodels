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
        The Vertex is always created with a single incoming particle and all the other
        particles are set as outgoing (through charge conjugation, if necessary).
        In addition the vertex is rearranged to a canonical form: the BSM particle
        with largest pdg is the incoming particle and if the incoming particle pdg
        is negative, the vertex is conjugated (so the incoming pdg is always positive).

        :param incoming: List of Particle objects entering the vertex (i.e. X for X -> a,b)
        :param outgoing: List of Particle objects exiting the vertex (i.e. a,b for X -> a,b)
        """

        GenericGraph.__init__(self)

        # Create a graph with a single incoming particle
        # (the BSM particle with largest pdg)
        ## First define all particles as incoming:
        all_particles = list(incoming) + [p.chargeConjugate() for p in outgoing]
        for p in all_particles:
            if not hasattr(p,'pdg'):
                raise SModelSError(f"VertexGraphs should only be created with particles with the pdg attribute (missing in {p})")
            if not hasattr(p,'isSM'):
                raise SModelSError(f"VertexGraphs should only be created with particles with the isSM attribute (missing in {p})")
        
        
        ## Sort all particles by BSM and pdg values
        ## (for MultiParticles use the largest pdg for sorting)
        all_particles = sorted(all_particles, 
                              key = lambda p: (not p.isSM,abs(p.getPdg())), reverse=True)
        if len(all_particles) > 0:
            in_particle = all_particles[0]
            out_particles = [p.chargeConjugate() for p in all_particles[1:]]

            # Now if in_particle only has negative PDGs, conjugate the vertex:
            if in_particle.getPdg() < 0:
                in_particle = in_particle.chargeConjugate()
                out_particles = [p.chargeConjugate() for p in out_particles[:]]

            # Finally sort outgoing particles again
            out_particles = sorted(out_particles, 
                                    key = lambda p: (not p.isSM,abs(p.getPdg())), 
                                    reverse=True)

            self.add_node(in_particle,0)
            self.add_nodes_from(out_particles)
            # Add an edge from the incoming particle to all the outgoing ones
            self.add_edges_from([(0,i) for i in range(1,len(out_particles)+1)])


    def __repr__(self):
        """
        Returns the string representation of the graph.
        """

        return str(self)
    
    def __str__(self):
        """
        Returns a string listing the graph nodes.
        """

        gStr = str(self.indexToNode(0)) + " > "
        gStr += ",".join([str(d) for d in self.daughters(0)])

        return gStr
    
    def __eq__(self, other):
        """
        Define vertex equality based on the vertex PDGs.

        :return: True/False
        """

        momA_pdg = self.indexToNode(0).pdg
        momB_pdg = other.indexToNode(0).pdg

        if momA_pdg != momB_pdg:
            return False
        
        daughtersA_pdg = sorted([p.pdg for p in self.daughters(0)])
        daughtersB_pdg = sorted([p.pdg for p in other.daughters(0)])

        if daughtersA_pdg != daughtersB_pdg:
            return False
        
        return True
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def cycle(self,n=1):
        """
        Rotate the vertex by n cycles .
        For n=1: vertex A > B, C goes to vertex C > A, B).
        For n=2: vertex A > B, C goes to vertex B > C, A)...

        :param n: Number of cycles

        :return: Rotated GraphVertex object
        """

        old_indices = self.nodeIndices()
        newNode_indices = old_indices[n:]
        newNode_indices += old_indices[:n]

        newVertex = VertexGraph()
        for i,nodeIndex in enumerate(newNode_indices):
            d = self.indexToNode(nodeIndex)
            # Charge conjugate
            # if a final state is moving to an initial state
            if i == 0 and nodeIndex != 0:
                d = d.chargeConjugate()
            # Charge conjugate
            # if an initial state is moving to a final state
            elif nodeIndex == 0 and i != 0:
                d = d.chargeConjugate()
            newVertex.add_node(d,i)
            # Add edge if the node is not the initial state 
            if i != 0:
                newVertex.add_edge(0,i)

        return newVertex
    
    def matchTo(self,other):
        """
        Check if it matches other, i.e. if the incoming particles match and
        if all of the outgoing particles match (irrespective of the ordering).
        The matching is based only on the Particle's pdgs and whether they are
        a SM or BSM particle.
        
        :param other: VertexGraph object.
        
        :return: a new VertexGraph object with the daughters ordered according to
                 other.
        """


        # Get node objects
        mom1 = self.indexToNode(0)
        mom2 = other.indexToNode(0)

        # Compare incoming particles
        comp = mom1.cmpProperties(mom2,
                                properties=['isSM','pdg'])
        if comp != 0:
            return None
        
        daughters1 = self.daughters(0)
        daughter2_indices = other.daughterIndices(0)
        daughters2 = other.indexToNode(daughter2_indices)
        # GO over all permutations of daughters1:
        for daughters1_perm in itertools.permutations(daughters1):
            daughters1_perm = list(daughters1_perm)
            matches = {i2 : None for i2 in daughter2_indices}
            for i2,d2 in zip(daughter2_indices,daughters2):
                for i1,d1 in enumerate(daughters1_perm):
                    comp = d1.cmpProperties(d2,
                                            properties=['isSM','pdg'])
                    if comp == 0:
                        matches[i2] = d1
                        daughters1_perm.pop(i1)
                        break
                # If no matches were found for d2
                # go to the next daughters1 permutation
                if matches[i2] is None:
                    break
                
            # If all daughters from other have been matched
            # there is no need to look for other permutations
            if all(d_match is not None for d_match in matches.values()):
                break
        
        # If after going over all permutations
        # there are still unmatched daughers, return None
        if any(d_match is None for d_match in matches.values()):
            return None

        # Create new vertex following the indices in other:
        matchedVertex = VertexGraph()
        matchedVertex.add_node(mom1,0)
        for i2,d1 in matches.items():
            matchedVertex.add_node(d1,i2)
            if i2 != 0:
                matchedVertex.add_edge(0,i2)

        return matchedVertex
