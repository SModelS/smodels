"""
.. module:: vertexGraph
   :synopsis: This is a base class for describing vertices as simple graphs.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


from smodels.base.exceptions import SModelSBaseError as SModelSError
from smodels.base.genericGraph import GenericGraph
from collections import OrderedDict
from itertools import product


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
        ## Sort all particles by BSM and pdg values
        all_particles = sorted(all_particles, 
                              key = lambda p: (not p.isSM,abs(p.pdg)), reverse=True)
        if len(all_particles) > 0:
            in_particle = all_particles[0]
            out_particles = [p.chargeConjugate() for p in all_particles[1:]]

            # Now if in_particle has a negative PDG, conjugate the vertex:
            if in_particle.pdg < 0:
                in_particle = in_particle.chargeConjugate()
                out_particles = [p.chargeConjugate() for p in out_particles[:]]

            # Finally sort outgoing particles again
            out_particles = sorted(out_particles, 
                                    key = lambda p: (not p.isSM,abs(p.pdg)), 
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
    
