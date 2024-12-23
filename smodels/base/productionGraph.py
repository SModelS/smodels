"""
.. module:: productionGraph
   :synopsis: This is a base class for describing production diagrams (hard scattering) as simple graphs.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


from smodels.base.exceptions import SModelSBaseError as SModelSError
from smodels.base.genericGraph import GenericGraph
from smodels.share.models.SMparticles import proton
from smodels.base.vertexGraph import VertexGraph
from smodels.base.inclusiveObjects import InclusiveValue
from smodels.base.particle import Particle,MultiParticle


# Used to construct generic z2-odd and z2-even particles:
anyBSM = Particle(label='anyBSM', isSM=False, isSelfConjugate=True, pdg = InclusiveValue())
anySM = Particle(label='anySM', isSM=True,  isSelfConjugate=True, pdg = InclusiveValue())
anyParticle = MultiParticle(label='*', particles=[anyBSM,anySM])


class ProductionGraph(GenericGraph):

    """
    A class for describing and manipulating production diagrams described as simple graphs.
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


def get2to1ProdGraph(model,initialStates=[proton,proton],finalState=anyBSM):
    """
    Get a list of 2 > 1 production graphs with the given initial states and the given
    final state using the vertices from the model.
    The initialStates and finalState arguments can be used to select the desired set of possible
    production diagrams.

    :param initialStates: List of 2 initial states (Particle or MultiParticle objects) for the hard scattering
    :param finalState: Final state (Particle or MultiParticle objects) for the hard scattering

    :return:
    """

    if len(initialStates) != 2:
        raise SModelSError(f"2to1 production graphs should have two initial states defined (instead of {initialStates})")
    
    if not isinstance(finalState,(Particle,MultiParticle)):
        raise SModelSError(f"Final state should be a Particle or a MultiParticle object and not {type(finalState)}")
    
    # Create the vertex defining the 2->1 process:
    prodV = VertexGraph(incoming=initialStates,outgoing=[finalState])

    # Get vertices from model
    modelVertices = model.vertices

    # Get matched vertices:
    matchedV = []
    for mVertex in modelVertices:
        matchV = mVertex.matchTo(prodV)
        if matchV is None:
            continue
        matchedV.append(matchV)

    return matchedV


def get2to2ProdGraph(model,initialStates=[proton,proton],
                     finalStates=[anyBSM,anyParticle],
                     channels=['s','t','u']):
    """
    Get a list of 2 > 2 production graphs with the given initial states and the given
    finals states using the vertices from the model.
    The initialStates, finalState and channels arguments can be used to select the desired set of possible
    production diagrams.

    :param initialStates: List of 2 initial states (Particle or MultiParticle objects) for the hard scattering
    :param finalStates: List of 2 final states (Particle or MultiParticle objects) for the hard scattering
    :param channels: List of strings specifying which channels to consider (possible channels are s,t, or u) 

    :return:
    """

    if len(initialStates) != 2:
        raise SModelSError(f"2to2 production graphs should have two initial states defined (instead of {initialStates})")
    if len(finalStates) != 2:
        raise SModelSError(f"2to2 production graphs should have two final states defined (instead of {finalStates})")


    matches = {}
    # Create the generic vertices defining the 2->2 s-channel process:
    if any('s' == channel for channel in channels):
        # Get vertices from model
        modelVertices = model.vertices
        inParticles = initialStates[:]
        outParticles = finalStates[:]
        prodL_schannel = VertexGraph(incoming=inParticles,outgoing=[anyParticle])
        matched_schannel = []
        for mVertexL in modelVertices:
            matched_left = mVertexL.matchTo(prodL_schannel)
            if matched_left is None:
                continue
            # Create the right vertex with the incoming particle given by
            #  same particle outgoing from the left vertex
            propParticle = matched_left.indexToNode(matched_left.outgoing_nodes[0])
            prodR_schannel = VertexGraph(incoming=[propParticle],outgoing=outParticles)
            for mVertexR in modelVertices:
                matched_right = mVertexR.matchTo(prodR_schannel)
                if matched_right is None:
                    continue
                matched_schannel.append((matched_left,matched_right))
        
        matches['s'] = matched_schannel

    # For a t-channel we just need to switch one initial state with one final state:
    if any('t' == channel for channel in channels):
        inParticles = [initialStates[0],finalStates[0]]
        outParticles = [initialStates[1],finalStates[1]]
        # With this change we can just recycle the s-channel calculation
        matches['t'] = get2to2ProdGraph(model,initialStates=inParticles,
                                        finalStates=outParticles,
                                        channels=['s'])['s']
    
    # For a u-channel we just need to another switch between one initial state and one final state:
    if any('u' == channel for channel in channels):
        inParticles = [initialStates[0],finalStates[1]]
        outParticles = [initialStates[1],finalStates[0]]
        # With this change we can just recycle the s-channel calculation
        matches['u'] = get2to2ProdGraph(model,initialStates=inParticles,
                                        finalStates=outParticles,
                                        channels=['s'])['s']

    return matches
    
