"""
.. module:: tree
   :synopsis: Classes used to construct trees (root directed graphs) describing
              simplified model topologies.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.theory.particle import Particle
from smodels.tools.inclusiveObjects import InclusiveValue

# Define a common inclusive particle object
# to be used with the inclusive node
IncluviseParticle = Particle(label='Inclusive')


class ParticleNode(object):
    """
    Simple wrapper for creating graphs with Particle objects.
    It is necessary because the same particle can appear multiple
    times within a tree, so
    Particle objects can not be directly used as nodes
    (since the same particle can not appear as distinct nodes)

    :ivar particle: Stores the Particle object
    :ivar nodeNumber: Node identifier
    :ivar nodeWeight: Stores the node weight
                      (1 for stable particles, BR for unstable
                      and production xsec for primary vertex)
    """

    def __init__(self, particle, nodeWeight=1.0,
                 canonName=None, isFinalState=False,
                 finalStates=None, isInclusive=False):
        self.particle = particle

        # Since ParticleNodes are identified by their numbering,
        # if it is not specifically assigned, automatically assign
        # a new number which does not overlap with any previous class instances
        self.canonName = canonName
        # Node weight:
        # (1 for stable particles, BR for unstable and
        # production xsec for primary vertex)
        self.nodeWeight = nodeWeight

        # Flag to tag nodes which should not be decayed
        self.isFinalState = isFinalState

        # Possibility to store final states generated by
        # the cascade decay of the node
        if isinstance(finalStates, list):
            self.finalStates = sorted(finalStates)
        else:
            self.finalStates = finalStates

        # Flag to identify as non-inclusive node:
        self.isInclusive = isInclusive

    def __hash__(self):
        return object.__hash__(self)

    def __cmp__(self, other):
        """
        Node comparison based on compareTo method

        :return: 0 if nodes are equal, 1 if self > other, -1 otherwise
        """

        return self.compareTo(other)

    def __lt__(self, other):
        return self.__cmp__(other) == -1

    def __gt__(self, other):
        return self.__cmp__(other) == 1

    def __eq__(self, other):
        return self.__cmp__(other) == 0

    def __ne__(self, other):
        return self.__cmp__(other) != 0

    def __str__(self):

        return str(self.particle)

    def __repr__(self):

        return self.__str__()

    def __getattr__(self, attr):
        """
        Returns the attribute from particle.

        :parameter attr: Attribute string

        :return: Attribute from self.particle
        """

        return getattr(self.particle, attr)

    def __getstate__(self):
        """
        Since getattr is defined, we must defined getstate
        for pickling/unpickling.
        """

        attrDict = self.__dict__

        return attrDict

    def __setstate__(self, state):
        """
        Since getattr is defined, we must defined getstate
        for pickling/unpickling.
        """

        self.__dict__.update(state)

    def __add__(self, other):
        """
        Adds two nodes. The properties of self are kept, except
        for the particle and nodeWeight, which are added with other.

        :param other: ParticleNode object

        :return: a copy of self with the particle combined with other.particle
                 and nodeWeight added.
        """

        newNode = self.copy()
        newNode.particle = self.particle + other.particle
        newNode.nodeWeight = self.nodeWeight + other.nodeWeight

        return newNode

    def compareTo(self, other):
        """
        Compare nodes accoring to  particles.

        :param other: ParticleNode or InclusiveParticleNode object

        :return: 1 if self > other, -1 if self < other and 0 if self == other
        """

        if other.isInclusive:
            return -other.compareTo(self)

        if not isinstance(other, ParticleNode):
            raise SModelSError("Can not compare node to %s" % type(other))

        if self.particle != other.particle:
            if self.particle > other.particle:
                return 1
            else:
                return -1

        return 0

    def equalTo(self, other):
        """
        Compare nodes accoring to their canonical name
        and particle.

        :param other: ParticleNode or InclusiveParticleNode object

        :return: True if nodes are equal, false otherwise
        """

        return (self.compareTo(other) == 0)

    def copy(self):
        """
        Makes a shallow copy of itself. The particle attribute
        shares the same object with the original copy.
        :return: ParticleNode object
        """

        newNode = ParticleNode(particle=self.particle,
                               canonName=self.canonName,
                               isInclusive=self.isInclusive,
                               isFinalState=self.isFinalState)

        if self.finalStates is not None:
            newNode.finalStates = self.finalStates[:]
        if isinstance(self.nodeWeight,(int,float,type(None))):
            newNode.nodeWeight = self.nodeWeight
        else:
            newNode.nodeWeight = self.nodeWeight.copy()

        return newNode


class InclusiveParticleNode(ParticleNode):
    """
    An inclusive ParticleNode class. It will return True when compared to any other
    ParticleNode object or InclusiveParticleNode object.

    :ivar particle: IncluviseParticle (dummy)
    :ivar nodeNumber: Node identifier
    :ivar nodeWeight: 1 (dummy value)
    :ivar finalStates: Allowed final states (final state nodes)
    """

    def __init__(self, particle=IncluviseParticle,
                 nodeNumber=None, nodeWeight=1.0,
                 canonName=InclusiveValue(), isFinalState=True,
                 finalStates=[], isInclusive=True):

        ParticleNode.__init__(self, particle=particle,
                              nodeNumber=nodeNumber, nodeWeight=nodeWeight,
                              canonName=canonName, isFinalState=isFinalState,
                              finalStates=finalStates, isInclusive=isInclusive)

    def compareTo(self, other, inclusive=True):
        """
        Compares only the finalStates attributes of self and other.
        All the final states in other have to match at least one final
        state in self.

        :param other: ParticleNode or InclusiveParticleNode object
        :param inclusive: If False, particles are required to be identical
                          (the inclusiveness of MultiParticles or InclusiveNodes
                          are not considered when comparing)

        :return: 0 (self == other), 1 (self > other), -1 (self < other)
        """

        fsA = self.finalStates[:]  # Sorted final states in A
        fsB = other.finalStates[:]  # Sorted final states in B

        # If inclusive=False, requires final states to be the same particle
        # and other node to be an InclusiveNode
        if not inclusive:
            if not other.isInclusive:
                return -1
            if len(fsA) != len(fsB):
                if len(fsA) > len(fsB):
                    return 1
                else:
                    return -1
            for iB, fs in enumerate(fsB):
                if fs is fsA[iB]:
                    continue
                if fs < fsA[iB]:
                    return 1
                else:
                    return -1
            return 0

        for fs in fsB:
            if fs in fsA:
                continue
            elif not other.isInclusive:
                return -1  # Always smaller than ParitcleNode
            elif fs > fsA[-1]:  # Define other > self if fs > largest state in self
                return -1
            else:
                return 1
        return 0

    def copy(self):
        """
        Makes a shallow copy of itself. The particle attribute
        shares the same object with the original copy.
        :return: ParticleNode object
        """

        newNode = InclusiveParticleNode(particle=self.particle,
                                        nodeWeight=self.nodeWeight,
                                        canonName=self.canonName,
                                        isFinalState=self.isFinalState,
                                        finalStates=self.finalStates,
                                        isInclusive=self.isInclusive)

        return newNode

    def longStr(self):
        """
        Returns a string representation of self containing
        its final states.

        :return: String
        """

        nodeStr = str(self)
        if self.finalStates:
            fStates = [str(ptc) for ptc in self.finalStates]
            nodeStr += '\n(%s)' % (','.join(fStates))
        return nodeStr

    def __getattr__(self, attr):
        """
        Returns None if does not contain the attribute.

        :parameter attr: Attribute string

        :return: Attribute value or None
        """
        if attr not in self.__dict__:
            return None