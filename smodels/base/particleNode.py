"""
.. module:: tree
   :synopsis: Classes used to construct trees (root directed graphs) describing
              simplified model topologies.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


from smodels.base.exceptions import SModelSBaseError as SModelSError
from smodels.base.particle import Particle
from smodels.base.inclusiveObjects import InclusiveValue

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
    """

    def __init__(self, particle,
                 isFinalState=False,
                 isInclusive=False, inclusiveList=False):

        self.particle = particle

        # Flag to tag nodes which should not be decayed
        self.isFinalState = isFinalState

        # Flag to identify as non-inclusive node:
        self.isInclusive = isInclusive

        # Flag to identify as non-inclusive node:
        self.inclusiveList = inclusiveList

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

        label = str(self.particle)
        if self.inclusiveList:
            label = '*'+label
        return label

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
        for the particle, which are added with other.

        :param other: ParticleNode object

        :return: a copy of self with the particle combined with other.particle
        """

        newNode = self.copy()
        newNode.particle = self.particle + other.particle

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
            raise SModelSError(f"Can not compare node to {type(other)}")

        cmp = self.particle.__cmp__(other.particle)

        return cmp

    def equalTo(self, other):
        """
        Compare nodes accoring to their particle.

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
                               isInclusive=self.isInclusive,
                               isFinalState=self.isFinalState)

        return newNode


class InclusiveParticleNode(ParticleNode):
    """
    An inclusive ParticleNode class. It will return True when compared to any other ParticleNode object or InclusiveParticleNode object.

    :ivar particle: IncluviseParticle (dummy)
    """

    def __init__(self, particle=IncluviseParticle):

        ParticleNode.__init__(self, particle=particle,
                              isInclusive=True)

    def compareTo(self, other):
        """
        Dummy method. Always return 0, since it will always match any node.

        :param other: ParticleNode or InclusiveParticleNode object

        :return: 0
        """

        return 0

    def copy(self):
        """
        Makes a shallow copy of itself. The particle attribute
        shares the same object with the original copy.
        :return: ParticleNode object
        """

        newNode = InclusiveParticleNode(particle=self.particle)

        return newNode

    def __getattr__(self, attr):
        """
        Returns None if does not contain the attribute.

        :parameter attr: Attribute string

        :return: Attribute value or None
        """

        if attr not in self.__dict__:
            return None
