"""
.. module:: particle
   :synopsis: Defines the particle, multiparticle and particle list classes, their methods and related functions

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

import weakref
from typing import Any
from collections.abc import Iterable


class Particle(object):
    """
    An instance of this class represents a single particle.
    The properties are: label, pdg, mass, electric charge, color charge, width
    """

    _instances : set = set()
    _lastID : int = 0
    _attr_cache : dict = {}  # frozen attribute tuple -> existing Particle instance

    @classmethod
    def _frozen_key(cls, attrDict: dict) -> tuple:
        """Build a hashable key from the normalized attribute dict."""
        def _hashable(v):
            try:
                hash(v)
                return v
            except TypeError:
                return repr(v)
        return tuple(sorted((k, _hashable(v)) for k, v in attrDict.items()))

    
    def __new__(cls, attributesDict: dict = {}, **kwargs: Any) -> 'Particle':
        """
        Creates a particle. If a particle with the exact same attributes have
        already been created return this particle instead.
        Assigns an ID to the instance using the class Particle._instance
        list. Reset the comparison dictionary.

        :param attributesDict: A dictionary with particle attributes (useful for pickling/unpickling).
                               Attributes can also be directly assigned using keyword arguments.

        Possible properties for arguments.
        isSM: True/False
        label: str, e.g. 'e-'
        pdg: number in pdg
        mass: mass of the particle
        echarge: electric charge as multiples of the unit charge
        colordim: color dimension of the particle
        spin: spin of the particle
        totalwidth: total width
        """

        if not kwargs and not attributesDict:
            raise ValueError("Particle object can not be created with empty attributes")

        attrDict = dict(attributesDict.items())
        attrDict.update(kwargs)
        attrDict = cls._normalizedAttrDict(attrDict)
        key = cls._frozen_key(attrDict)
        cached = cls._attr_cache.get(key)
        if cached is not None:
            return cached

        newParticle = super(Particle, cls).__new__(cls)
        for attr, value in attrDict.items():
            setattr(newParticle, attr, value)
        newParticle._id = Particle.getID()
        newParticle._comp = {newParticle._id: 0}
        newParticle._isStable = None
        newParticle._isPrompt = None
        cls._instances.add(weakref.ref(newParticle))
        cls._attr_cache[key] = newParticle
        return newParticle

    def __getnewargs__(self) -> tuple:
        """
        Required for unpickling the object.
        When loading the pickled object, it will call __new__ with the
        arguments returned by this method.
        """

        attrDict = self.__class__._normalizedAttrDict(self.__dict__)
        return (attrDict,)

    def __getstate__(self) -> dict:
        """
        Makes sure the particle ID and comparison matrix
        are not stored in the pickle file, so they are dynamically assigned
        when the pickle file is loaded.
        """

        attrDict = self.__class__._normalizedAttrDict(self.__dict__)

        return attrDict

    def __setstate__(self, state: dict):
        """
        Dummy function, since all the initialization and attribute
        setting is handled by __new__.
        """
        pass

    def __hash__(self) -> int:
        """
        Return the object address. Required for using weakref
        """
        return self._id
    
    @classmethod
    def _normalizedAttrDict(cls, data: dict) -> dict:
        """Drop runtime-only cache fields from attribute comparisons."""

        attrDict = dict(data.items())
        attrDict.pop('_id', None)
        attrDict.pop('_comp', None)
        attrDict.pop('_isStable', None)
        attrDict.pop('_isPrompt', None)
        return attrDict


    @classmethod
    def getinstances(cls) -> list:
        """Return all live instances of this Particle class (pruning dead weakrefs)."""
        dead = set()
        instances = []
        for ref in Particle._instances:
            obj = ref()
            if obj is not None:
                instances.append(obj)
            else:
                dead.add(ref)
        Particle._instances -= dead
        return instances

    @classmethod
    def getID(cls) -> int:
        """Return a unique integer ID for a new particle instance."""
        Particle._lastID += 1
        return Particle._lastID

    def __cmp__(self, other: 'Particle') -> int:
        """
        Compares particle with other.
        The comparison is based on the particle properties.

        :param other:  particle to be compared (Particle or MultiParticle object)

        :return: -1 if particle < other, 1 if particle > other and 0 if particle == other
        """

        # First check if we have already compared to this object
        if other._id in self._comp:
            return self._comp[other._id]
        elif self._id in other._comp:
            return -other._comp[self._id]

        cmpProp = self.cmpProperties(other)  # Objects have not been compared yet.
        self._comp[other._id] = cmpProp
        other._comp[self._id] = -cmpProp
        return cmpProp

    def __lt__(self, p2: 'Particle') -> bool:
        return self.__cmp__(p2) == -1

    def __gt__(self, p2: 'Particle') -> bool:
        return self.__cmp__(p2) == 1

    def __eq__(self, p2: object) -> bool:
        return self.__cmp__(p2) == 0

    def __ne__(self, p2: object) -> bool:
        return self.__cmp__(p2) != 0

    def __str__(self) -> str:
        if hasattr(self, 'label'):
            return self.label
        else:
            return ''

    def __repr__(self) -> str:
        return self.__str__()

    def __add__(self, other: 'Particle') -> 'Particle':
        """
        Define addition of two Particle objects
        or a Particle object and a MultiParticle object.
        The result is a MultiParticle object containing
        both particles.
        """

        if isinstance(other, MultiParticle):
            return other.__add__(self)
        elif self.contains(other):
            return self
        elif other.contains(self):
            return other
        else:
            combined = MultiParticle(particles=[self, other])
            return combined

    def __radd__(self, other: 'Particle') -> 'Particle':
        return self.__add__(other)

    def __iadd__(self, other: 'Particle') -> 'Particle':
        return self.__add__(other)

    def describe(self) -> str:
        """Return a string representation of all particle attributes."""
        return str(self.__dict__)

    def eqProperties(self, other: 'Particle',
                     properties: list[str] = ['isSM', 'spin', 'colordim', 'eCharge', 'mass', 'totalwidth']) -> bool:
        """
        Check if particle has the same properties (default is spin, colordim and eCharge)
        as other. Only compares the attributes which have been defined in both objects.

        :param other: a Particle or MultiParticle object
        :param properties: list with properties to be compared. Default is spin, colordim and eCharge

        :return: True if all properties are the same, False otherwise.
        """

        if self.cmpProperties(other, properties=properties) == 0:
            return True
        else:
            return False

    def cmpProperties(self, other: 'Particle',
                      properties: list[str] = ['isSM', 'spin', 'colordim', 'eCharge', 'mass', 'totalwidth']) -> int:
        """
        Compare properties (default is isSM, spin, colordim, eCharge, mass and totalwidth).
        Return 0 if properties are equal, -1 if self < other and 1 if self > other.
        Only compares the attributes which have been defined in both objects.
        The comparison is made in hierarchical order, following the order
        defined by the properties list.

        :param other: a Particle or MultiParticle object
        :param properties: list with properties to be compared. Default is spin, colordim and eCharge

        :return: 0 if properties are equal, -1 if self < other and 1 if self > other.
        """

        if isinstance(other, (MultiParticle)):
            return -1 * other.cmpProperties(self, properties=properties)

        for prop in properties:
            if not hasattr(self, prop) or not hasattr(other, prop):
                continue
            x = getattr(self, prop)
            y = getattr(other, prop)
            if x == y:
                continue
            if x > y:
                return 1
            else:
                return -1

        return 0

    def copy(self) -> 'Particle':
        """
        Make a copy of self with a distinct ID.

        :return: A Particle object identical to self, except for its ID and comparison dict
        """

        cls = self.__class__
        newParticle = object.__new__(cls)
        for attr, value in self.__dict__.items():
            setattr(newParticle, attr, value)
        newParticle._id = cls.getID()
        newParticle._comp = {newParticle._id: 0}
        cls._instances.add(weakref.ref(newParticle))
        # copy() intentionally creates a NEW particle (different ID),
        # so do NOT add to _attr_cache.

        return newParticle

    def chargeConjugate(self, label: str | None = None) -> 'Particle':
        """
        Returns the charge conjugate particle (flips the sign of eCharge).
        If it has a pdg property also flips its sign.
        If label is None, the charge conjugate name is defined as the original name plus "~" or
        if the original name ends in "+" ("-"), it is replaced by "-" ("+").

        :parameter label: If defined, defines the label of the charge conjugated particle.

        :return: the charge conjugate particle (Particle object)
        """

        particleAttr = dict(self.__dict__.items())
        for attr, value in particleAttr.items():
            if attr in ['pdg', 'eCharge'] and isinstance(value, (float, int)):
                particleAttr[attr] = -1*value
            if attr == 'label':
                if value[-1] == '+':
                    particleAttr[attr] = value[:-1]+'-'
                elif value[-1] == '-':
                    particleAttr[attr] = value[:-1]+'+'
                elif value[-1] == '~':
                    particleAttr[attr] = value[:-1]
                else:
                    particleAttr[attr] = value+'~'

        # Overwrite default labelling
        if label is not None:
            particleAttr['label'] = label
        pConjugate = self.__class__(**particleAttr)

        return pConjugate

    def isNeutral(self) -> bool:
        """
        Return True if the particle is electrically charged and color neutral.
        If these properties have not been defined, return True.

        :return: True/False
        """

        if hasattr(self, 'eCharge') and self.eCharge != 0:
            return False
        if hasattr(self, 'colordim') and self.colordim != 1:
            return False

        return True

    def isMET(self) -> bool:
        """
        Checks if the particle can be considered as MET.
        If the _isInvisible attribute has not been defined, it will return True/False is isNeutral() = True/False.
        Else it will return the _isInvisible attribute.

        :return: True/False
        """

        if hasattr(self, '_isInvisible'):
            return self._isInvisible
        else:
            return self.isNeutral()

    def isPrompt(self) -> bool:
        """
        Checks if the particle decays promptly.
        If _isPrompt has been set, return its value,
        otherwise set to True if self.totalwidth == inf.

        :return: True/False
        """
        if self._isPrompt is None:
            self._isPrompt = (self.totalwidth.asNumber() == float('inf'))

        return self._isPrompt

    def isStable(self) -> bool:
        """
        Checks if the particle is stable.
        If _isStable has been set, return its value,
        otherwise set to True if self.totalwidth == 0.

        :return: True/False
        """

        if self._isStable is None:
            self._isStable = (self.totalwidth.asNumber() == 0)

        return self._isStable

    def contains(self, particle: 'Particle') -> bool:
        """
        If particle is a Particle object check if self and particle are the same object.

        :param particle: Particle or MultiParticle object

        :return: True/False
        """

        if self is particle:
            return True
        else:
            return False


class InvisibleParticle(Particle):

    """
    An instance of this class represents an invisible particle created by invisible compression,
    so it does not correspond to a physical BSM particle. 
    It behaves like Particle, but object de-duplication and ID generation
    are restricted to InvisibleParticle (or subclasses of it), without looking at
    Particle._instances, which improves performance.
    """

    _instances : set = set()
    _attr_cache : dict = {}  # separate cache for InvisibleParticle
    
    def __new__(cls, attributesDict: dict = {}, **kwargs: Any) -> 'InvisibleParticle':
        """Create an InvisibleParticle, enforcing invisibility on all instances."""
        attrDict = dict(attributesDict.items())
        attrDict.update(kwargs)
        # Enforce invisibility for all InvisibleParticle instances.
        attrDict['_isInvisible'] = True
        attrDict = cls._normalizedAttrDict(attrDict)
        key = cls._frozen_key(attrDict)
        cached = cls._attr_cache.get(key)
        if cached is not None:
            return cached

        newParticle = super(Particle, cls).__new__(cls)
        for attr, value in attrDict.items():
            setattr(newParticle, attr, value)
        newParticle._id = cls.getID()
        newParticle._comp = {newParticle._id: 0}
        newParticle._isStable = None
        newParticle._isPrompt = None
        cls._instances.add(weakref.ref(newParticle))
        cls._attr_cache[key] = newParticle
        return newParticle

    @classmethod
    def getinstances(cls) -> list:
        """Return all live instances of this InvisibleParticle class."""
        dead = set()
        instances = []
        for ref in InvisibleParticle._instances:
            obj = ref()
            if obj is not None:
                instances.append(obj)
            else:
                dead.add(ref)
        InvisibleParticle._instances -= dead
        return instances




class MultiParticle(Particle):

    """ An instance of this class represents a list of particle object to allow for inclusive expressions such as jets.
        The properties are: label, pdg, mass, electric charge, color charge, width
    """

    _instances : set = set()

    def __new__(cls, label: str | None = None, particles: list = [], attributesDict: dict = {}, **kwargs: Any) -> 'MultiParticle':
        """
        Creates a multiparticle. If a multiparticle with the exact same particles
        already been created return this multiparticle instead.
        Assigns an ID to the isntance using the class Particle._instance
        list. Reset the comparison dictionary.

        :param label: Label for the MultiParticle (string)
        :param particles: List of Particle or MultiParticle objects (list)
        :param attributesDict: A dictionary with particle attributes (useful for pickling/unpickling).
                               Attributes can also be directly assigned using keyword arguments.
        """

        particles = sorted(particles)
        if not label:
            label = "/".join(sorted(list(set([p.label for p in particles]))))
        attrDict = dict(attributesDict.items())
        attrDict.update(kwargs)
        attrDict = cls._normalizedAttrDict(attrDict)
        for obj in MultiParticle.getinstances()[:]:
            # Directly compare attributes, except for particles,label,id and _comp
            objAttr = cls._normalizedAttrDict(obj.__dict__)
            objAttr.pop('label', None)
            objAttr.pop('particles', None)
            if objAttr != attrDict:
                continue
            pListB = obj.particles
            if len(particles) != len(pListB):
                continue
            if any(pA is not pListB[i] for i, pA in enumerate(particles)):
                continue
            return obj

        newMultiParticle = super(Particle, cls).__new__(cls)
        for attr, value in attrDict.items():
            setattr(newMultiParticle, attr, value)
        newMultiParticle.particles = particles[:]
        newMultiParticle.label = label
        newMultiParticle._id = Particle.getID()
        newMultiParticle._comp = {newMultiParticle._id: 0}
        newMultiParticle._comp.update(dict([[ptc._id, 0] for ptc in particles]))
        MultiParticle._instances.add(weakref.ref(newMultiParticle))
        return newMultiParticle

    @classmethod
    def getinstances(cls) -> list:
        """Return all live instances of this MultiParticle class."""
        dead = set()
        instances = []
        for ref in MultiParticle._instances:
            obj = ref()
            if obj is not None:
                instances.append(obj)
            else:
                dead.add(ref)
        MultiParticle._instances -= dead

        return instances

    def __getnewargs__(self) -> tuple:
        """
        Required for unpickling the object.
        When loading the pickled object, it will call __new__ with the
        arguments returned by this method.
        """

        attrDict = self.__class__._normalizedAttrDict(self.__dict__)
        attrDict.pop('label', None)
        attrDict.pop('particles', None)

        return (self.label, self.particles, attrDict)

    def __getstate__(self) -> dict:
        """
        Makes sure the particle ID and comparison matrix
        are not stored in the pickle file, so they are dynamically assigned
        when the pickle file is loaded.
        """

        attrDict = self.__class__._normalizedAttrDict(self.__dict__)

        return attrDict

    def __setstate__(self, state: dict):
        """
        Dummy function, since all the initialization and attribute
        setting is handled by __new__.
        """
        pass

    def __getattr__(self, attr: str):
        """
        If MultiParticle does not have attribute, return a list
        if the attributes of each particle in self.particles.
        If not all particles have the attribute, it will raise an exception.
        If all particles share a common attribute, a single value
        will be returned.

        :parameter attr: Attribute string

        :return: Attribute or list with the attribute values in self.particles
        """

        try:
            values = [getattr(particle, attr) for particle in self.particles]
            if all(type(x) == type(values[0]) for x in values):
                if all(x == values[0] for x in values):
                    return values[0]
            return values
        except (AttributeError, IndexError, TypeError) as e:
            raise AttributeError(e)

    def cmpProperties(self, other: 'Particle',
                      properties: list[str] = ['isSM', 'spin', 'colordim', 'eCharge', 'mass', 'totalwidth']) -> int:
        """
        Compares the properties in self with the ones in other.
        If other is a Particle object, checks if any of the particles in self matches
        other.
        If other is a MultiParticle object, checks if any particle in self matches
        any particle in other.
        If self and other are equal returns 0,
        else returns the result of comparing the first particle of self with other.

        :param other: a Particle or MultiParticle object
        :param properties: list with properties to be compared. Default is spin, colordim and eCharge

        :return: 0 if properties are equal, -1 if self < other and 1 if self > other.
        """

        # Check if any of its particles matches other
        if isinstance(other, Particle):
            otherParticles = [other]
        elif isinstance(other, MultiParticle):
            otherParticles = other.particles

        for otherParticle in otherParticles:
            if any(p.cmpProperties(otherParticle, properties) == 0 for p in self.particles):
                return 0

        cmpv = self.particles[0].cmpProperties(otherParticles[0], properties)
        return cmpv

    def __add__(self, other: 'Particle') -> 'Particle':
        """
        Define addition of two Particle objects
        or a Particle object and a MultiParticle object.
        The result is a MultiParticle object containing
        both particles.
        """

        if other is self or self.contains(other):  # Check if other is self or a subset of self
            return self
        # Check if self is a subset of other
        if other.contains(self):
            return other
        elif isinstance(other, MultiParticle):
            addParticles = [ptc for ptc in other.particles if not self.contains(ptc)]
        elif isinstance(other, Particle):
            addParticles = [other]

        combinedParticles = self.particles + addParticles
        combined = MultiParticle(particles=combinedParticles)

        return combined

    def __radd__(self, other: 'Particle') -> 'Particle':
        return self.__add__(other)

    def __iadd__(self, other: 'Particle') -> 'Particle':
        return self.__add__(other)

    def getPdgs(self) -> list:
        """
        pdgs appearing in MultiParticle
        :return: list of pgds of particles in the MultiParticle
        """

        pdgs = [particle.pdg for particle in self.particles]

        return pdgs

    def getLabels(self) -> list:
        """
        labels appearing in MultiParticle
        :return: list of labels of particles in the MultiParticle
        """

        labels = [particle.label for particle in self.particles]

        return labels

    def isNeutral(self) -> bool:
        """
        Return True if ALL particles in particle list are neutral.

        :return: True/False
        """

        neutral = all(particle.isNeutral() for particle in self.particles)
        return neutral

    def isMET(self) -> bool:
        """
        Checks if all the particles in self can be considered as MET.

        :return: True/False
        """

        met = all(particle.isMET() for particle in self.particles)
        return met

    def contains(self, particle: 'Particle') -> bool:
        """
        Check if MultiParticle contains the Particle object or MultiParticle object.

        :param particle: Particle or MultiParticle object

        :return: True/False
        """

        if isinstance(particle, MultiParticle):
            checkParticles = particle.particles
        else:
            checkParticles = [particle]

        for otherParticle in checkParticles:
            hasParticle = False
            for p in self.particles:
                if p.contains(otherParticle):
                    hasParticle = True
            if not hasParticle:
                return False

        return True
