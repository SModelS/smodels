"""
.. module:: particle
   :synopsis: Defines the particle, multiparticle and particle list classes, their methods and related functions

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

import itertools,weakref

class Particle(object):
    """
    An instance of this class represents a single particle.
    The properties are: label, pdg, mass, electric charge, color charge, width
    """

    _instances = set()
    _lastID = 0

    def __new__(cls,attributesDict={}, **kwargs):
        """
        Creates a particle. If a particle with the exact same attributes have
        already been created return this particle instead.
        Assigns an ID to the instance using the class Particle._instance
        list. Reset the comparison dictionary.

        :param attributesDict: A dictionary with particle attributes (useful for pickling/unpickling).
                               Attributes can also be directly assigned using keyword arguments.

        Possible properties for arguments.
        Z2parity: int, +1 or -1
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
        attrDict.pop('_id',None)
        attrDict.pop('_comp',None)
        for obj in Particle.getinstances():
            if not isinstance(obj,Particle):
                continue
            objAttr = dict(obj.__dict__.items())
            objAttr.pop('_id',None)
            objAttr.pop('_comp',None)
            if objAttr != attrDict:
                continue
            return obj

        newParticle = super(Particle, cls).__new__(cls)
        for attr,value in attrDict.items():
            setattr(newParticle,attr,value)
        newParticle._id = Particle.getID()
        newParticle._comp = {newParticle._id : 0}
        Particle._instances.add(weakref.ref(newParticle))
        return newParticle

    def __getnewargs__(self):
        """
        Required for unpickling the object.
        When loading the pickled object, it will call __new__ with the
        arguments returned by this method.
        """

        attrDict = dict(self.__dict__.items())
        #Make sure pickled/unpickled objects do no store ID nor comparison dict
        attrDict.pop('_id',None)
        attrDict.pop('_comp',None)
        return (attrDict,)

    def __getstate__(self):
        """
        Makes sure the particle ID and comparison matrix
        are not stored in the pickle file, so they are dynamically assigned
        when the pickle file is loaded.
        """

        attrDict = dict(self.__dict__.items())
        #Make sure pickled objects do no store ID nor comparison dict
        attrDict.pop('_id',None)
        attrDict.pop('_comp',None)

        return attrDict

    def __setstate__(self,state):
        """
        Dummy function, since all the initialization and attribute
        setting is handled by __new__.
        """
        pass

    def __hash__(self):
        """
        Return the object address. Required for using weakref
        """
        return self._id

    @classmethod
    def getinstances(cls):
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
    def getID(cls):
        if len(Particle.getinstances()) == 0:
            Particle._lastID = 0
        else:
            Particle._lastID += 1
        return Particle._lastID

    def __cmp__(self,other):
        """
        Compares particle with other.
        The comparison is based on the particle properties.

        :param other:  particle to be compared (Particle or MultiParticle object)

        :return: -1 if particle < other, 1 if particle > other and 0 if particle == other
        """

        #First check if we have already compared to this object
        if other._id in self._comp:
            return self._comp[other._id]
        elif self._id in other._comp:
            return -other._comp[self._id]

        cmpProp = self.cmpProperties(other) #Objects have not been compared yet.
        self._comp[other._id] = cmpProp
        other._comp[self._id] = -cmpProp
        return cmpProp

    def __lt__( self, p2 ):
        return self.__cmp__(p2) == -1

    def __gt__( self, p2 ):
        return self.__cmp__(p2) == 1

    def __eq__( self, p2 ):
        return self.__cmp__(p2) == 0

    def __ne__( self, p2 ):
        return self.__cmp__(p2) != 0

    def __str__(self):
        if hasattr(self, 'label'):
            return self.label
        else: return ''

    def __repr__(self):
        return self.__str__()

    def __add__(self, other):
        """
        Define addition of two Particle objects
        or a Particle object and a MultiParticle object.
        The result is a MultiParticle object containing
        both particles.
        """

        if isinstance(other,MultiParticle):
            return other.__add__(self)
        elif self.contains(other):
            return self
        elif other.contains(self):
            return other
        else:
            combined = MultiParticle(particles= [self,other])
            return combined

    def __radd__(self,other):
        return self.__add__(other)

    def __iadd__(self,other):
        return self.__add__(other)


    def describe(self):
        return str(self.__dict__)

    def eqProperties(self,other,
                     properties = ['Z2parity','spin','colordim','eCharge','mass','totalwidth']):
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

    def cmpProperties(self,other,
                      properties = ['Z2parity','spin','colordim','eCharge','mass','totalwidth']):
        """
        Compare properties (default is spin, colordim and eCharge).
        Return 0 if properties are equal, -1 if self < other and 1 if self > other.
        Only compares the attributes which have been defined in both objects.
        The comparison is made in hierarchical order, following the order
        defined by the properties list.

        :param other: a Particle or MultiParticle object
        :param properties: list with properties to be compared. Default is spin, colordim and eCharge

        :return: 0 if properties are equal, -1 if self < other and 1 if self > other.
        """

        if isinstance(other,(MultiParticle)):
            return -1*other.cmpProperties(self,properties=properties)

        for prop in properties:
            if not hasattr(self,prop) or not hasattr(other,prop):
                continue
            x = getattr(self,prop)
            y = getattr(other,prop)
            if x == y:
                continue
            if x > y:
                return 1
            else:
                return -1

        return 0

    def copy(self):
        """
        Make a copy of self with a distinct ID.

        :return: A Particle object identical to self, except for its ID and comparison dict
        """

        newParticle = object.__new__(Particle)
        for attr,value in self.__dict__.items():
            setattr(newParticle,attr,value)
        newParticle._id = Particle.getID()
        newParticle._comp = {newParticle._id : 0}
        Particle._instances.add(weakref.ref(newParticle))

        return newParticle

    def chargeConjugate(self,label=None):
        """
        Returns the charge conjugate particle (flips the sign of eCharge).
        If it has a pdg property also flips its sign.
        If label is None, the charge conjugate name is defined as the original name plus "~" or
        if the original name ends in "+" ("-"), it is replaced by "-" ("+").

        :parameter label: If defined, defines the label of the charge conjugated particle.

        :return: the charge conjugate particle (Particle object)
        """

        particleAttr = dict(self.__dict__.items())
        for attr,value in particleAttr.items():
            if attr in ['pdg','eCharge'] and isinstance(value,(float,int)):
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

        #Overwrite default labelling
        if label is not None:
            particleAttr['label'] = label
        pConjugate = Particle(**particleAttr)

        return pConjugate

    def isNeutral(self):
        """
        Return True if the particle is electrically charged and color neutral.
        If these properties have not been defined, return True.

        :return: True/False
        """

        if hasattr(self,'eCharge') and self.eCharge != 0:
            return False
        if hasattr(self,'colordim') and self.colordim != 1:
            return False

        return True

    def isMET(self):
        """
        Checks if the particle can be considered as MET.
        If the _isInvisible attribute has not been defined, it will return True/False is isNeutral() = True/False.
        Else it will return the _isInvisible attribute.

        :return: True/False
        """

        if hasattr(self,'_isInvisible'):
            return self._isInvisible
        else:
            return self.isNeutral()

    def isPrompt(self):
        """
        Checks if the particle decays promptly (e.g. totalwidth = inf).

        :return: True/False
        """

        return self.totalwidth.asNumber() == float('inf')

    def isStable(self):
        """
        Checks if the particle is stable (e.g. totalwidth = 0).

        :return: True/False
        """

        return self.totalwidth.asNumber() == 0.

    def contains(self,particle):
        """
        If particle is a Particle object check if self and particle are the same object.

        :param particle: Particle or MultiParticle object

        :return: True/False
        """

        if self is particle:
            return True
        else:
            return False


class MultiParticle(Particle):

    """ An instance of this class represents a list of particle object to allow for inclusive expressions such as jets.
        The properties are: label, pdg, mass, electric charge, color charge, width
    """

    def __new__(cls,label=None,particles=[],attributesDict={},**kwargs):
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
            label = "/".join([p.label for p in particles])
        attrDict = dict(attributesDict.items())
        attrDict.update(kwargs)
        attrDict.pop('_id',None)
        attrDict.pop('_comp',None)
        for obj in Particle.getinstances()[:]:
            if not isinstance(obj,MultiParticle):
                continue
            #Directly compare attributes, except for particles,label,id and _comp
            objAttr = dict(obj.__dict__.items())
            objAttr.pop('_id',None)
            objAttr.pop('_comp',None)
            objAttr.pop('label',None)
            objAttr.pop('particles',None)
            if objAttr != attrDict:
                continue
            pListB = obj.particles
            if len(particles) != len(pListB):
                continue
            if any(pA is not pListB[i] for i,pA in enumerate(particles)):
                continue
            return obj

        newMultiParticle = super(Particle, cls).__new__(cls)
        for attr,value in attrDict.items():
            setattr(newMultiParticle,attr,value)
        newMultiParticle.particles = particles[:]
        newMultiParticle.label = label
        newMultiParticle._id = Particle.getID()
        newMultiParticle._comp = {newMultiParticle._id : 0}
        newMultiParticle._comp.update(dict([[ptc._id,0] for ptc in particles]))
        Particle._instances.add(weakref.ref(newMultiParticle))
        return newMultiParticle

    def __getnewargs__(self):
        """
        Required for unpickling the object.
        When loading the pickled object, it will call __new__ with the
        arguments returned by this method.
        """

        attrDict = dict(self.__dict__.items())
        attrDict.pop('label',None)
        attrDict.pop('particles',None)
        #Make sure pickled/unpickled objects do no store ID nor comparison dict
        attrDict.pop('_id',None)
        attrDict.pop('_comp',None)

        return (self.label,self.particles,attrDict)

    def __getstate__(self):
        """
        Makes sure the particle ID and comparison matrix
        are not stored in the pickle file, so they are dynamically assigned
        when the pickle file is loaded.
        """

        attrDict = dict(self.__dict__.items())
        #Make sure pickled objects do no store ID nor comparison dict
        attrDict.pop('_id',None)
        attrDict.pop('_comp',None)

        return attrDict

    def __setstate__(self,state):
        """
        Dummy function, since all the initialization and attribute
        setting is handled by __new__.
        """
        pass

    def __getattr__(self,attr):
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
            values = [getattr(particle,attr) for particle in self.particles]
            if all(type(x) == type(values[0]) for x in values):
                if all(x == values[0] for x in values):
                    return values[0]
            return values
        except (AttributeError,IndexError,TypeError) as e: ## FIXME redundant?
            raise AttributeError(e)

    def cmpProperties(self,other,
                      properties = ['Z2parity','spin','colordim','eCharge','mass','totalwidth']):
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

        #Check if any of its particles matches other
        if isinstance(other,Particle):
            otherParticles = [other]
        elif isinstance(other,MultiParticle):
            otherParticles = other.particles

        for otherParticle in otherParticles:
            if any(p.cmpProperties(otherParticle,properties) == 0 for p in self.particles):
                return 0

        cmpv = self.particles[0].cmpProperties(otherParticles[0],properties)
        return cmpv

    def __add__(self, other):
        """
        Define addition of two Particle objects
        or a Particle object and a MultiParticle object.
        The result is a MultiParticle object containing
        both particles.
        """

        if other is self or self.contains(other): #Check if other is self or a subset of self
            return self
        #Check if self is a subset of other
        if other.contains(self):
            return other
        elif isinstance(other,MultiParticle):
            addParticles = [ptc for ptc in other.particles if not self.contains(ptc)]
        elif isinstance(other,Particle):
            addParticles = [other]

        combinedParticles = self.particles + addParticles
        combined = MultiParticle(particles = combinedParticles)

        return combined

    def __radd__(self,other):
        return self.__add__(other)

    def __iadd__(self,other):

        addParticles = []
        if isinstance(other,MultiParticle):
            addParticles = [ptc for ptc in other.particles if not self.contains(ptc)]
        elif isinstance(other,Particle):
            if not self.contains(other):
                addParticles = [other]
        if addParticles:
            self.particles += addParticles[:]
            self.particles = sorted(self.particles)
            #Since the multiparticle changed, reset comparison tracking:
            self._comp = {self._id : 0}
            for ptc in self.particles:
                self._comp[ptc._id] = 0

        return self

    def getPdgs(self):
        """
        pdgs appearing in MultiParticle
        :return: list of pgds of particles in the MultiParticle
        """

        pdgs = [particle.pdg for particle in self.particles]

        return pdgs

    def getLabels(self):
        """
        labels appearing in MultiParticle
        :return: list of labels of particles in the MultiParticle
        """

        labels = [particle.label for particle in self.particles]

        return labels

    def isNeutral(self):
        """
        Return True if ALL particles in particle list are neutral.

        :return: True/False
        """

        neutral = all(particle.isNeutral() for particle in self.particles)
        return neutral

    def isMET(self):
        """
        Checks if all the particles in self can be considered as MET.

        :return: True/False
        """

        met = all(particle.isMET() for particle in self.particles)
        return met

    def contains(self,particle):
        """
        Check if MultiParticle contains the Particle object or MultiParticle object.

        :param particle: Particle or MultiParticle object

        :return: True/False
        """

        if isinstance(particle,MultiParticle):
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


class ParticleList(object):
    """
    Simple class to store a list of particles.
    """

    _instances = set()
    _lastID = 0


    def __new__(cls,particles):
        """
        Creates a particle list. If a list with the exact same particles have
        already been created return this list instead.
        Assigns an ID to the instance using the class ParticleList._instance
        list. Reset the comparison dictionary.

        :param particles: List of Particle or MultiParticle objects (list)
        """

        pList = sorted(particles)
        for obj in ParticleList.getinstances():
            if len(obj) != len(pList):
                continue
            if any(ptc is not obj.particles[iptc] for iptc,ptc in enumerate(pList)):
                continue
            return obj

        newList = super(ParticleList, cls).__new__(cls)
        newList.particles = pList[:]
        newList._id = ParticleList.getID()
        newList._comp = {newList._id : 0}
        ParticleList._instances.add(weakref.ref(newList))
        return newList

    def __getnewargs__(self):
        """
        Required for unpickling the object.
        When loading the pickled object, it will call __new__ with the
        arguments returned by this method.
        """
        return (self.particles,)

    def __getstate__(self):
        """
        Makes sure the particle list ID and comparison matrix
        are not stored in the pickle file, so they are dynamically assigned
        when the pickle file is loaded.
        """

        attrDict = dict(self.__dict__.items())
        #Make sure pickled objects do no store ID nor comparison dict
        attrDict.pop('_id',None)
        attrDict.pop('_comp',None)

        return attrDict

    def __hash__(self):
        return self._id

    @classmethod
    def getinstances(cls):
        dead = set()
        instances = []
        for ref in ParticleList._instances:
            obj = ref()
            if obj is not None:
                instances.append(obj)
            else:
                dead.add(ref)
        ParticleList._instances -= dead
        return instances

    @classmethod
    def getID(cls):
        if len(ParticleList.getinstances()) == 0:
            ParticleList._lastID = 0
        else:
            ParticleList._lastID += 1
        return ParticleList._lastID

    def __cmp__(self,other):
        """
        Compares two particle lists irrespective of the particle ordering.

        :param other:  particle list to be compared (ParticleList object)

        :return: -1 if self < other, 1 if self > other, 0 is self == other
        """

        #First check if we have already compared to this object
        if other._id in self._comp:
            return self._comp[other._id]

        if len(self) != len(other):
            comp = len(self) > len(other)
            if comp:
                comp = 1
            else:
                comp = -1
            self._comp[other._id] = comp
            other._comp[self._id] = comp
            return comp

        #Compare even final states irrespective of ordering:
        for particles in itertools.permutations(self.particles):
            particles = list(particles)
            if particles == other.particles:
                self._comp[other._id] = 0
                other._comp[self._id] = 0
                return 0

        comp = self.particles > other.particles
        if comp:
            comp = 1
        else:
            comp = -1
        self._comp[other._id] = comp
        other._comp[self._id] = -comp

        return comp

    def __lt__( self, p2 ):
        return self.__cmp__(p2) == -1

    def __gt__( self, p2 ):
        return self.__cmp__(p2) == 1

    def __eq__( self, p2 ):
        return self.__cmp__(p2) == 0

    def __ne__( self, p2 ):
        return self.__cmp__(p2) != 0

    def __iter__(self):
        return iter(self.particles)

    def __getitem__(self, i):
        return self.particles[i]

    def __setitem__(self, i, value):
        self.particles[i] = value

    def __len__(self):
        return len(self.particles)

    def __str__(self):
        return str(self.particles)

    def __repr__(self):
        return self.__str__()
