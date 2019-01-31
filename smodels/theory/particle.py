"""
.. module:: particle
   :synopsis: Defines the particle class and particle list class, their methods and related functions

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

import copy,itertools,weakref

class Particle(object):
    """
    An instance of this class represents a single particle. 
    The properties are: label, pdg, mass, electric charge, color charge, width 
    """

    def __init__(self, **kwargs):
        """ 
        
        Initializes the particle.
        Possible properties: 
        Z2parity: int, +1 or -1
        label: str, e.g. 'e-'
        pdg: number in pdg
        mass: mass of the particle
        echarge: electric charge as multiples of the unit charge
        colordim: color dimension of the particle 
        spin: spin of the particle
        width: total width
        decays: possible decays in pyslha format e.g. [ 0.5[1000022, -1], 0.5[1000022, -2] ]
                 
        """  

        self._static = False
        self.id = id(self)
        self._comp = {self.id : 0}
        for attr,value in kwargs.items():
            if not attr == '_static':
                setattr(self,attr,value)
                
        #Leave the static attribute for last:
        if '_static' in kwargs:
            self._static = kwargs['_static']

    def __cmp__(self,other):
        """
        Compares particle with other.
        The comparison is based on the particle properties.
        If the particles differ they are sorted according to their label.
        If in addition the particles differ, but have the same label, they
        are sorted according to their properties.
        
        :param other:  particle to be compared (Particle object)
        
        :return: -1 if self.label < other.label, 0 if self == other, +1, if self.label > other.label.
        """    
        
        if not isinstance(other,(MultiParticle,Particle)):
            return +1

        #First check if we have already compared to this object
        idOther = other.id        
        if idOther in self._comp:
            return self._comp[idOther]

        idSelf = self.id
        if idSelf in other._comp:
            return -other._comp[idSelf]

        cmpProp = self.cmpProperties(other) #Objects have not been compared yet.
        if not self._static:
            self._comp[idOther] = cmpProp
        if not other._static:
            other._comp[idSelf] = -cmpProp
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
    
    def __setattr__(self,attr,value):
        """
        Override setattr method.
        If the _static attribute is True, will not
        change the particle attribute.
        """
          
        if attr == '_static':
            self.__dict__[attr] = value
        elif self._static is False:            
            self.__dict__[attr] = value
     
    def __setstate__(self, state):
        """
        Override setstate method. Required for pickling.
        """

        self._static = False
        self.__dict__.update(state)

    def __add__(self, other):
        """
        Define addition of two Particle objects
        or a Particle object and a MultiParticle object.
        The result is a MultiParticle object containing
        both particles.
        """

        if not isinstance(other,(MultiParticle,Particle)):
            raise TypeError("Can only add particle objects")
        elif isinstance(other,MultiParticle):
            return other.__add__(self)
        elif self.contains(other):
            return self
        elif other.contains(self):
            return other
        else:
            combined = MultiParticle(label = 'multiple', particles= [self,other])
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
        Make a copy of self (using deepcopy)
        
        :return: A Particle object identical to self
        """

        return copy.deepcopy(self)

    def chargeConjugate(self,label=None):
        """
        Returns the charge conjugate particle (flips the sign of eCharge).        
        If it has a pdg property also flips its sign.
        If label is None, the charge conjugate name is defined as the original name plus "~" or
        if the original name ends in "+" ("-"), it is replaced by "-" ("+")

        :parameter label: If defined, defines the label of the charge conjugated particle.

        :return: the charge conjugate particle (Particle object)
        """
        
        pConjugate = self.copy()
        pConjugate.id = id(pConjugate)
        pConjugate._static = False #Temporarily set it to False to change attributes
        pConjugate._comp = {pConjugate.id : 0}
                    
        if hasattr(pConjugate, 'pdg') and pConjugate.pdg:
            pConjugate.pdg *= -1       
        if hasattr(pConjugate, 'eCharge') and pConjugate.eCharge:
            pConjugate.eCharge *= -1    
        if hasattr(pConjugate, 'label'):                
            if pConjugate.label[-1] == "+":
                pConjugate.label = pConjugate.label[:-1] + "-"
            elif pConjugate.label[-1] == "-":
                pConjugate.label = pConjugate.label[:-1] + "+"
            elif pConjugate.label[-1] == "~":
                pConjugate.label = pConjugate.label[:-1]
            else:
                pConjugate.label += "~"            
        
        if not label is None:
            pConjugate.label = label
            
        pConjugate._static = self._static #Restore the initial state

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
        If the isMET attribute has not been defined, it will return True/False is isNeutral() = True/False.
        Else it will return the isMET attribute.
        
        :return: True/False
        """
        
        if hasattr(self,'_isMET'):
            return self._isMET
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

    """ An instance of this class represents a list of particle object to allow for inclusive expresions such as jet. 
        The properties are: label, pdg, mass, electric charge, color charge, width 
    """
    
    def __init__(self, label, particles, **kwargs):

        """ 
        Initializes the particle list.
        """        

        self._static = False
        self.label = label
        self.particles = particles
        Particle.__init__(self,**kwargs)
        self._comp = dict([[self.id,0]] + [[ptc.id,0] for ptc in particles])

    def __getattribute__(self,attr):
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
            return super(MultiParticle,self).__getattribute__(attr) #Python2
        except:
            pass

        try:
            return super().__getattribute__(attr) #Python3
        except:
            pass
         
        try:
            values = [getattr(particle,attr) for particle in self.particles]
            if all(type(x) == type(values[0]) for x in values):
                if all(x == values[0] for x in values):
                    return values[0]
            return values
        except:
            raise AttributeError
            
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

    def __getstate__(self):
        """
        Override getstate method. Required for pickling.
        """  
        return self.__dict__

    def __setstate__(self, state):
        """
        Override setstate method. Required for pickling.
        """

        self._static = False
        self.__dict__.update(state)

    def __add__(self, other):
        """
        Define addition of two Particle objects
        or a Particle object and a MultiParticle object.
        The result is a MultiParticle object containing
        both particles.
        """

        if not isinstance(other,(MultiParticle,Particle)):
            raise TypeError("Can not add a Particle object to %s" %type(other))
        elif other is self or self.contains(other): #Check if other is self or a subset of self
            return self
        #Check if self is a subset of other
        if other.contains(self):
            return other        
        elif isinstance(other,MultiParticle):
            addParticles = other.particles
        elif isinstance(other,Particle):
            addParticles = [other]

        combinedParticles = [ptc for ptc in addParticles if not self.contains(ptc)]
        combinedParticles += self.particles[:]

        combined = MultiParticle(label = 'multiple', particles = combinedParticles)

        return combined
    
    def __radd__(self,other):
        return self.__add__(other)
    
    def __iadd__(self,other):
        
        if isinstance(other,MultiParticle):
            self.particles += [ptc for ptc in other.particles if not self.contains(ptc)]
        elif isinstance(other,Particle):
            if not self.contains(other):
                self.particles.append(other)
                if not self._static:
                    self._comp[other.id] = 0

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

        if not isinstance(particle,(Particle,MultiParticle)):
            raise False
        elif isinstance(particle,MultiParticle):
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

    def __init__(self,particles):

        self.particles = particles[:]
        self.id = id(self)
        self._comp = {self.id : 0}
        self._static = False
        self._instances.add(weakref.ref(self))

    def __hash__(self):
        return self.id

    @classmethod
    def getVertex(cls,particles):
        pList = sorted(particles)
        for obj in ParticleList.getinstances():
            if len(obj) != len(pList):
                continue
            if all(ptc is obj.particles[iptc] for iptc,ptc in enumerate(pList)):
                return obj
        return ParticleList(pList)

    @classmethod
    def getNinstances(cls):
        dead = set([])
        for ref in cls._instances:
            obj = ref()
            if obj is None:
                dead.add(ref)
        cls._instances -= dead
        return len(cls._instances)

    @classmethod
    def getinstances(cls):
        dead = set()
        for ref in cls._instances:
            obj = ref()
            if obj is not None:
                yield obj
            else:
                dead.add(ref)
        cls._instances -= dead


    def identical(self,other):        
        if len(self) != len(other):
            return False
        if any((ptc is not other.particles[iptc]) for iptc,ptc in enumerate(self.particles)):
            return False
        
        return True
        
    def __cmp__(self,other):
        """
        Compares two particle lists irrespective of the particle ordering.
        
        :param other:  particle list to be compared (ParticleList object)
        
        :return: -1 if self < other, 1 if self > other, 0 is self == other 
        """    
        
        if not isinstance(other,ParticleList):
            raise ValueError

        #First check if we have already compared to this object
        idOther = other.id  
        if idOther in self._comp:
            return self._comp[idOther]

        idSelf = self.id
        if idSelf in other._comp:
            return -other._comp[idSelf]

        #Compare even final states irrespective of ordering:
        for particles in itertools.permutations(self.particles):
            particles = list(particles)
            if particles == other.particles:
                self._comp[other.id] = 0
                other._comp[self.id] = 0
                return 0
        
        comp = self.particles > other.particles
        if comp:
            comp = 1
        else:
            comp = -1
        if not self._static:
            self._comp[idOther] = comp
        if not other._static:
            other._comp[idSelf] = -comp
            
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
        return str([str(ptc) for ptc in self.particles]) 
        
    def __repr__(self):
        return self.__str__()        

    def resetComp(self):
        self._comp = {self.id : 0}        