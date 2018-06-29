"""
.. module:: particle
   :synopsis: Defines the particle class and particle list class, their methods and related functions

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

import copy

class Particle(object):
    """
    An instance of this class represents a single particle. 
    The properties are: label, pdg, mass, electric charge, color charge, width 
    """

    def __init__(self, **kwargs):
        """ 
        
        Initializes the particle.
        Possible properties: 
        Z2parity: str, 'even' or 'odd'
        label: str, e.g. 'e-'
        pdg: number in pdg
        mass: mass of the particle
        echarge: electric charge as multiples of the unit charge
        colordim: color dimension of the particle 
        spin: spin of the particle
        width: total width
        decays: possible decays in pyslha format e.g. [ 0.5[1000022, -1], 0.5[1000022, -2] ]
                 
        """  

        for key,value in kwargs.items():
            self.addProperty(key,value)

    def __cmp__(self,other):
        """
        Compares particle with other.
        The comparison is done based only on the label.
        :param other:  particle to be compared (Particle object)
        
        :return: -1 if self < other, 0 if self == other, +1, if self > other.
        """    
        
        if not isinstance(other,(ParticleList,Particle,ParticleWildcard)):
            return +1
        
        return self.cmpProperties(other, properties=['Z2parity','label']) 
                 

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
    
    def describe(self):
        return str(self.__dict__)
    
    
    def addProperty(self,label,value,overwrite=False):
        """
        Add property with label and value.
        If overwrite = False and property already exists, it will not be overwritten
        :parameter label: property label (e.g. "mass", "label",...)
        :parameter value: property value (e.g. 171*GeV, "t+",...)
        """
        
        if overwrite:
            setattr(self, label, value)
        elif not hasattr(self, label) or getattr(self, label) is None:
            setattr(self, label, value)    

    def eqProperties(self,other, properties = ['spin','colordim','eCharge']):
        """
        Check if particle has the same properties (default is spin, colordim and eCharge)
        as other. Only compares the attributes which have been defined in both objects.
        
        :param other: a Particle, ParticleList or ParticleWildCard object
        :param properties: list with properties to be compared. Default is spin, colordim and eCharge
        
        :return: True if all properties are the same, False otherwise.
        """
        
        if self.cmpProperties(other, properties=properties) == 0:
            return True
        else:
            return False
            
    def cmpProperties(self,other, properties = ['spin','colordim','eCharge']):
        """
        Compare properties (default is spin, colordim and eCharge).
        Return 0 if properties are equal, -1 if self < other and 1 if self > other.
        Only compares the attributes which have been defined in both objects.
        The comparison is made in hierarchical order, following the order
        defined by the properties list.
        
        :param other: a Particle, ParticleList or ParticleWildCard object
        :param properties: list with properties to be compared. Default is spin, colordim and eCharge
        
        :return: 0 if properties are equal, -1 if self < other and 1 if self > other.
        """

        if isinstance(other,(ParticleList,ParticleWildcard)):
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
            
        return pConjugate

        
    def isNeutral(self):
        """
        Return True if the particle label is MET or if it color and charge
        neutral.
        
        :return: True/False
        """
        
        if self.eCharge == 0 and self.colordim == 0:
            return True
        
        return False
    
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
    
        

class ParticleList(object):

    """ An instance of this class represents a list of particle object to allow for inclusive expresions such as jet. 
        The properties are: label, pdg, mass, electric charge, color charge, width 
    """
    
    def __init__(self, label, particles, **kwargs):

        """ 
        Initializes the particle list.
        """        
                
        self.label = label
        self.particles = particles
        
        for key,value in kwargs.items():
            self.addProperty(key,value)
        
        
     
    def __cmp__(self,other):
        """
        Compares the list with another list or a single particle.
        :param other: particle list (ParticleList) or single particle (Particle) to be compared
        :return: If other is a single particle, returns:
                    -1 if all particles in self < other, 
                     0 if any particle in self == other,
                    +1, if any particle in self > other.
                If other is a particle list, returns 
                -1 if self < other, 0 if self == other, +1, if self > other.     
        """        
        
        if isinstance(other,ParticleList):
            if self.particles != other.particles:
                comp = self.particles > other.particles
                if comp:
                    return 1
                else:
                    return -1 
            else:
                return 0
        
        if isinstance(other,Particle):
            if other in self.particles:
                return 0
            else:
                for p in self.particles:
                    if p > other: return +1
                return -1
            
    def __eq__( self, other ):
        return self.__cmp__(other) == 0
        
    def __ne__( self, other ):
        return self.__cmp__(other) != 0
    
    def __lt__( self, other ):
        return self.__cmp__(other) == -1

    def __gt__( self, other ):
        return self.__cmp__(other) == 1

    def __str__(self):         
        return self.label
    
    def __repr__(self):
        return self.__str__()    
    
    def __getattribute__(self,name):
        """
        If ParticleList does not have attribute, return a list
        if the attributes of each particle in self.particles.
        If not all particles have the attribute, it will raise an exception.
        If all particles share a common attribute, a single value
        will be returned.
        
        :parameter name: Attribute string
        
        :return: Attribute or list with the attribute values in self.particles
        """
        
        try:
            return object.__getattribute__(self, name)
        except:
            values = [getattr(particle,name) for particle in self.particles]
            if all(type(x) == type(values[0]) for x in values):
                if all(x == values[0] for x in values):
                    return values[0]
            return values
        
    def addProperty(self,label,value,overwrite=False):
        """
        Add property with label and value.
        If overwrite = False and property already exists, it will not be overwritten
        :parameter label: property label (e.g. "mass", "label",...)
        :parameter value: property value (e.g. 171*GeV, "t+",...)
        """
        
        if overwrite:
            setattr(self, label, value)
        elif not hasattr(self, label) or getattr(self, label) is None:
            setattr(self, label, value)    

    def eqProperties(self,other, properties = ['spin','colordim','eCharge']):
        """
        Check if particle has the same properties (default is spin, colordim and eCharge)
        as other. Only compares the attributes which have been defined in both objects.
        
        :param other: a Particle, ParticleList or ParticleWildCard object
        :param properties: list with properties to be compared. Default is spin, colordim and eCharge
        
        :return: True if all properties are the same, False otherwise.
        """
        
        if self.cmpProperties(other, properties=properties) == 0:
            return True
        else:
            return False

        
    def cmpProperties(self,other, properties = ['spin','colordim','eCharge']):
        """
        Compares the properties in self with the ones in other (default is spin, colordim and eCharge).
        Return 0 if properties are equal, -1 if self < other and 1 if self > other.
        Only compares the attributes which have been defined in both objects.
        The comparison is made in hierarchical order, following the order
        defined by the properties list.
        
        If the property in the ParticleList is a list, check if any value of the list
        matches the property of other. If both properties are lists, check if the
        lists share at least one common element.
        
        :param other: a Particle, ParticleList or ParticleWildCard object
        :param properties: list with properties to be compared. Default is spin, colordim and eCharge
        
        :return: 0 if properties are equal, -1 if self < other and 1 if self > other.
        """
        
        if isinstance(other,ParticleWildcard):
            return -1*other.cmpProperties(self,properties=properties)
        
        for prop in properties:            
            if not hasattr(self,prop) or not hasattr(other,prop):
                continue
            
            x = getattr(self,prop)
            y = getattr(other,prop)
            
            if x == y:
                continue
            
            if not isinstance(x,list) and not isinstance(y,list): #Compare single values:
                if x > y:
                    return 1
                else:
                    return -1            
            elif isinstance(x,list) and not isinstance(y,list): #Compare list to single value
                if any(xval == y for xval in x):
                    continue
                elif any(xval > y for xval in x):
                    return 1
                else:
                    return -1
            
            elif not isinstance(x,list) and isinstance(y,list): #Compare list to single value
                if any(yval == x for yval in y):
                    continue
                elif any(yval > x for yval in y):
                    return -1
                else:
                    return 1
            else:
                if any(xval in y for xval in x): #See if any element in the two lists match
                    continue
                elif x > y:
                    return 1
                else:
                    return -1
            
        return 0




    def describe(self):
        return str(self.__dict__)  
        

    def getPdgs(self):
        """
        pdgs appearing in ParticleList
        :return: list of pgds of particles in the ParticleList
        """
        
        pdgs = [particle.pdg for particle in self.particles]
        
        return pdgs
        
    def getLabels(self):
        """
        labels appearing in ParticleList
        :return: list of labels of particles in the ParticleList
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



class ParticleWildcard(Particle):
    """
    A particle wildcard class. It will return True when compared to any other Particle or ParticleList object.
    The only exception is if the ParticleWildcard has a defined parity. In this case it will only be equal
    to particles with the same parity.
    """
    
    def __init__(self,**kwargs):
        if not kwargs or not 'label' in kwargs:
            kwargs['label'] = '*'
        Particle.__init__(self,**kwargs)
        
    def __repr__(self):
        return self.__str__()



    def cmpProperties(self,other, properties = ['Z2parity']):
        """
        Returns 0 if the Z2parities match, otherwise
        -1 if self.Z2parity < other.Z2parity and 1 if self.Z2parity > other.Z2parity.
        
        :param other: a Particle, ParticleList or ParticleWildCard object
        :param properties: (dummy) list with properties to be compared. It is not used.
        
        :return: 0 if properties are equal, -1 if self < other and 1 if self > other.
        """
        
        prop = 'Z2parity'
                
        if not hasattr(self,prop) or not hasattr(other,prop):
            return 0
        else:
            x = getattr(self,prop)
            y = getattr(other,prop)
            if isinstance(y,list):
                if x in y:
                    return 0
                elif any(yval > x for yval in y):
                    return -1
                else:
                    return 1
            else:
                if x == y:
                    return 0
                elif x > y:
                    return 1
                else:
                    return -1
