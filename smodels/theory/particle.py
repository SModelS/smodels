"""
.. module:: particle
   :synopsis: Defines the particle class and particle list class, their methods and related functions

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

from smodels.tools.physicsUnits import GeV


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
        The comparison is done based on the label.
        :param other:  particle to be compared (Particle object)
        :return: -1 if self < other, 0 if self == other, +1, if self > other.

        Comparison to '*' always returns True. (= Particle Wild Card)
        """    
        
        if not isinstance(other,(ParticleList,Particle,ParticleWildcard)):
            return +1
        
        if isinstance(other,ParticleList):
            return other.__cmp__(self)
        
        if isinstance(other,ParticleWildcard):
            return other.__cmp__(self)
        
        #Make sure particles with distinct Z2 parity are distinct
        if hasattr(self, 'Z2parity') and hasattr(other, 'Z2parity'):
            if self.Z2parity != other.Z2parity:
                comp = self.Z2parity > other.Z2parity
                if comp:
                    return 1
                else:
                    return -1
                
        if hasattr(self, 'label') and hasattr(other, 'label'):
            if self.label != other.label:            
                comp = self.label > other.label
                if comp:
                    return 1
                else:
                    return -1       
            
        return 0  #Particles are equal
                 

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


    def copy(self):
        """
        Make a copy of self.
        
        :return: A Particle object identical to self
        """

        p = Particle()
        for name in dir(self):
            if '__' in name:
                continue
            setattr(p,name,getattr(self,name))

        return p

    def chargeConjugate(self):
        """
        Returns the charge conjugate particle (flips the sign of eCharge).        
        If it has a pdg property also flips its sign.
        The charge conjugate name is defined as the original name plus "~" or
        if the original name ends in "+" ("-"), it is replaced by "-" ("+")

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
                
        return pConjugate

        
    def isNeutral(self):
        """
        Return True if the particle label is MET or if it color and charge
        neutral.
        
        :return: True/False
        """
        
        if self.label == 'MET':
            return True
        if self.eCharge == 0 and self.colordim == 0:
            return True
        
        return False
        
        


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
            
    def __eq__( self, p2 ):
        return self.__cmp__(p2) == 0
        
    def __ne__( self, p2 ):
        return self.__cmp__(p2) != 0
    
    def __lt__( self, p2 ):
        return self.__cmp__(p2) == -1

    def __gt__( self, p2 ):
        return self.__cmp__(p2) == 1

    def __str__(self):         
        return self.label
    
    def __repr__(self):
        return self.__str__()    
    
    def __getattribute__(self,name):
        """
        If ParticleList does not have attribute, return a list
        if the attributes of each particle in self.particles.
        If all particles share a common attribute, a single value
        will be return.
        
        :parameter name: Attribute string
        
        :return: Attribute or list with the attribute values in self.particles
        """
        
        try:
            return object.__getattribute__(self, name)
        except:
            values = [getattr(particle,name) for particle in self.particles]
            if len(list(set(values))) == 1:
                values = values[0]
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
        
        neutrals = [particle.isNeutral() for particle in self.particles]
        return sum(neutrals) == len(self.particles)


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

    def __cmp__(self,other):
        #Make sure particles with distinct Z2 parity are distinct
        if hasattr(self, 'Z2parity') and hasattr(other, 'Z2parity') and self.Z2parity != other.Z2parity:
            comp = self.Z2parity > other.Z2parity
            if comp:
                return 1
            else:
                return -1        
        elif isinstance(other,(Particle,ParticleList,ParticleWildcard)):
            return 0
        else:
            return -1

    def __eq__(self,other):
        return self.__cmp__(other) == 0  
    
    def __ne__(self,other):
        return self.__cmp__(other) != 0    
        