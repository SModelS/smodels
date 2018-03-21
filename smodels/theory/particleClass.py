"""
.. module:: particles
   :synopsis: Defines the particle class and particle list class, their methods and related functions

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

from smodels.tools.smodelsLogging import logger
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from copy import copy


class Particles(object):
	
    """ An instance of this class represents a single particle. 
        The properties are: label, pdg, mass, electric charge, color charge, width 
    """

    def __init__(self, Z2parity, label, pdg, mass, eCharge, colordim, spin, width, branches):
	    
        """ 
        Initializes the particle.
        All properties must be defined, if not known set to None.

        :parameter label: string, e.g. 'e-'
        :parameter pdg: number in pdg
        :parameter mass: in MeV
        :parameter charge: as multiples of the unit charge
        :parameter color charge: 
        :parameter width: total width 

        """

        self.Z2parity = Z2parity
        self.label = label
        self.pdg = pdg
        self.mass = mass
        self.eCharge = eCharge
        self.colordim = colordim
        self.spin = spin
        self.width = width
        self.branches = branches
        
    def __cmp__(self,other):
        """
        Compares the branch with other.        
        The comparison is made based on .
        OBS: The particles inside each vertex MUST BE sorted (see branch.sortParticles())         
        :param other:  branch to be compared (Branch object)
        :return: -1 if self < other, 0 if self == other, +1, if self > other.
        """
        
        if self.vertnumb != other.vertnumb:
            comp = self.vertnumb > other.vertnumb
            if comp: return 1
            else: return -1
        elif self.vertparts != other.vertparts:
            comp = self.vertparts > other.vertparts
            if comp: return 1
            else: return -1
        elif self.particles != other.particles:    
            comp = self.particles > other.particles
            if comp: return 1
            else: return -1   
        else:
            m1m2eq, compm = compareBSMparticles( self.BSMparticles, other.BSMparticles )                 
            if not m1m2eq:              
                if compm: return 1
                else: return -1
                
            else:
                return 0  #Branches are equal 

    def __cmp__(self,other): 
        """
        Compares particle with other.
        The comparison is made based on the label.
        :param other:  particle to be compared (Particle object)
        :return: -1 if self < other, 0 if self == other, +1, if self > other.
                
        Allows for u, d, s quarks to be considered equal. Allows for all neutrinos to be considered equal.
        """    
        if self.label != other.label:
        
            if ( self.label == 'qd' or self.label == 'qu' or self.label == 'qs' ) and ( other.label == 'qd' or other.label == 'qu' or other.label == 'qs' ):
                return 0      
                
            if ( self.label == 'nue' or self.label == 'numu' or self.label == 'nuta' ) and ( other.label == 'nue' or other.label == 'numu' or other.label == 'nuta' ):    
                        return 0                  
        
            comp = self.label > other.label
            if comp: return 1
            else: return -1
   
        else:
            return 0  #Particles are equal         

    def __lt__( self, p2 ):
        return self.__cmp__ ( p2 ) == -1

    def __eq__( self, p2 ):
        return self.__cmp__ ( p2 ) == 0
        


    def __str__(self): 
        """
        Returns all properties of the particle
        """
        return str(self.__dict__ )




    def chargeConjugate(self):
        """
        Returns the charge conjugate particle (flips the sign of eCharge).
        If it has a _pid property also flips its sign.
        The charge conjugate name is defined as the original name plus "*" or
        if the original name ends in "+" ("-"), it is replaced by "-" ("+")

        :return: the charge conjugate particle (Particle object)
        """

        pConjugate = copy(self)
        
        pConjugate.eCharge *= -1
        pConjugate.pdg *= -1

        if pConjugate.label[-1] == "+":
            pConjugate.label = pConjugate.label[:-1] + "-"
        elif pConjugate.label[-1] == "-":
            pConjugate.label = pConjugate.label[:-1] + "+"

                
        return pConjugate




class ParticleList(object):

    """ An instance of this class represents a list of particle object to allow for inclusive expresions such as jet. 
        The properties are: label, pdg, mass, electric charge, color charge, width 
    """
    
    def __init__(self, label, particles):

        """ 
        Initializes the particle list.

        :parameter label: string, e.g. 'jet'
        :parameter particles: particle objects to be added to list

        """        
        self.label = label
        self.particles = particles
        
    
    def __eq__(self, other): 
        return self.label == other.label    
    
    
    def __str__(self): 
        """
        :return: label and particles of the ParticleList
        """        
        return str(self.__dict__ )    
        
        
    def getPdgs(self):
        """
        pdgs appearing in list
        :return: list of pgds of particles in the ParticleList
        """        
        pdgs = [particle.pdg for particle in self.particles]
        return pdgs
        
        
    def getNames(self):
        """
        names appearing in list
        :return: list of names of particles in the ParticleList
        """
        names = [particle.label for particle in self.particles]
        return names
        
        
        




        
        
def particleInList(particle, particleList, inlist=True):
    """ 
    checks whether particle is in particleList
    :param particle: particle either as Particle class object or as str
    :param particleList: list of particles, can contain Particle objects and ParticleList objects, type must be list
    :param inlist: if set to False only checks labels of objects in list and ignores Particles contained in ParticleLists
    :return: True/False (particle is in the list/ not in the list)
    """
    
    found = False
    
    for p in particleList:        
        if particle == p: 
            found = True
    
    if not inlist: return found
            
    if hasattr(p,'particles') and found==False:
        for plist in particleList:
            for ptc in p.particles:
                if particle == ptc: 
                    found = True

    if found: return True
    else: return False 



def sortParticleList(ptcList):
    """
    sorts a list of particle or particle list objects by their label
    :param ptcList: particle list to be sorted
    :return: sorted list of particles
    """
    
    newPtcList = sorted(ptcList, key=lambda x: x.label) 
    return newPtcList
    
         


