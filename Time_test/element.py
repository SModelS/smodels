#!/usr/bin/env python

"""
.. module:: Theory.Element
   :synopsis: missing
    
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
    
"""
from ParticleNames import simParticles
from branch import Branch
import copy,sys

class Element(object):
    """
    Abstract class. DO NOT INSTANTIATE!
    """


    def __init__(self, info=None):
        """
        Constructor
        """
        self.branches = [Branch(),Branch()]
        self.weight = None
                
        if info:
            if type(info) == type('string'):   #Creates element from particle string
                st = info.replace(" ","").replace("'","")
                st = st[st.find("[[["):st.find("]]]")+3]
                b1=st[2:st.find("]],[[")+1]
                b2=st[st.find("]],[[")+4:st.find("]]]")+1]
                self.branches[0] = copy.deepcopy(Branch(b1))
                self.branches[1] = copy.deepcopy(Branch(b2))
            elif type(info) == type([]) and type(info[0]) == type(Branch()): #Creates element from branch pair
                for ib,branch in enumerate(info): self.branches[ib] = copy.deepcopy(branch) 


    def __eq__ (self, other):
        return self.isEqual(other)
    
    def __str__ ( self ):
        """ returns the canonical name of the element, e.g. [[jet],[jet]] """
        particleString = str(self.getParticles()).replace(" ","").replace("'","")
        return particleString

    def copy(self):
      elementcopy = Element()
      elementcopy.branches = self.branches
      elementcopy.weight = self.weight

      return elementcopy

    def isEqual(self,other,useDict=True):
        """ Compare two Elements
        If all masses and particles are equal, returns True,
        otherwise returns False
        If order = False, test both branch orderings (for an element doublet only) 
        If useDict=True, allow for inclusive particle labels"""
                        
        if type(self) != type(other): return False
        mass = self.getMasses()
        massA = other.getMasses()
                
        if self.particlesMatch(other,order=True,useDict=useDict) and mass == massA: return True               
        other_b = other.switchBranches()            
        mass_b = other_b.getMasses()
        if self.particlesMatch(other_b,order=True,useDict=useDict) and mass_b == massA: return True

        return False


    def isEqual2(self,other,useDict=True):
        """ Compare two Elements
        If all masses and particles are equal, returns True,
        otherwise returns False
        If order = False, test both branch orderings (for an element doublet only) 
        If useDict=True, allow for inclusive particle labels"""
                        
        if type(self) != type(other): return False
        mass = self.getMasses()
        massA = other.getMasses()
        ptcs = self.getParticles()
        ptcsA = other.getParticles()                
        if simParticles(ptcsA,ptcs,useDict) and mass == massA: return True               
        ptcs_b = [ptcsA[1],ptcsA[0]]
        mass_b = [massA[1],massA[0]]
        if simParticles(ptcs_b,ptcs,useDict) and mass_b == mass: return True

        return False




    def particlesMatch(self,other,order=False,useDict=True):
        """ Compare two Elements
        If particles match, returns True,
        otherwise returns False
        If order = False, test both branch orderings (for an element doublet only) 
        If useDict=True, allow for inclusive particle labels"""
        
        if type(self) != type(other): return False
        ptcs = self.getParticles()
        ptcsA = other.getParticles()
        if simParticles(ptcsA,ptcs,useDict): return True
        if not order:
            ptcs_b = other.switchBranches().getParticles()   
            if simParticles(ptcsA,ptcs_b,useDict): return True
        
        return False    

    def switchBranches(self):
        """ If the element contains a pair of branches, switches them"""
        newEl = self.copy()
        if len(self.branches) == 2: newEl.branches = [newEl.branches[1],newEl.branches[0]]
        return newEl

    def getParticles(self):
        """
        Returns the array of particles in the element
        """
        ptcarray = []
        for branch in self.branches: ptcarray.append(branch.particles)
        return ptcarray
    
    def getMasses(self):
        """
        Returns the array of masses in the element
        """
        massarray = []
        for branch in self.branches: massarray.append(branch.masses)
        return massarray
    
    def getDaughters(self):
        """
        Return the pair of daughter IDs (can be None, if the element does not have a definite daughter)    
        """        
        return (self.branches[0].daughterID,self.branches[1].daughterID)
    
    def getMothers(self):
        """
        Return the pair of mother IDs (can be None, if the element does not have a mother daughter)    
        """        
        return (self.branches[0].momID,self.branches[1].momID)
    
    def getEinfo(self):
        """
        Get global topology info from particle string.
        """
        vertnumb = []
        vertparts = []
        for branch in self.branches:
            vertnumb.append(len(branch.masses))
            vertparts.append([len(ptcs) for ptcs in branch.particles])
            if len(vertparts[len(vertparts)-1]) == vertnumb[len(vertnumb)-1]-1:
                vertparts[len(vertparts)-1].append(0)  #Append 0 for stable LSP
        return {"vertnumb" : vertnumb, "vertparts" : vertparts}
    

