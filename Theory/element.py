#!/usr/bin/env python

"""
.. module:: Theory.Element
   :synopsis: missing
    
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
    
"""
from ParticleNames import PtcDic, Reven, simParticles
from branch import Branch
import copy
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class Element(object):
    """
    Element class. Holds a pair of branches and the element weight (cross-section*BR)
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
    
    def __ne__ (self, other):
        return not self.isEqual(other)    
    
    def __str__ ( self ):
        """ returns the canonical name of the element, e.g. [[jet],[jet]] """
        particleString = str(self.getParticles()).replace(" ","").replace("'","")
        return particleString


    def isEqual(self,other,order=False,useDict=True):
        """ Compare two Elements
        If all masses and particles are equal, returns True,
        otherwise returns False
        If order = False, test both branch orderings (for an element doublet only) 
        If useDict=True, allow for inclusive particle labels"""
                           
        if type(self) != type(other): return False
        mass = self.getMasses()
        massA = other.getMasses()
                   
        if self.particlesMatch(other,order=True,useDict=useDict) and mass == massA: return True               
        if not order:
            other_b = other.switchBranches()            
            mass_b = other_b.getMasses()
            if self.particlesMatch(other_b,order=True,useDict=useDict) and mass == mass_b: return True
   
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
    
    def copy(self):
        """
        Creates a copy of itself (faster than deepcopy)
        """
        elcopy = Element()
        elcopy.branches = self.branches
        elcopy.weight = self.weight
        
        return elcopy
    
    def setMasses(self,mass,same_order=True,oppos_order=False):
        """
        Sets the element masses to the input mass array.
        If same_order, set the masses to the same branch ordering.
        If oppos_order, set the masses to the opposite branch ordering.
        If both same_order and oppos_order, set the masses to the smaller of the two orderings
        """
                
        if same_order and oppos_order:
            newmass = smallerMass(mass)
        elif same_order:
            newmass = mass
        elif oppos_order:
            newmass = [mass[1],mass[0]]
        else:
            logger.error('[Element.setMasses]: called with no possible ordering!')
            return False
        
        self.branches[0].masses = copy.deepcopy(newmass[0])
        self.branches[1].masses = copy.deepcopy(newmass[1])

    def switchBranches(self):
        """ If the element contains a pair of branches, switches them"""
#         newEl = copy.deepcopy(self)
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
    
    def getLength(self):
        """ Returns the maximum of the two branch lengths """
        
        return max(self.branches[0].getLength(),self.branches[1].getLength())
    
    def hasTopInList(self,elementList):
        """
        Checks if the element topology matches any of the topologies in the element list
        """
        if type(elementList) != type([]) or len(elementList) == 0: return False
        for element in elementList:
            if type(element) != type(self): continue
            info1 = self.getEinfo()
            info2 = element.getEinfo()
            info2_b = element.switchBranches().getEinfo()
            if info1 == info2 or info1 == info2_b: return True
        
        return False
        
      
    def isInList(self, listOfElements, igmass=False, useDict=True):
        """ checks if the element is present in the element list. If igmass=False also check if
        the analysis has the element mass array"""

        for el in listOfElements:
            if igmass:
                if self.particlesMatch(el,useDict): return True
            else:
                if self.isEqual(el,useDict): return True

        return False
    
    def checkConsistency(self):
        """
        Checks if the particles defined in the element exist and are consistent with the element info
        """

        info = self.getEinfo()
        for ib,branch in enumerate(self.branches):            
            for iv,vertex in enumerate(branch.particles):
                if len(vertex) != info['vertparts'][ib][iv]:
                    logger.error("[Element.checkConsistency]: Wrong syntax")
                    return False
                for ptc in vertex:
                    if not ptc in Reven.values() and not PtcDic.has_key(ptc):
                        logger.error("[Element.checkConsistency]: unknown particle"+ptc)
                        return False
        return True    
    
    
    def compressElement(self,DoCompress,DoInvisible,minmassgap):
        """ Keep compressing they original element and the derived ones till they can be compressed no more.
            Returns a list with the compressed elements."""


        added = True
        newElements = [copy.deepcopy(self)]        
#Keep compressing the new topologies generated so far until no new compressions can happen:
        while added:
            added = False
    #Check for mass compressed topologies   
            if DoCompress:
                for element in newElements:             
                    newel = element.massCompress(minmassgap)
                    if newel and not newel.hasTopInList(newElements): #Avoids double counting (conservative)
                        newElements.append(newel) 
                        added = True
      
    #Check for invisible compressed topologies (look for effective LSP, such as LSP + neutrino = LSP')      
            if DoInvisible:
                for element in newElements:
                    newel = element.invisibleCompress()
                    if newel and not newel.hasTopInList(newElements): #Avoids double counting (conservative)
                        newElements.append(newel) 
                        added = True

                         
        newElements.pop(0)  #Remove original element                            
        return newElements
    
         
    def massCompress(self,mingap):
        """ if two masses in this topology are degenerate, returns a compressed copy of the element. If no compression is
        possible, return None. """
        
        newelement = copy.deepcopy(self)
        vertnumb = self.getEinfo()["vertnumb"]
        if max(vertnumb) < 2: return None   #Nothing to be compressed           
#Loop over branches
        for ib,branch in enumerate(self.branches):
            if vertnumb[ib] < 2: continue
            masses = branch.masses
            for ivertex in range(vertnumb[ib]-1):
                if abs(masses[ivertex]-masses[ivertex+1]) < mingap:
                    newelement.branches[ib].particles.pop(ivertex)
                    newelement.branches[ib].masses.pop(ivertex)

        if newelement.isEqual(self):
            return None
        else:
            return newelement
        
    def invisibleCompress(self):
        """ compress cascade decays ending with neutrinos + daughter If no compression is
        possible, return None. """
        
        
        newelement = copy.deepcopy(self)        
        vertnumb = self.getEinfo()["vertnumb"]
        if max(vertnumb) < 2: return None   #Nothing to be compressed        
#Loop over branches
        for ib,branch in enumerate(self.branches):
            if vertnumb[ib] < 2: continue
            particles = branch.particles            
            for ivertex in range(vertnumb[ib]-2,-1,-1):
                if particles[ivertex].count('nu') == len(particles[ivertex]):
                    newelement.branches[ib].masses.pop(ivertex+1)
                    newelement.branches[ib].particles.pop(ivertex)
                else:
                    break

                    
        if newelement.isEqual(self):
            return None
        else:
            return newelement


def smallerMass(mass1,mass2):
    """
    Uses an ordering criterium (machine-independent) to select the smaller of the two mass arrays
    """
    mass1List = []
    mass2List = []
    if mass1 == mass2: return mass1
    try:    
        for branch in mass1: mass1List.extend(branch)
        for branch in mass2: mass2List.extend(branch)
        if len(mass1List) == len(mass2List):
            for im,m1 in enumerate(mass1List):
                if m1 < mass2List[im]: return mass1
                if m1 > mass2List[im]: return mass2
    except:
        pass
    
    logger.error('[smallerMasses]: invalid input')
    return False
        
        