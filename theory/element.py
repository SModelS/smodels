"""
.. module:: theory.element
   :synopsis: missing
    
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
    
"""

from ParticleNames import PtcDic, Reven, simParticles, elementsInStr
from branch import Branch
import crossSection
import logging

logger = logging.getLogger(__name__)


class Element(object):
    """
    Element class. Holds a pair of branches and the element weight
    (cross-section*BR).
    
    """
    def __init__(self, info=None):
        """
        Constructor.
        
        """
        self.branches = [Branch(),Branch()]
        self.weight = crossSection.XSectionList()
                
        if info:
            if type(info) == type(str()):   #Creates element from particle string
                elements = elementsInStr(info)
                if not elements or len(elements) > 1:
                    logging.error("Wrong input string "+info)
                    return False                
                else:
                    el = elements[0]
                    branches = elementsInStr(el[1:-1])
                    if not branches or len(branches) != 2:
                        logging.error("Wrong input string "+info)
                        return False 
                    self.branches = []
                    for branch in branches: self.branches.append(Branch(branch))
            elif type(info) == type([]) and type(info[0]) == type(Branch()): #Creates element from branch pair
                for ib,branch in enumerate(info): self.branches[ib] = branch.copy() 


    def __eq__ (self, other):
        return self.isEqual(other)
    
    
    def __ne__ (self, other):
        return not self.isEqual(other)    
    
    
    def __str__ ( self ):
        """
        Returns the canonical name of the element, e.g. [[jet],[jet]].
        
        """
        particleString = str(self.getParticles()).replace(" ","").replace("'","")
        return particleString


    def isEqual(self,other,order=False,useDict=True):
        """
        Compare two Elements. If all masses and particles are equal, returns
        True, otherwise returns False. If order = False, test both branch
        orderings (for an element doublet only). If useDict=True, allow for
        inclusive particle labels.
        
        """                           
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
        """
        Compare two Elements. If particles match, returns True, otherwise
        returns False. If order = False, test both branch orderings (for an
        element doublet only). If useDict=True, allow for inclusive particle
        labels.
        
        """        
        if type(self) != type(other): return False
        ptcs = self.getParticles()
        ptcsA = other.getParticles()
        if simParticles(ptcs,ptcsA,useDict): return True
        if not order:
            ptcsB = other.switchBranches().getParticles()   
            if simParticles(ptcs,ptcsB,useDict): return True
        
        return False
    
    
    def copy(self):
        """
        Creates a copy of itself (faster than deepcopy).
        
        """
        newel = Element()
        newel.branches = []
        for branch in self.branches: newel.branches.append(branch.copy())                
        newel.weight = self.weight.copy()
        return newel
    
    
    def setMasses(self,mass,same_order=True,oppos_order=False):
        """
        Sets the element masses to the input mass array. If same_order, set
        the masses to the same branch ordering. If oppos_order, set the masses
        to the opposite branch ordering. If both same_order and oppos_order,
        set the masses to the smaller of the two orderings.
        
        """                
        if same_order and oppos_order:
            newmass = smallerMass(mass)
        elif same_order:
            newmass = mass
        elif oppos_order:
            newmass = [mass[1],mass[0]]
        else:
            logger.error('Called with no possible ordering!')
            return False
        if len(newmass) != len(self.branches):
            logger.error('Called with wrong number of mass branches!')
            return False
               
        for i,mass in enumerate(newmass): self.branches[i].masses = mass[:]


    def switchBranches(self):
        """
        If the element contains a pair of branches, switches them.
                
        """
        newEl = self.copy()
        if len(self.branches) == 2: newEl.branches = [newEl.branches[1],newEl.branches[0]]
        return newEl


    def getParticles(self):
        """
        Returns the array of particles in the element.   
             
        """
        ptcarray = []
        for branch in self.branches: ptcarray.append(branch.particles)
        return ptcarray
    
    
    def getMasses(self):
        """
        Returns the array of masses in the element.   
             
        """
        massarray = []
        for branch in self.branches: massarray.append(branch.masses)
        return massarray
    
    
    def getDaughters(self):
        """
        Return the pair of daughter IDs (can be None, if the element does not
        have a definite daughter).      
           
        """        
        return (self.branches[0].daughterID,self.branches[1].daughterID)
    
    
    def getMothers(self):
        """
        Return the pair of mother IDs (can be None, if the element does not.
        have a mother daughter).
        
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
        """
        Returns the maximum of the two branch lengths.
        
        """        
        return max(self.branches[0].getLength(),self.branches[1].getLength())
    
    
    def hasTopInList(self,elementList):
        """
        Checks if the element topology matches any of the topologies in the
        element list.
        
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
        """
        Checks if the element is present in the element list. Ifigmass=False
        also check if the analysis has the element mass array.
        
        """
        for el in listOfElements:
            if igmass:
                if self.particlesMatch(el,useDict): return True
            else:
                if self.isEqual(el,useDict): return True

        return False
    
    
    def checkConsistency(self):
        """
        Checks if the particles defined in the element exist and are consistent
        with the element info.
        
        """
        info = self.getEinfo()
        for ib,branch in enumerate(self.branches):            
            for iv,vertex in enumerate(branch.particles):
                if len(vertex) != info['vertparts'][ib][iv]:
                    logger.error("Wrong syntax")
                    return False
                for ptc in vertex:
                    if not ptc in Reven.values() and not PtcDic.has_key(ptc):
                        logger.error("Unknown particle"+ptc)
                        return False
        return True    
    
    
    def compressElement(self,DoCompress,DoInvisible,minmassgap):
        """
        Keep compressing they original element and the derived ones till they
        can be compressed no more. Returns a list with the compressed elements.
        
        """
        added = True
        newElements = [self.copy()]   
        # Keep compressing the new topologies generated so far until no new compressions can happen:
        while added:
            added = False
            # Check for mass compressed topologies   
            if DoCompress:
                for element in newElements:             
                    newel = element.massCompress(minmassgap)
                    if newel and not newel.hasTopInList(newElements): #Avoids double counting (conservative)
                        newElements.append(newel) 
                        added = True
      
            # Check for invisible compressed topologies (look for effective LSP, such as LSP + neutrino = LSP')      
            if DoInvisible:
                for element in newElements:
                    newel = element.invisibleCompress()
                    if newel and not newel.hasTopInList(newElements): #Avoids double counting (conservative)
                        newElements.append(newel) 
                        added = True

                         
        newElements.pop(0)  # Remove original element                            
        return newElements
    
         
    def massCompress(self,mingap):
        """
        If two masses in this topology are degenerate, returns a compressed
        copy of the element. If no compression is possible, return None.
        
        """        
        newelement = self.copy()
        vertnumb = self.getEinfo()["vertnumb"]
        if max(vertnumb) < 2: return None   # Nothing to be compressed           
        # Loop over branches
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
        """
        Compress cascade decays ending with neutrinos + daughter If no
        compression is possible, return None.
        """
        newelement = self.copy()        
        vertnumb = self.getEinfo()["vertnumb"]
        if max(vertnumb) < 2: return None   # Nothing to be compressed        
        # Loop over branches
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
    Uses an ordering criterium (machine-independent) to select the smaller of
    the two mass arrays.
    
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
    
    logger.error('Invalid input')
    return False
