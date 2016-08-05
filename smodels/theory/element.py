"""
.. module:: theory.element
   :synopsis: Module holding the Element class and its methods.
    
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
    
"""

from smodels.theory.particleNames import elementsInStr
from smodels.theory.branch import Branch
from smodels.theory import crossSection
from smodels.particles import rEven, ptcDic
import logging
from smodels.theory.exceptions import SModelSTheoryError as SModelSError

logger = logging.getLogger(__name__)


class Element(object):
    """
    An instance of this class represents an element.    
    This class possesses a pair of branches and the element weight
    (cross-section * BR).
    
    :ivar branches: list of branches (Branch objects)
    :ivar weight: element weight (cross-section * BR)
    :ivar motherElements: only for elements generated from a parent element
                          by mass compression, invisible compression,etc.
                          Holds a pair of (whence, mother element), where
                          whence describes what process generated the element    
    """
    def __init__(self, info=None):
        """
        Initializes the element. If info is defined, tries to generate
        the element using it.
        
        :parameter info: string describing the element in bracket notation
                         (e.g. [[[e+],[jet]],[[e-],[jet]]])
        """
        self.branches = [Branch(), Branch()]
        self.weight = crossSection.XSectionList()
        self.motherElements = []
        self.elID = 0
        self.covered = False
        self.tested = False 
        
        if info:
            # Create element from particle string
            if type(info) == type(str()):
                elements = elementsInStr(info)
                if not elements or len(elements) > 1:
                    nel = 0
                    if elements:
                        nel = len(elements)
                    logging.error("Malformed input string. Number of elements "
                                  "is %d (expected 1) in: ``%s''", nel, info)
                    return False
                else:
                    el = elements[0]
                    branches = elementsInStr(el[1:-1])
                    if not branches or len(branches) != 2:
                        logging.error("Malformed input string. Number of "
                                      "branches is %d (expected 2) in: ``%s''",
                                      len(branches), info)
                        return False
                    self.branches = []                    
                    for branch in branches:
                        self.branches.append(Branch(branch))
            # Create element from branch pair
            elif type(info) == type([]) and type(info[0]) == type(Branch()):
                for ib, branch in enumerate(info):
                    self.branches[ib] = branch.copy()
        
    
    def __cmp__(self,other):
        """
        Compares the element with other.        
        The comparison is made based on branches.
        OBS: The elements and the branches must be sorted! 
        :param other:  element to be compared (Element object)
        :return: -1 if self < other, 0 if self == other, +1, if self > other.
        """
        
        #Compare branches:
        if self.branches != other.branches:            
            comp = self.branches > other.branches
            if comp: return 1
            else: return -1
        else:
            return 0     


    def __hash__(self):
        return object.__hash__(self)


    def __str__(self):
        """
        Create the element bracket notation string, e.g. [[[jet]],[[jet]].
        
        :returns: string representation of the element (in bracket notation)    
        """
        particleString = str(self.getParticles()).replace(" ", "").\
                replace("'", "")
        return particleString
    
    def sortBranches(self):
        """
        Sort branches. The smallest branch is the first one.
        See the Branch object for definition of branch size and comparison
        """
        
        #First make sure each branch is individually sorted 
        #(particles in each vertex are sorted)
        for br in self.branches:
            br.sortParticles()
        #Now sort branches
        self.branches = sorted(self.branches)


    def particlesMatch(self, other, branchOrder=False):
        """
        Compare two Elements for matching particles only.
        Allow for inclusive particle labels (such as the ones defined in particles.py).
        If branchOrder = False, check both branch orderings.
        
        :parameter other: element to be compared (Element object)
        :parameter branchOrder: If False, check both orderings, otherwise
                                check the same branch ordering
        :returns: True, if particles match; False, else;        
        """
        
        
        if type(self) != type(other):
            return False
        
        if len(self.branches) != len(other.branches):
            return False

        #Check if particles inside each branch match in the correct order
        branchMatches = []
        for ib,br in enumerate(self.branches):
            branchMatches.append(br.particlesMatch(other.branches[ib]))
        if sum(branchMatches) == 2:
                return True
        elif branchOrder:
            return False
        else:       
        #Now check for opposite order
            for ib,br in enumerate(self.switchBranches().branches):
                if not br.particlesMatch(other.branches[ib]):
                    return False

            return True


    def copy(self):
        """
        Create a copy of self.        
        Faster than deepcopy.     
        
        :returns: copy of element (Element object)   
        """
        newel = Element()
        newel.branches = []
        for branch in self.branches:
            newel.branches.append(branch.copy())
        newel.weight = self.weight.copy()
        newel.motherElements = self.motherElements[:]
        newel.elID = self.elID
        return newel


    def setMasses(self, mass, sameOrder=True, opposOrder=False):
        """
        Set the element masses to the input mass array.
        
        
        :parameter mass: list of masses ([[masses for branch1],[masses for branch2]])
        :parameter sameOrder: if True, set the masses to the same branch ordering
                              If True and opposOrder=True, set the masses to the
                              smaller of the two orderings.
        :parameter opposOrder: if True, set the masses to the opposite branch ordering.
                               If True and sameOrder=True, set the masses to the
                               smaller of the two orderings.             
        """
        if sameOrder and opposOrder:
            newmass = sorted(mass)
        elif sameOrder:
            newmass = mass
        elif opposOrder:
            newmass = [mass[1], mass[0]]
        else:
            logger.error("Called with no possible ordering")            
            raise SModelSError()
        if len(newmass) != len(self.branches):
            logger.error("Called with wrong number of mass branches")
            raise SModelSError()

        for i, mass in enumerate(newmass):
            self.branches[i].masses = mass[:]


    def switchBranches(self):
        """
        Switch branches, if the element contains a pair of them.
        
        :returns: element with switched branches (Element object)                
        """
        
        newEl = self.copy()
        if len(self.branches) == 2:
            newEl.branches = [newEl.branches[1], newEl.branches[0]]
        return newEl


    def getParticles(self):
        """
        Get the array of particles in the element.
        
        :returns: list of particle strings                
        """
        
        ptcarray = []
        for branch in self.branches:
            ptcarray.append(branch.particles)
        return ptcarray


    def getMasses(self):
        """
        Get the array of masses in the element.    
        
        :returns: list of masses (mass array)            
        """
        massarray = []
        for branch in self.branches:
            massarray.append(branch.masses)
        return massarray

    def getPIDs(self):
        """
        Get the list of IDs (PDGs of the intermediate states appearing the cascade decay), i.e.
        [  [[pdg1,pdg2,...],[pdg3,pdg4,...]] ].
        The list might have more than one entry if the element combines different pdg lists:
        [  [[pdg1,pdg2,...],[pdg3,pdg4,...]],  [[pdg1',pdg2',...],[pdg3',pdg4',...]], ...]
        
        :returns: list of PDG ids
        """
        
        pids = []
        for ipid,PIDlist in enumerate(self.branches[0].PIDs):            
            for ipid2,PIDlist2 in enumerate(self.branches[1].PIDs):
                pids.append([self.branches[0].PIDs[ipid],self.branches[1].PIDs[ipid2]])
        
        return pids

    def getDaughters(self):
        """
        Get a pair of daughter IDs (PDGs of the last intermediate 
        state appearing the cascade decay), i.e. [ [pdgLSP1,pdgLSP2] ]    
        Can be a list, if the element combines several daughters:
        [ [pdgLSP1,pdgLSP2],  [pdgLSP1',pdgLSP2']] 
        
        :returns: list of PDG ids
        """
        
        pids = self.getPIDs()
        daughterPIDs = []
        for pidlist in pids:
            daughterPIDs.append([pidlist[0][-1],pidlist[1][-1]])
            
        return daughterPIDs
    
    def getMothers(self):
        """
        Get a pair of mother IDs (PDGs of the first intermediate 
        state appearing the cascade decay), i.e. [ [pdgMOM1,pdgMOM2] ]    
        Can be a list, if the element combines several mothers:
        [ [pdgMOM1,pdgMOM2],  [pdgMOM1',pdgMOM2']] 
        
        :returns: list of PDG ids
        """
        
        pids = self.getPIDs()
        momPIDs = []
        for pidlist in pids:
            momPIDs.append([pidlist[0][0],pidlist[1][0]])
            
        return momPIDs    


    def getEinfo(self):
        """
        Get topology info from particle string.
        
        :returns: dictionary containing vertices and number of final states information  
        """
        
        vertnumb = []
        vertparts = []
        for branch in self.branches:
            vertparts.append([len(ptcs) for ptcs in branch.particles])
            vertnumb.append(len(branch.particles))
                
        return {"vertnumb" : vertnumb, "vertparts" : vertparts}


    def _getLength(self):
        """
        Get the maximum of the two branch lengths.    
        
        :returns: maximum length of the element branches (int)    
        """
        return max(self.branches[0].getLength(), self.branches[1].getLength())


    def checkConsistency(self):
        """
        Check if the particles defined in the element exist and are consistent
        with the element info.
        
        :returns: True if the element is consistent. Print error message
                  and exits otherwise.
        """
        info = self.getEinfo()
        for ib, branch in enumerate(self.branches):
            for iv, vertex in enumerate(branch.particles):
                if len(vertex) != info['vertparts'][ib][iv]:
                    logger.error("Wrong syntax")
                    raise SModelSError()
                for ptc in vertex:
                    if not ptc in rEven.values() and not ptc in ptcDic:
                        logger.error("Unknown particle. Add " + ptc + " to smodels/particle.py")
                        raise SModelSError()
        return True

    
    def compressElement(self, doCompress, doInvisible, minmassgap):
        """
        Keep compressing the original element and the derived ones till they
        can be compressed no more.
        
        :parameter doCompress: if True, perform mass compression
        :parameter doInvisible: if True, perform invisible compression
        :parameter minmassgap: value (in GeV) of the maximum 
                               mass difference for compression
                               (if mass difference < minmassgap, perform mass compression)
        :returns: list with the compressed elements (Element objects)        
        """
        added = True
        newElements = [self.copy()]
        # Keep compressing the new topologies generated so far until no new
        # compressions can happen:
        while added:
            added = False
            # Check for mass compressed topologies
            if doCompress:
                for element in newElements:
                    newel = element.massCompress(minmassgap)
                    # Avoids double counting (conservative)
                    if newel and not newel.hasTopInList(newElements):
                        newElements.append(newel)
                        added = True

            # Check for invisible compressed topologies (look for effective
            # LSP, such as LSP + neutrino = LSP')
            if doInvisible:
                for element in newElements:
                    newel = element.invisibleCompress()
                    # Avoids double counting (conservative)
                    if newel and not newel.hasTopInList(newElements):
                        newElements.append(newel)
                        added = True

        newElements.pop(0)  # Remove original element
        return newElements

    def massCompress(self, minmassgap):
        """
        Perform mass compression.
        
        :parameter minmassgap: value (in GeV) of the maximum 
                               mass difference for compression
                               (if mass difference < minmassgap -> perform mass compression)
        :returns: compressed copy of the element, if two masses in this
                  element are degenerate; None, if compression is not possible;        
        """
        
        masses = self.getMasses()
        massDiffs = []
        #Compute mass differences in each branch
        for massbr in masses:
            massDiffs.append([massbr[i]-massbr[i+1] for i in range(len(massbr)-1)])
        #Compute list of vertices to be compressed in each branch            
        compVertices = []
        for ibr,massbr in enumerate(massDiffs):
            compVertices.append([])
            for iv,massD in enumerate(massbr):            
                if massD < minmassgap: compVertices[ibr].append(iv)
        if not sum(compVertices,[]): return None #Nothing to be compressed
        else:
            newelement = self.copy()
            newelement.motherElements = [ ("mass", self.copy()) ]
            for ibr,compbr in enumerate(compVertices):
                if compbr:            
                    new_branch = newelement.branches[ibr]
                    ncomp = 0
                    for iv in compbr:
                        new_branch.masses.pop(iv-ncomp)
                        new_branch.particles.pop(iv-ncomp)
                        ncomp +=1
                    new_branch.setInfo() 

        newelement.sortBranches()
        return newelement
    

    def hasTopInList(self, elementList):
        """
        Check if the element topology matches any of the topologies in the
        element list.
        
        :parameter elementList: list of elements (Element objects)
        :returns: True, if element topology has a match in the list, False otherwise.        
        """
        if type(elementList) != type([]) or len(elementList) == 0:
            return False
        for element in elementList:
            if type(element) != type(self):
                continue
            info1 = self.getEinfo()
            info2 = element.getEinfo()
            info2B = element.switchBranches().getEinfo()
            if info1 == info2 or info1 == info2B:
                return True
        return False


    def invisibleCompress(self):
        """
        Perform invisible compression.
        
        :returns: compressed copy of the element, if element ends with invisible
                  particles; None, if compression is not possible
        """
        newelement = self.copy()
        newelement.motherElements = [ ("invisible", self.copy()) ]

        # Loop over branches
        for ib, branch in enumerate(self.branches):
            particles = branch.particles
            if not branch.particles:
                continue # Nothing to be compressed
            #Go over the branch starting at the end and remove invisible vertices: 
            for ivertex in reversed(range(len(particles))):
                if particles[ivertex].count('nu') == len(particles[ivertex]):
                    newelement.branches[ib].masses.pop(ivertex+1)
                    newelement.branches[ib].particles.pop(ivertex)
                else:
                    break
            newelement.branches[ib].setInfo()

        newelement.sortBranches()
        if newelement == self:
            return None
        else:            
            return newelement


    def combineMotherElements ( self, el2 ):
        """
        Combine mother elements from self and el2 into self
        
        :parameter el2: element (Element Object)  
        """
        if len(self.motherElements)==0: 
            # no mothers? then you yourself are mother!
            tmp=self.copy()
            self.motherElements.append ( ("combine", tmp) )
        for m in el2.motherElements:
            self.motherElements.append ( (m[0], m[1].copy()) )
        if len(el2.motherElements)==0: 
            # no mothers? then yo yourself are mother now
            tmp=el2.copy()
            self.motherElements.append ( ("combine", tmp) )


    def combinePIDs(self,el2):
        """
        Combine the PIDs of both elements. If the PIDs already appear in self,
        do not add them to the list.
        
        :parameter el2: element (Element Object) 
        """
           
        elPIDs = self.getPIDs()
        newelPIDs = el2.getPIDs()
        for pidlist in newelPIDs:                    
            if not pidlist in elPIDs:
                self.branches[0].PIDs.append(pidlist[0])
                self.branches[1].PIDs.append(pidlist[1])
