"""
.. module:: theory.element
   :synopsis: missing TODO write synopsis
    
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
    
"""

from smodels.theory.particleNames import ptcDic, rEven, simParticles, elementsInStr
from smodels.theory.branch import Branch
from smodels.theory import crossSection
from smodels.theory.printer import Printer
import logging

logger = logging.getLogger(__name__)


class Element(Printer):
    """
    An instance of this class represents an element.
    
    This class possesses a pair of branches and the element weight
    (cross-section * BR).
    
    """
    def __init__(self, info=None):
        self.branches = [Branch(), Branch()]
        self.weight = crossSection.XSectionList()
        self.motherElements = []
        """ Elements that arise from compression have mother elements.
            Mother elements are pairs of ( whence, element ),
            'whence' describing what the element is from 
            (mass compression, invisible compression, etc),
            while 'element' is the actual object.
            If element is not due to compression, 
            then list remains empty.
        """
        
        if info:
            # Create element from particle string
            if type(info) == type(str()):
                elements = elementsInStr(info)
                if not elements or len(elements) > 1:
                    nel = 0
                    if elements:
                        nel = len(elements)
                    logging.error("Malformed input string. Number of elements "
                                  "is %d, expected 1: in ``%s''", nel, info)
                    return False
                else:
                    el = elements[0]
                    branches = elementsInStr(el[1:-1])
                    if not branches or len(branches) != 2:
                        logging.error("Malformed input string. Number of "
                                      "branches is %d, expected 2: in ``%s''",
                                      len(branches), info)
                        return False
                    self.branches = []
                    for branch in branches:
                        self.branches.append(Branch(branch))
            # Create element from branch pair
            elif type(info) == type([]) and type(info[0]) == type(Branch()):
                for ib, branch in enumerate(info):
                    self.branches[ib] = branch.copy()

    def combineMotherElements ( self, el2 ):
        """ combine mother elements from self and el2 into self """
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

    def __eq__(self, other):
        return self.isEqual(other)


    def __hash__(self):
        return object.__hash__(self)


    def __ne__(self, other):
        return not self.isEqual(other)


    def __str__(self):
        """
        Create the canonical name of the element, e.g. [[jet],[jet]].
        
        """
        particleString = str(self.getParticles()).replace(" ", "").\
                replace("'", "")
        return particleString


    def isEqual(self, other, order=False, useDict=True):
        """
        Compare two Elements for equality.
        
        If order == False, test both branch orderings (for an element doublet
        only). If useDict == True, allow for nclusive particle labels.
        
        :returns: True, if all masses and particles are equal; False, else;        
        
        """
        if type(self) != type(other):
            return False
        mass = self.getMasses()
        massA = other.getMasses()

        if self.particlesMatch(other, order=True, useDict=useDict) \
                and mass == massA:
            return True
        if not order:
            otherB = other.switchBranches()
            massB = otherB.getMasses()
            if self.particlesMatch(otherB, order=True, useDict=useDict) \
                    and mass == massB:
                return True

        return False


    def particlesMatch(self, other, order=False, useDict=True):
        """
        Compare two Elements for matching particles.
        
        If order == False, test both branch orderings (for an element doublet
        only). If useDict == True, allow for inclusive particle abels.
        
        :returns: True, if particles match; False, else;
        
        """
        if type(self) != type(other):
            return False
        ptcs = self.getParticles()
        ptcsA = other.getParticles()
        if simParticles(ptcs, ptcsA, useDict):
            return True
        if not order:
            ptcsB = other.switchBranches().getParticles()
            if simParticles(ptcs, ptcsB, useDict):
                return True

        return False


    def copy(self):
        """
        Create a copy of self.
        
        Faster than deepcopy.
        
        """
        newel = Element()
        newel.branches = []
        for branch in self.branches:
            newel.branches.append(branch.copy())
        newel.weight = self.weight.copy()
        import copy
        # if len(self.motherElements)>=0:
        ## newel.motherElements = copy.deepcopy( self.motherElements ) 
        newel.motherElements = self.motherElements[:]
        return newel


    def setMasses(self, mass, sameOrder=True, opposOrder=False):
        """
        Set the element masses to the input mass array.
        
        If sameOrder == True, set the masses to the same branch ordering. If
        opposOrder == True, set the masses to the opposite branch ordering. If
        both sameOrder == True and opposOrder == True, set the masses to the
        smaller of the two orderings.
        
        """
        if sameOrder and opposOrder:
            if mass[0] == _smallerMass(mass[0], mass[1]):
                newmass = mass
            else:
                newmass = [mass[1], mass[0]]
        elif sameOrder:
            newmass = mass
        elif opposOrder:
            newmass = [mass[1], mass[0]]
        else:
            logger.error("Called with no possible ordering")
            return False
        if len(newmass) != len(self.branches):
            logger.error("Called with wrong number of mass branches")
            return False

        for i, mass in enumerate(newmass):
            self.branches[i].masses = mass[:]


    def switchBranches(self):
        """
        Switch branches, if the element contains a pair of them.
                
        """
        newEl = self.copy()
        if len(self.branches) == 2:
            newEl.branches = [newEl.branches[1], newEl.branches[0]]
        return newEl


    def getParticles(self):
        """
        Get the array of particles in the element.   
             
        """
        ptcarray = []
        for branch in self.branches:
            ptcarray.append(branch.particles)
        return ptcarray


    def getMasses(self):
        """
        Get the array of masses in the element.   
             
        """
        massarray = []
        for branch in self.branches:
            massarray.append(branch.masses)
        return massarray


    def getDaughters(self):
        """
        Get a pair of daughter IDs.
        
        Can be None, if the element does not have a definite daughter.      
           
        """
        return (self.branches[0].daughterID, self.branches[1].daughterID)


    def getMothers(self):
        """
        Get a pair of mother IDs.
        
        Can be None, if the element does not have a mother daughter.
        
        """
        return (self.branches[0].momID, self.branches[1].momID)


    def getEinfo(self):
        """
        Get global topology info from particle string.
        
        """
        vertnumb = []
        vertparts = []
        for branch in self.branches:
            vertnumb.append(len(branch.masses))
            vertparts.append([len(ptcs) for ptcs in branch.particles])
            if len(vertparts[len(vertparts) - 1]) == \
                    vertnumb[len(vertnumb) - 1] - 1:
                # Append 0 for stable LSP
                vertparts[len(vertparts) - 1].append(0)
        return {"vertnumb" : vertnumb, "vertparts" : vertparts}


    def _getLength(self):
        """
        Get the maximum of the two branch lengths.
        
        """
        return max(self.branches[0].getLength(), self.branches[1].getLength())


    def isInList(self, listOfElements, igmass=False, useDict=True):
        """
        Check if the element is present in the element list.
        
        If igmass == False also check if the analysis has the element mass
        array.
        
        """
        for el in listOfElements:
            if igmass:
                if self.particlesMatch(el, useDict):
                    return True
            else:
                if self.isEqual(el, useDict):
                    return True

        return False


    def checkConsistency(self):
        """
        Check if the particles defined in the element exist and are consistent
        with the element info.
        
        """
        info = self.getEinfo()
        for ib, branch in enumerate(self.branches):
            for iv, vertex in enumerate(branch.particles):
                if len(vertex) != info['vertparts'][ib][iv]:
                    logger.error("Wrong syntax")
                    return False
                for ptc in vertex:
                    if not ptc in rEven.values() and not ptc in ptcDic:
                        logger.error("Unknown particle" + ptc)
                        return False
        return True


    def compressElement(self, doCompress, doInvisible, minmassgap):
        """
        Keep compressing they original element and the derived ones till they
        can be compressed no more.
        
        :returns: list with the compressed elements
        
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


    def massCompress(self, mingap):
        """
        Perform mass compression.
        
        :returns: compressed copy of the element, if two masses in this
                  topology are degenerate; None, if compression is not possible;
        
        """
        newelement = self.copy()
        newelement.motherElements = [ ("mass", self.copy()) ]

        vertnumb = self.getEinfo()["vertnumb"]
        # Nothing to be compressed
        if max(vertnumb) < 2:
            return None
        # Loop over branches
        for ib, branch in enumerate(self.branches):
            if vertnumb[ib] < 2:
                continue
            masses = branch.masses
            for ivertex in range(vertnumb[ib] - 1):
                if abs(masses[ivertex] - masses[ivertex + 1]) < mingap:
                    newelement.branches[ib].particles[ivertex] = None
                    newelement.branches[ib].masses[ivertex] = None
            while newelement.branches[ib].particles.count(None) > 0:
                iNone = newelement.branches[ib].particles.index(None)
                newelement.branches[ib].particles.pop(iNone)
                newelement.branches[ib].masses.pop(iNone)


        if newelement.isEqual(self):
            return None
        else:
            return newelement


    def hasTopInList(self, elementList):
        """
        Check if the element topology matches any of the topologies in the
        element list.
        
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
        Compress cascade decays ending with neutrinos and daughter.
        
        If no compression is possible, return None.
        
        """
        newelement = self.copy()
        newelement.motherElements = [ ("invisible", self.copy()) ]

        vertnumb = self.getEinfo()["vertnumb"]
        # Nothing to be compressed
        if max(vertnumb) < 2:
            return None
        # Loop over branches
        for ib, branch in enumerate(self.branches):
            if vertnumb[ib] < 2:
                continue
            particles = branch.particles
            for ivertex in range(vertnumb[ib] - 2, -1, -1):
                if particles[ivertex].count('nu') == len(particles[ivertex]):
                    newelement.branches[ib].masses.pop(ivertex + 1)
                    newelement.branches[ib].particles.pop(ivertex)
                else:
                    break

        if newelement.isEqual(self):
            return None
        else:
            return newelement

    def formatData(self):
        """
        Select data preparation method through dynamic binding.
        """
        return Printer.formatElementData(self)


def _smallerMass(mass1, mass2):
    """
    Select the smaller of two mass arrays.
    
    Use an ordering criterion (machine-independent) for selection.
    
    """
    mass1List = []
    mass2List = []
    if mass1 == mass2:
        return mass1
    try:
        for branch in mass1:
            mass1List.extend(branch)
        for branch in mass2:
            mass2List.extend(branch)
        if len(mass1List) == len(mass2List):
            for im, m1 in enumerate(mass1List):
                if m1 < mass2List[im]:
                    return mass1
                if m1 > mass2List[im]:
                    return mass2
    except:  # TODO: except what?
        pass

    logger.error("Invalid input")
    return False
