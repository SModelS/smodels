"""
.. module:: element
   :synopsis: Module holding the Element class and its methods.
    
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
    
"""

from smodels.theory.auxiliaryFunctions import elementsInStr
from smodels.theory.branch import Branch, InclusiveBranch
from smodels.theory import crossSection
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger

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
    def __init__(self, info=None, finalState=None):
        """
        Initializes the element. If info is defined, tries to generate
        the element using it.
        
        :parameter info: string describing the element in bracket notation
                         (e.g. [[[e+],[jet]],[[e-],[jet]]])
                         
        :parameter finalState: list containing the final state labels for each branch
                         (e.g. ['MET', 'HSCP'] or ['MET','MET'])
                         
        """
        self.branches = [Branch(), Branch()]
        self.weight = crossSection.XSectionList() # gives the weight for all decays promptly
        self.decayLabels = []
        self.motherElements = [("original", self)]
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
                    logger.error("Malformed input string. Number of elements "
                                  "is %d (expected 1) in: ``%s''", nel, info)
                    return None
                else:
                    el = elements[0]
                    branches = elementsInStr(el[1:-1])
                    if not branches or len(branches) != 2:
                        logger.error("Malformed input string. Number of "
                                      "branches is %d (expected 2) in: ``%s''",
                                      len(branches), info)
                        return None
                    self.branches = []
                    if finalState:
                        if len(finalState) != len(branches):
                            raise SModelSError("Number of final states (%i) does not match number of branches (%i)" 
                                               %(len(finalState),len(branches)))
                    else:
                        finalState = [None]*len(branches)                  
                    for ibr,branch in enumerate(branches):                        
                        if branch == '[*]':
                            self.branches.append(InclusiveBranch(finalState[ibr]))                           
                        else:
                            self.branches.append(Branch(branch,finalState[ibr])) 

            # Create element from branch pair
            elif isinstance(info,list) and all(isinstance(x,(Branch,InclusiveBranch)) for x in info):                
                self.branches = [br.copy() for br in info]
            else:
                raise SModelSError("Can not create element from input type %s" %type(info))
        
        self.setEinfo()
        
    
    def __cmp__(self,other):
        """
        Compares the element with other.        
        The comparison is made based on branches.
        OBS: The elements and the branches must be sorted! 
        :param other:  element to be compared (Element object)
        :return: -1 if self < other, 0 if self == other, +1, if self > other.
        """

        if not isinstance(other,Element):
            return -1

        #Compare branches:
        if self.branches != other.branches:          
            comp = self.branches > other.branches
            if comp: return 1
            else: return -1
        else:
            return 0     

    def __eq__(self,other):
        return self.__cmp__(other)==0

    def __lt__(self,other):
        return self.__cmp__(other)<0

    def __hash__(self):
        return object.__hash__(self)


    def __str__(self):
        """
        Create the element bracket notation string, e.g. [[[jet]],[[jet]].
        
        :returns: string representation of the element (in bracket notation)    
        """
        
        elStr = "["+",".join([str(br) for br in self.branches])+"]"
        elStr = elStr.replace(" ", "").replace("'", "")
        return elStr
    
    def __repr__(self):
        
        return self.__str__()

    def toStr(self):
        """
        Returns a string with the element represented in bracket notation,
        including the final states, e.g. [[[jet]],[[jet]] (MET,MET)
        """
        
        elStr = str(self)+' '+str(tuple(self.getFinalStates())).replace("'","")
        
        return elStr
    
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
        Allow for inclusive particle labels (such as the ones defined in finalStateParticles.py).
        If branchOrder = False, check both branch orderings.
        
        :parameter other: element to be compared (Element object)
        :parameter branchOrder: If False, check both orderings, otherwise
                                check the same branch ordering
        :returns: True, if particles match; False, else;        
        """

        if not isinstance(other,Element):
            return False

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

        #Allows for derived classes (like inclusive classes)
        newel = self.__class__()
        newel.branches = []
        for branch in self.branches:
            newel.branches.append(branch.copy())
        newel.weight = self.weight.copy()
        newel.motherElements = self.motherElements[:]
        newel.elID = self.elID
        return newel



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
        Get the array of even particle objects in the element.
        
        :returns: list of Particle objects                
        """
        
        particles = [branch.evenParticles for branch in self.branches]

        return particles
    
    def getFinalStates(self):
        """
        Get the array of final state (last BSM particle) particle objects in the element.
        
        :returns: list of Particle objects
        """
        
        fsparticles = [branch.BSMparticles[-1] for branch in self.branches]

        return fsparticles


    def getMasses(self):
        """
        Get the array of BSM masses in the element.    
        
        :returns: list of masses (mass array)            
        """        
        
        massarray = [branch.getMasses() for branch in self.branches]

        return massarray
        
    def getBSMparticles(self):
        """
        Get the list of BSM particles appearing the cascade decay,
        including the last (stable) one.

        :returns: list of Particle or ParticleList objects
        """    
        
        BSMparticles = []
        for branch in self.branches:
            BSMparticles.append([particle for particle in branch.BSMparticles])      

        return BSMparticles

    def getPIDs(self):
        """
        Get the list of IDs (PDGs of the intermediate states appearing the cascade decay), i.e.
        [  [[pdg1,pdg2,...],[pdg3,pdg4,...]] ].
        
        :returns: list of PDG ids
        """

        BSMpids = []
        for br in self.getBSMparticles():
            BSMpids.append([particle.pdg for particle in br])

        return BSMpids

    def getDaughters(self):
        """
        Get the list of daughter (last/stable BSM particle in the decay) PDGs. 
        Can be a nested list, if the element combines several daughters:
        [ [pdgDAUG1,pdgDAUG2],  [pdgDAUG1',pdgDAUG2']] 

        
        :returns: list of PDG ids
        """
        
        daughterPIDs = []
        for branch in self.branches():
            daughterPIDs.append(branch.BSMparticles[-1].pdg)
            
        return daughterPIDs
    
    def getMothers(self):
        """
        Get list of mother PDGs.    
        Can be a nested list, if the element combines several mothers:
        [ [pdgMOM1,pdgMOM2],  [pdgMOM1',pdgMOM2']] 
        
        :returns: list of PDG ids
        """                        

        momPIDs = []
        for branch in self.branches:
            momPIDs.append(branch.BSMparticles[0].pdg)        
        return momPIDs    
        
    def getMotherPIDs(self):
        """
        Get PIDs of all mothers.
        
        :returns: list of mother PDG ids
        """
        
        # get a list of all original mothers, i.e. that only consists of first generation elements
        allmothers = self.motherElements
        while not all(mother[0] == 'original' for mother in allmothers):                        
            for im, mother in enumerate(allmothers):
                if mother[0] != 'original':
                    allmothers.extend( mother[1].motherElements )
                    allmothers.pop(im)                                    

        # get the PIDs of all mothers      
        motherPids = []
        for mother in allmothers:  
            motherPids.append( mother[1].getPIDs() ) 

        branch1 = []
        branch2 = []
        for motherPid in motherPids:
            branch1.append(motherPid[0])
            branch2.append(motherPid[1])
        for pid1 in branch1:
            for pid2 in branch2:
                pids = [pid1,pid2]
                if not pids in motherPids: motherPids.append(pids)  
   
        return motherPids     
        


    def getEinfo(self):
        """
        Get element topology info from branch topology info.
        
        :returns: dictionary containing vertices and number of final states information  
        """
        
        vertnumb = []
        vertparts = []
        for branch in self.branches:
            if branch.vertnumb is None:
                branch.setInfo()
            bInfo = branch.getInfo()
            vertparts.append(bInfo['vertparts'])
            vertnumb.append(bInfo['vertnumb'])
                
        return {"vertnumb" : vertnumb, "vertparts" : vertparts}

    def setEinfo(self):
        """
        Set topology info for each branch.  
        """
        
        for branch in self.branches:
            branch.setInfo()


    def _getLength(self):
        """
        Get the maximum of the two branch lengths.    
        
        :returns: maximum length of the element branches (int)    
        """
        return max(self.branches[0].getLength(), self.branches[1].getLength())


    def checkConsistency(self):
        """
        Check if the particles defined in the element are consistent
        with the element info.
        
        :returns: True if the element is consistent. Print error message
                  and exits otherwise.
        """
        info = self.getEinfo()
        for ib, branch in enumerate(self.branches):
            for iv, vertex in enumerate(branch.evenParticles):
                if len(vertex) != info['vertparts'][ib][iv]:
                    logger.error("Wrong syntax")
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
        
        if not doCompress and not doInvisible:
            return []
        
        added = True
        newElements = [self]
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
    
    
    def removeVertex(self,ibr,iv):
        """
        Remove vertex iv located in branch ibr.
        The "vertex-mother" in BSMparticles and (SM) particles in the vertex
        are removed from the branch. The vertex index corresponds
        to the BSM decay (iv = 0 will remove the first BSM particle,...)
        
        :parameter ibr: Index of branch (int)
        :parameter iv: Index of vertex in branch ibr (int)
        
        """
        
        self.branches[ibr].removeVertex(iv)

    def massCompress(self, minmassgap):
        """
        Perform mass compression.
        
        :parameter minmassgap: value (in GeV) of the maximum 
                               mass difference for compression
                               (if mass difference < minmassgap -> perform mass compression)
        :returns: compressed copy of the element, if two masses in this
                  element are degenerate; None, if compression is not possible;        
        """

        newelement = self.copy()
        newelement.motherElements = [("mass", self)]
        #Loop over branches and look for small mass differences 
        for ibr,branch in enumerate(newelement.branches):
            #Get mass differences      
     
            massDiffs = []
            for i,mom in enumerate(branch.BSMparticles[:-1]):
                if isinstance(branch.BSMparticles[i+1].mass, list):
                    pmass = min(branch.BSMparticles[i+1].mass)
                else: pmass = branch.BSMparticles[i+1].mass
                massDiffs.append(mom.mass - pmass)
            #Get vertices which have deltaM < minmassgap:
            removeVertices = [iv for iv,m in enumerate(massDiffs) if m < minmassgap]
            #Remove vertices till they exist:
            while removeVertices:
                newelement.removeVertex(ibr,removeVertices[0])
                branch = newelement.branches[ibr]
                massDiffs = [mom.mass - branch.BSMparticles[i+1].mass for i,mom in enumerate(branch.BSMparticles[:-1])]
                removeVertices = [iv for iv,m in enumerate(massDiffs) if m < minmassgap]
                
        for ibr,branch in enumerate(newelement.branches):
            if branch.vertnumb != self.branches[ibr].vertnumb:
                newelement.sortBranches()
                return newelement
            
        #New element was not compressed, return None
        return None
    

    def invisibleCompress(self):
        """
        Perform invisible compression.
        
        :returns: compressed copy of the element, if element ends with invisible
                  particles; None, if compression is not possible
        """
        
        newelement = self.copy()
        newelement.motherElements = [("invisible", self)]

        # Loop over branches
        for branch in newelement.branches:
            if not branch.evenParticles:
                continue
            #Check if the last decay should be removed:            
            neutralSM = all(ptc.isMET() for ptc in branch.evenParticles[-1])
            neutralBSM = branch.BSMparticles[-2].isMET() 
            if neutralBSM and neutralSM:
                removeLastVertex = True
            else:
                removeLastVertex = False
                
            while len(branch.BSMparticles) > 1 and removeLastVertex:
                bsmMom = branch.BSMparticles[-2]
                branch.removeVertex(len(branch.BSMparticles)-2)
                #For invisible compression, keep the mother of the vertex and not the daughter:
                branch.BSMparticles[-1] = bsmMom
                #Re-check if the last decay should be removed:
                if not branch.evenParticles:
                    continue
                neutralSM = all(ptc.isMET() for ptc in branch.evenParticles[-1])
                neutralBSM = branch.BSMparticles[-2].isMET()
                if neutralBSM and neutralSM:
                    removeLastVertex = True
                else:
                    removeLastVertex = False

        for ibr,branch in enumerate(newelement.branches):
            if branch.vertnumb != self.branches[ibr].vertnumb:
                newelement.sortBranches()
                return newelement
            
        #New element was not compressed, return None
        return None


    def combineWith(self, other):
        """
        Combine two elements. Should only be used if the elements
        are considered as equal.
        The elements branches and weights are combined as well as their mothers. 
        
        :parameter other: element (Element Object)  
        """
        
        if self != other:
            raise SModelSError("Asked to combine distinct elements")

        self.motherElements += other.motherElements[:]
        self.weight.combineWith(other.weight)
        for ibr,branch in enumerate(self.branches):
            branch.combineWith(other.branches[ibr])
    
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

