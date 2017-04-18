"""
.. module:: branch
   :synopsis: Module holding the branch class and methods.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

from smodels.theory.particleNames import simParticles, elementsInStr, getFinalStateLabel
from smodels.tools.physicsUnits import fb, MeV
from smodels.particles import rEven, ptcDic, finalStates
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger

class Branch(object):
    """
    An instance of this class represents a branch.    
    A branch-element can be constructed from a string (e.g., ('[b,b],[W]').
    
    :ivar masses: list of masses for the intermediate states
    :ivar particles: list of particles (strings) for the final states
    :ivar PIDs: a list of the pdg numbers of the intermediate states appearing in the branch.
                If the branch represents more than one possible pdg list, PIDs will correspond
                to a nested list (PIDs = [[pid1,pid2,...],[pidA,pidB,...])
    :ivar maxWeight: weight of the branch (XSection object)
    """
    def __init__(self, info=None, finalState=None):
        """
        Initializes the branch. If info is defined, tries to generate
        the branch using it.
        
        :parameter info: string describing the branch in bracket notation
                         (e.g. [[e+],[jet]])

        :parameter finalState: final state label string for the branch
                         (e.g. 'MET' or 'HSCP')                         
        """
        self.masses = []
        self.particles = []
        self.PIDs = []
        self.maxWeight = None
        self.vertnumb = None
        self.vertparts = None
        self.stable = False
        self.finalState = None
        if type(info) == type(str()):
            branch = elementsInStr(info)
            if not branch or len(branch) > 1:
                logger.error("Wrong input string " + info)
                raise SModelSError()
            else:
                branch = branch[0]
                vertices = elementsInStr(branch[1:-1])
                for vertex in vertices:
                    ptcs = vertex[1:-1].split(',')
                    # Syntax check:
                    for ptc in ptcs:
                        if not ptc in rEven.values() \
                                and not ptc in ptcDic:
                            logger.error("Unknown particle. Add " + ptc + " to smodels/particle.py")
                            raise SModelSError()
                    self.particles.append(ptcs)
            self.vertnumb = len(self.particles)
            self.vertparts = [len(v) for v in self.particles]
        
        self.setFinalState(finalState)


    def __str__(self):
        """
        Create the branch bracket notation string, e.g. [[e+],[jet]].
        
        :returns: string representation of the branch (in bracket notation)    
        """
        st = str(self.particles).replace("'", "")
        st = st.replace(" ", "")
        return st

    def __cmp__(self,other):
        """
        Compares the branch with other.        
        The comparison is made based on .
        OBS: The particles inside each vertex MUST BE sorted (see branch.sortParticles())         
        :param other:  branch to be compared (Branch object)
        :return: -1 if self < other, 0 if self == other, +1, if self > other.
        """
        
        # print ( type ( self.masses[0] ) )
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
        elif not self.masses == other.masses:
            comp = self.masses > other.masses
            if comp: return 1
            else: return -1
        elif self.finalState != other.finalState:
            comp = self.finalState > other.finalState
            if comp: return 1
            else: return -1            
        else:
            return 0  #Branches are equal

    def sortParticles(self):
        """
        Sort the particles inside each vertex
        """
        
        for iv,vertex in enumerate(self.particles):
            self.particles[iv] = sorted(vertex)

    def setInfo(self):
        """
        Defines the number of vertices (vertnumb) and number of
        particles in each vertex (vertpats) properties, if they have not
        been defined yet.
        """

        self.vertnumb = len(self.particles)
        self.vertparts = [len(v) for v in self.particles]
        
    def setFinalState(self,finalState=None):
        """
        If finalState = None, define the branch final state according to the PID of the
        last R-odd particle appearing in the cascade decay.
        Else set the final state to the finalState given
        :parameter finalState: String defining the final state
        """
        
        if finalState:
            if not finalState in finalStates:
                logger.error("Final state %s has not been defined. Add it to particles.py." %finalState)
                raise SModelSError
            else:
                self.finalState = finalState
        #If PIDs have been defined, use it:     
        elif self.PIDs:
            fStates = set()
            for pidList in self.PIDs:
                fStates.add(getFinalStateLabel(pidList[-1]))
            
            if len(fStates) != 1:
                logger.error("Error obtaining the final state for branch %s" %self)
                raise SModelSError
            else:
                self.finalState = list(fStates)[0]
        #Else do nothing
        else:
            self.finalState = None
    
    def particlesMatch(self, other):
        """
        Compare two Branches for matching particles, 
        allow for inclusive particle labels (such as the ones defined in particles.py).
        Includes the final state in the comparison. 
        
        :parameter other: branch to be compared (Branch object)
        :returns: True if branches are equal (particles and masses match); False otherwise.              
        """

        #First of all check final state:
        if self.finalState != other.finalState:
            return False

        #If particles are identical, avoid further checks
        if self.particles == other.particles:
            return True        
       
        if not isinstance(other,Branch):
            return False        
        
        #Make sure number of vertices and particles have been defined
        self.setInfo()
        other.setInfo()
        if self.vertnumb != other.vertnumb:
            return False
        if self.vertparts != other.vertparts:
            return False

        for iv,vertex in enumerate(self.particles):
            if not simParticles(vertex,other.particles[iv]):
                return False
              
        return True
   

    def copy(self):
        """
        Generate an independent copy of self.        
        Faster than deepcopy.
        
        :returns: Branch object
        """
        
        #Allows for derived classes (like wildcard classes)
        newbranch = self.__class__()
        newbranch.masses = self.masses[:]
        newbranch.particles = self.particles[:]
        newbranch.finalState = self.finalState
        newbranch.PIDs = []
        newbranch.stable = self.stable
        self.setInfo()
        newbranch.vertnumb = self.vertnumb
        newbranch.vertparts = self.vertparts[:]
        for pidList in self.PIDs:
            newbranch.PIDs.append(pidList[:])
        if not self.maxWeight is None:
            newbranch.maxWeight = self.maxWeight.copy()
        return newbranch


    def getLength(self):
        """
        Returns the branch length (number of R-odd particles).
        
        :returns: length of branch (number of R-odd particles)
        """
        
        return len(self.masses)

    def __lt__( self, b2 ):
        return self.__cmp__ ( b2 ) == -1

    def __eq__( self, b2 ):
        return self.__cmp__ ( b2 ) == 0

    def _addDecay(self, br, massDictionary):
        """
        Generate a new branch adding a 1-step cascade decay        
        This is described by the br object, with particle masses given by
        massDictionary.
        
        :parameter br: branching ratio object (see pyslha). Contains information about the decay.
        :parameter massDictionary: dictionary containing the masses for all intermediate states.
        :returns: extended branch (Branch object). False if there was an error.
        """
        newBranch = self.copy()
        newparticles = []
        newmass = []
        
        if len(self.PIDs) != 1:
            logger.error("During decay the branch should \
                            not have multiple PID lists!")
            return False
        

        for partID in br.ids:
            # Add R-even particles to final state
            if partID in rEven:
                newparticles.append(rEven[partID])
            else:
                # Add masses of non R-even particles to mass vector
                newmass.append(massDictionary[partID])
                newBranch.PIDs[0].append(partID)

        if len(newmass) > 1:
            logger.warning("Multiple R-odd particles in the final state: " +
                           str(br.ids))
            return False

        if newparticles:
            newBranch.particles.append(sorted(newparticles))
        if newmass:
            newBranch.masses.append(newmass[0])
        if not self.maxWeight is None:
            newBranch.maxWeight = self.maxWeight * br.br
        #If there are no daughters, assume branch is stable
        if not br.ids:
            newBranch.stable = True

        return newBranch


    def decayDaughter(self, brDictionary, massDictionary):
        """
        Generate a list of all new branches generated by the 1-step cascade
        decay of the current branch daughter.
        
        :parameter brDictionary: dictionary with the decay information
                                 for all intermediate states (values are br objects, see pyslha)
        :parameter massDictionary: dictionary containing the masses for all intermediate states.
        :returns: list of extended branches (Branch objects). Empty list if daughter is stable or
                  if daughterID was not defined.
        """

        if len(self.PIDs) != 1:
            logger.error("Can not decay branch with multiple PID lists")
            return False                
        if not self.PIDs[0][-1]:
            # Do nothing if there is no R-odd daughter (relevant for RPV decays
            # of the LSP)
            self.stable = True
            return [self]
        #If decay table is not defined, assume daughter is stable:
        if not self.PIDs[0][-1] in brDictionary:
            self.stable = True
            return [self]
        # List of possible decays (brs) for R-odd daughter in branch        
        brs = brDictionary[self.PIDs[0][-1]]       
        newBranches = []
        for br in brs:
            if not br.br:
                continue  #Skip zero BRs            
            # Generate a new branch for each possible decay (including "stable decays")
            newBranches.append(self._addDecay(br, massDictionary))
        
        if not newBranches:
            # Daughter is stable, there are no new branches
            self.stable = True
            return [self]
        else:                       
            return newBranches



class BranchWildcard(Branch):
    """
    A branch wildcard class. It will return True when compared to any other branch object.
    """
    
    def __init__(self):
        Branch.__init__(self)
        self.masses = ListWildcard()
        self.particles =  ListWildcard()
        self.PIDs = ListWildcard()
        self.vertnumb = IntWildcard()
        self.vertparts = ListWildcard()
        self.finalState = ListWildcard()        
        
    def __str__(self):
        return '[*]'

    def __cmp__(self,other):
        if isinstance(other,Branch):
            return 0
        else:
            return -1

    def __eq__(self,other):
        return self.__cmp__(other) == 0
    
    def setInfo(self):
        """
        Defines the number of vertices (vertnumb) and number of
        particles in each vertex (vertpats) properties, if they have not
        been defined yet.
        """

        self.vertnumb = IntWildcard()
        self.vertparts = ListWildcard()
        
class IntWildcard(int):
    """
    A integer wildcard class. It will return True when compared to any other integer object.
    """
    
    def __init__(self):
        int.__init__(self)
        
    def __str__(self):
        return '*'    

    def __cmp__(self,other):
        if isinstance(other,int):
            return 0
        else:
            return -1

    def __eq__(self,other):
        return self.__cmp__(other) == 0          
    
class ListWildcard(list):    
    """
    A list wildcard class. It will return True when compared to any other list object.
    """
    
    def __init__(self):
        int.__init__(self)
        
    def __str__(self):
        return '[*]'    

    def __cmp__(self,other):
        if isinstance(other,list):
            return 0
        else:
            return -1

    def __eq__(self,other):
        return self.__cmp__(other) == 0  

def decayBranches(branchList, brDictionary, massDictionary,
                  sigcut=0. *fb):
    """
    Decay all branches from branchList until all unstable intermediate states have decayed.
    
    :parameter branchList: list of Branch() objects containing the initial mothers
    :parameter brDictionary: dictionary with the decay information
                             for all intermediate states (values are br objects, see pyslha).
                            It may also contain information about long-lived particles.
    :parameter massDictionary: dictionary containing the masses for all intermediate states.
    :parameter promptDictionary: optional dictionary with the fraction of prompt and non-prompt decays.
                        Allows to deal with quasi-stable or long-lived particles
                        If not given, all particles are considered to
                        always decay promptly or to be stable.    
    :parameter sigcut: minimum sigma*BR to be generated, by default sigcut = 0.
                   (all branches are kept)
    :returns: list of branches (Branch objects)    
    """

    #Check for stable branches
    stableBranches = [branch for branch in branchList if branch.stable]
    unstableBranches = [branch for branch in branchList if not branch.stable]
        
    if not unstableBranches:
        #All branches have been decayed:
        finalBranches = sorted(stableBranches, key=lambda branch: branch.PIDs)
        return finalBranches
    else:
        #Decayed unstable branches
        newBranches = []
        for branch in unstableBranches:
            newBranches += [br for br in branch.decayDaughter(brDictionary, massDictionary) 
                           if br.maxWeight >= sigcut]
        
        return decayBranches(newBranches+stableBranches,brDictionary, massDictionary,sigcut)
