"""
.. module:: branch
   :synopsis: Module holding the branch class and methods.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

import sys
from smodels.theory.particleNames import simParticles, elementsInStr
from smodels.tools.physicsUnits import fb, MeV
from smodels.particleClass import SMpdgs, SMparticles, BSMpdgs, ptcDict, getObjectFromPdg, getObjectFromName
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger


class Branch(object):
    """
    An instance of this class represents a branch.    
    A branch-element can be constructed from a string (e.g., ('[b,b],[W]').
    
    :ivar masses: list of masses for the intermediate states
    :ivar particles: list of particles (Particles objects) for the final states
    :ivar PIDs: a list of the pdg numbers of the intermediate states appearing in the branch.
                If the branch represents more than one possible pdg list, PIDs will correspond
                to a nested list (PIDs = [[pid1,pid2,...],[pidA,pidB,...]])
    :ivar maxWeight: weight of the branch (XSection object)
    """
    def __init__(self, info=None):
        """
        Initializes the branch. If info is defined, tries to generate
        the branch using it.
        
        :parameter info: string describing the branch in bracket notation
                         (e.g. [[e+],[jet]])
        """

        self.intmParticles = []
        self.particles = []
        
        self.maxWeight = None
        self.vertnumb = None
        self.vertparts = None
        
        if type(info) == type(str()):
            branch = elementsInStr(info)
            if not branch or len(branch) > 1:
                logger.error("Wrong input string " + info)
                raise SModelSError()
            else:
                branch = branch[0]
                vertices = elementsInStr(branch[1:-1])
                for vertex in vertices:
                    particleNames = vertex[1:-1].split(',')
                    for name in particleNames:
                        if type(getObjectFromName(name)) == list:
                            ptcs = [particle for particle in getObjectFromName(name)]
                        else: ptcs = getObjectFromName(name) 
                    # Syntax check:
                    for ptc in particleNames:
                        if not ptc in SMparticles \
                                and not ptc in ptcDict:
                            logger.error("Unknown particle. Add " + ptc + " to smodels/particleClass.py")
                            raise SModelSError()

                    self.particles.append(ptcs)

            self.vertnumb = len(self.particles)
            self.vertparts = [len(v) for v in self.particles]
        

    def __str__(self):
        """
        Create the branch bracket notation string, e.g. [[e+],[jet]].
        
        :returns: string representation of the branch (in bracket notation)    
        """
        
        st = str([[particle.label for particle in particleList ] for particleList in self.particles ]).replace("'", "")
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
        elif not self.intmParticles == other.intmParticles: #  masses
            comp = self.intmParticles > other.intmParticles
            if comp: return 1
            else: return -1
        else:
            return 0  #Branches are equal

    def sortParticles(self):
        """
        Sort the particles inside each vertex
        """
        
        for iv,vertex in enumerate(self.particles):
            self.particles[iv] = sorted(vertex, key=lambda x: x.label.lower())

    def setInfo(self):
        """
        Defines the number of vertices (vertnumb) and number of
        particles in each vertex (vertpats) properties, if they have not
        been defined yet.
        """

        self.vertnumb = len(self.particles)
        self.vertparts = [len(v) for v in self.particles]

    
    def particlesMatch(self, other):
        """
        Compare two Branches for matching particles, 
        allow for inclusive particle labels (such as the ones defined in particleClass.py)
        
        :parameter other: branch to be compared (Branch object)
        :returns: True if branches are equal (particles and masses match); False otherwise.              
        """
        
        if type (other) != type(self):
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
        newbranch = Branch()
        newbranch.particles = self.particles[:]
        newbranch.intmParticles = []
        self.setInfo()
        newbranch.vertnumb = self.vertnumb
        newbranch.vertparts = self.vertparts[:]
        for pidList in self.intmParticles:
            newbranch.intmParticles.append(pidList[:])
        if not self.maxWeight is None:
            newbranch.maxWeight = self.maxWeight.copy()
        return newbranch


    def getLength(self):
        """
        Returns the branch length (number of R-odd particles).
        
        :returns: length of branch (number of R-odd particles)
        """
        
        return len(self.intmParticles[0]) 

    def __lt__( self, b2 ):
        return self.__cmp__ ( b2 ) == -1

    def __eq__( self, b2 ):
        return self.__cmp__ ( b2 ) == 0

    def _addDecay(self, br):
        """
        Generate a new branch adding a 1-step cascade decay        
        This is described by the br object, with particle masses given by BSMList.
        
        :parameter br: Decay object (see pyslha). Contains information about the decay.
        :returns: extended branch (Branch object). False if there was an error.
        """
      
        newBranch = self.copy()
        newparticles = []
        newintmParticles = []

        if len(self.intmParticles[0]) != 1: #if len(self.PIDs[0]) != 1:
            logger.error("During decay the branch should not have multiple PID lists!")
            return False   

        for partID in br.ids:
            # Add R-even particles to final state
            if partID in SMpdgs:
                newparticles.append(getObjectFromPdg(partID))
            else:
                # Add non R-even particles to intermediate state    ;masses of non R-even particles to mass vector
                newintmParticles.append(getObjectFromPdg(partID))

        if len(newintmParticles) > 1: # newmass
            logger.warning("Multiple R-odd particles in the final state: " +
                           str(br.ids))
            return False

        if newparticles:
            newBranch.particles.append(sorted(newparticles, key=lambda x: x.label.lower()))

        if newintmParticles:
            newBranch.intmParticles[0].extend(newintmParticles)

        if not self.maxWeight is None:
            newBranch.maxWeight = self.maxWeight * br.br
        
        return newBranch


    def decayDaughter(self):
        """
        Generate a list of all new branches generated by the 1-step cascade
        decay of the current branch daughter.
        :returns: list of extended branches (Branch objects). Empty list if daughter is stable or
                  if daughterID was not defined.
        """   
        if len(self.intmParticles) != 1: # PIDs
            logger.error("Can not decay branch with multiple PID lists")
            return False                
        if not self.intmParticles[0][-1]: # PIDs
            # Do nothing if there is no R-odd daughter (relevant for RPV decays
            # of the LSP)
            return []
        #If decay table is not defined, assume daughter is stable:
        if not self.intmParticles[0][-1].pdg in BSMpdgs: return [] #PIDs without .pdg
        # List of possible decays (brs) for R-odd daughter in branch        
        brs = self.intmParticles[0][-1].branches  #getObjectFromPdg(self.PIDs[0][-1]).branches 
        if len(brs) == 0:
            # Daughter is stable, there are no new branches
            return []
        newBranches = []

        for br in brs:
            if not br.br: continue  #Skip zero BRs
            # Generate a new branch for each possible decay
            newBranches.append(self._addDecay(br))

        return newBranches


def decayBranches(branchList, sigcut=0. *fb):
    """
    Decay all branches from branchList until all unstable intermediate states have decayed.
    
    :parameter branchList: list of Branch() objects containing the initial mothers
    :parameter sigcut: minimum sigma*BR to be generated, by default sigcut = 0.
                   (all branches are kept)
    :returns: list of branches (Branch objects)    
    """

    finalBranchList = []
    while branchList:

        # Store branches after adding one step cascade decay
        newBranchList = []
        for inbranch in branchList:
            if sigcut.asNumber() > 0. and inbranch.maxWeight < sigcut:
                # Remove the branches above sigcut and with length > topmax
                continue
            # Add all possible decays of the R-odd daughter to the original
            # branch (if any)
            newBranches = inbranch.decayDaughter()
            if newBranches:
                # New branches were generated, add them for next iteration
                newBranchList.extend(newBranches)
            else:
                # All particles have already decayed, store final branch
                finalBranchList.append(inbranch)
        # Use new branches (if any) for next iteration step
        branchList = newBranchList
        
    #Sort list by initial branch PID:
    finalBranchList = sorted(finalBranchList, key=lambda branch:  branch.intmParticles[0][0].pdg)  
    return finalBranchList


