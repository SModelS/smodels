"""
.. module:: branch
   :synopsis: Module holding the branch class and methods.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
        
"""

from smodels.theory.auxiliaryFunctions import elementsInStr,simParticles
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import fb
from smodels.theory.particle import ParticleList,ParticleWildcard
from smodels.tools.wildCards import ValueWildcard,ListWildcard
from smodels.experiment.finalStateParticles import finalStates



class Branch(object):
    """
    An instance of this class represents a branch.    
    A branch-element can be constructed from a string (e.g., ('[b,b],[W]').

    :ivar particles: list of particles (Particle objects) for the final states
    :ivar BSMparticles: a list of the intermediate states particles appearing in the branch.
                If the branch represents more than one possible particle list, BSMparticles will correspond
                to a nested list (BSMparticles = [[particle1, particle2,...],[particleA, particleB,...]])
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
        
        self.particles = []
        self.BSMparticles = []
        
        self.maxWeight = None
        self.vertnumb = None
        self.vertparts = None
        if isinstance(info,str):
            branch = elementsInStr(info)      
            if not branch or len(branch) > 1:
                raise SModelSError("Wrong input string " + info)
            else:                
                branch = branch[0]
                vertices = elementsInStr(branch[1:-1])        
                for vertex in vertices:
                    bsmParticle = finalStates.getParticlesWith(label='anyBSM')
                    if not bsmParticle:
                        raise SModelSError("Final state anyBSM has not been defined in finalStateParticles.py")
                    if len(bsmParticle) != 1:
                        raise SModelSError("Ambiguos defintion of label %s in finalStates" %bsmParticle[0].label)          
                    self.BSMparticles.append(bsmParticle[0])
                    particleNames = vertex[1:-1].split(',')
                    ptcs = []
                    for pname in particleNames:
                        smParticle = finalStates.getParticlesWith(label=pname)
                        if not smParticle:
                            raise SModelSError("Final state %s has not been defined in finalStateParticles.py " %pname)
                        ptcs.append(smParticle[0])
                    self.particles.append(ptcs)

            self.vertnumb = len(self.particles)
            self.vertparts = [len(v) for v in self.particles]
        if finalState:
            bsmParticle = finalStates.getParticlesWith(label=finalState)
            if not bsmParticle:
                raise SModelSError("Final state %s has not been defined in finalStateParticles.py")
            if len(bsmParticle) != 1:
                raise SModelSError("Ambiguos defintion of label %s in finalStateParticles.py" %finalState)
            self.BSMparticles.append(bsmParticle[0])

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
        The comparison is made based on vertnumb, vertparts, particles, masses of BSM particles
        and the last BSM particle appearing in the cascade decay.
        OBS: The particles inside each vertex MUST BE sorted (see branch.sortParticles())         
        :param other:  branch to be compared (Branch object)
        :return: -1 if self < other, 0 if self == other, +1, if self > other.
        """


        if not isinstance(other,(Branch,BranchWildcard)):
            return -1
        
        elif isinstance(other,BranchWildcard):
            return other.__cmp__(self)
        
        
        if self.vertnumb != other.vertnumb:
            comp = self.vertnumb > other.vertnumb
            if comp: return 1
            else: return -1
        elif not self.vertparts == other.vertparts:
            comp = self.vertparts > other.vertparts
            if comp: return 1
            else: return -1
            
        #Compare SM final states by label and Z2parity:
        for iv,vertex in enumerate(self.particles):            
            for iptc,particle in enumerate(vertex):
                comp = particle.cmpProperties(other.particles[iv][iptc],properties=['Z2parity','label'])
                if comp:
                    return comp 

        #Compare BSM states by Z2parity and mass:
        for iptc,bsmParticle in enumerate(self.BSMparticles):
            comp = bsmParticle.cmpProperties(other.BSMparticles[iptc],properties=['Z2parity','mass'])
            if comp:
                return comp
            
        #Compare final BSM state by Z2parity and quantum numbers:
        comp = self.BSMparticles[-1].cmpProperties(other.BSMparticles[-1],
                                                   properties=['Z2parity','colordim','eCharge'])
        if comp:
            return comp 

        return 0  #Branches are equal    
        

    def __lt__( self, b2 ):
        return self.__cmp__(b2) == -1

    def __eq__( self, b2 ):
        return self.__cmp__(b2) == 0
          
    def __ne__( self, b2 ):
        return not self.__cmp__(b2) == 0

        
    def getMasses(self):
        """
        Return list with masses of the BSM particles appearing in the branch
        
        :return: List with masses
        """          
        bsmMasses = [bsm.mass for bsm in self.BSMparticles]
        
        return bsmMasses             


    def sortParticles(self):
        """
        Sort the particles inside each vertex
        """
        
        for iv,vertex in enumerate(self.particles):
            self.particles[iv] = sorted(vertex, key=lambda x: x.label)


    def setInfo(self):
        """
        Defines the number of vertices (vertnumb) and number of
        particles in each vertex (vertpats) properties, if they have not
        been defined yet.
        """

        self.vertnumb = len(self.particles)
        self.vertparts = [len(v) for v in self.particles]
        
        
    def removeVertex(self,iv):
        """
        Remove vertex iv.
        The "vertex-mother" in BSMparticles and (SM) particles in the vertex
        are removed from the branch. The vertex index corresponds
        to the BSM decay (iv = 0 will remove the first BSM particle,...)
        
        :parameter iv: Index of vertex in branch (int)
        
        """
        
        self.BSMparticles = self.BSMparticles[:iv] + self.BSMparticles[iv+1:]
        self.particles = self.particles[:iv] + self.particles[iv+1:]
        self.setInfo()
        

    
    def particlesMatch(self, other):
        """
        Compare two Branches for matching particles.
        Only the final state particles (self.particles and self.BSMparticles[-1])
        are used for comparison.
        Allow for inclusive particle labels (such as the ones defined in particleDefinitions.py)

        :parameter other: branch to be compared (Branch object)
        
        :returns: True if branches are equal (particles and last BSM particle match); False otherwise.              
        """
        
       
        if not isinstance(other,Branch):
            return False        

        #Make sure number of vertices and particles have been defined
        self.setInfo()
        other.setInfo()
        if self.vertnumb != other.vertnumb:
            return False
        
        if self.vertparts != other.vertparts:
            return False
        
        #If particles are identical, avoid further checks
        if self.particles != other.particles:
            for iv,vertex in enumerate(self.particles):
                if not simParticles(vertex,other.particles[iv]):
                    return False

        #If the final state BSM particles have the same quantum numbers: 
        if not self.BSMparticles[-1].eqProperties(other.BSMparticles[-1],properties=['colordim','eCharge']):
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
        newbranch.particles = self.particles[:]
        newbranch.BSMparticles = self.BSMparticles[:]
        self.setInfo()
        newbranch.vertnumb = self.vertnumb
        newbranch.vertparts = self.vertparts[:]
        if not self.maxWeight is None:
            newbranch.maxWeight = self.maxWeight.copy()
        return newbranch


    def getLength(self):
        """
        Returns the branch length (number of R-odd particles).
        
        :returns: length of branch (number of R-odd particles)
        """
        
        return len(self.BSMparticles)
    
    
    def combineWith(self,other):
        """
        Combines itself with the other branch.
        Should only be used if both branches are considered equal.
        The BSM particles appearing in both branches are combined
        into ParticleList objects.
        
        :parameter other: branch (Branch Object)        
        """

        if self != other:
            raise SModelSError("Asked to combine distinct branches")
        
        for iptc,bsm in enumerate(other.BSMparticles):
            
            if bsm is self.BSMparticles[iptc]:
                continue #Particles are the same (do nothing)
            
            #Else create a particle list with particles from both
            if not isinstance(self.BSMparticles[iptc],ParticleList):
                bsmList = ParticleList(label = 'BSM (combined)', particles=[self.BSMparticles[iptc]])
            else:
                bsmList = self.BSMparticles[iptc]
                
            #Now combine particles from both branches:
            if not isinstance(bsm,ParticleList):
                bsmList.particles.append(bsm)
            else:
                bsmList.particles += bsm.particles
            
            #Finally, remove duplicates:
            newList = []
            for ptc in bsmList.particles:
                if any(ptc is x for x in newList):
                    continue
                newList.append(ptc)
            bsmList.particles = newList
            
            self.BSMparticles[iptc] = bsmList
                        


    def _addDecay(self, decay):
        """
        Generate a new branch adding a 1-step cascade decay        
        This is described by the br object, with particle masses given by BSMList.
        
        :parameter decay: Decay object (see pyslha). Contains information about the decay.
        :returns: extended branch (Branch object). False if there was an error.
        """
        
        newBranch = self.copy()      
        particles = [ptc for ptc in decay.daughters]
        oddParticles = [p for p in particles if p.Z2parity == 'odd']
        evenParticles = [p for p in particles if p.Z2parity == 'even']
        
        if len(oddParticles) != 1:
            logger.warning("Decay %s does not preserve Z2 and will be ignored" %str(decay))
            return False
        
        newBranch.BSMparticles.append(oddParticles[0])
        evenParticles = sorted(evenParticles, key=lambda x: x.label.lower())
        newBranch.particles.append(evenParticles)
        
        if not self.maxWeight is None:
            newBranch.maxWeight =  self.maxWeight*decay.br             
        
        newBranch.setInfo()
        return newBranch


    def decayDaughter(self):
        """
        Generate a list of all new branches generated by the 1-step cascade
        decay of the current branch daughter.
        :returns: list of extended branches (Branch objects). Empty list if daughter is stable or
                  if daughterID was not defined.
        """   
        
        if not self.BSMparticles or not self.BSMparticles[-1].decays: 
            return False

        newBranches = []
        for decay in self.BSMparticles[-1].decays:
            if not decay or not decay.br:
                continue  #Skip decay = None and zero BRs
            # Generate a new branch for each possible decay:
            newBr = self._addDecay(decay)
            if newBr:
                newBranches.append(newBr)
        
        if not newBranches:
            return False
        else:                       
            return newBranches



class BranchWildcard(Branch):
    """
    A branch wildcard class. It will return True when compared to any other branch object
    with the same final state.
    """
    
    def __init__(self,finalState=None):
        Branch.__init__(self)
        self.masses = ListWildcard()
        self.particles =  ListWildcard()
        if finalState:
            bsmParticle = finalStates.getParticlesWith(label=finalState)
            if not bsmParticle:
                raise SModelSError("Final state %s has not been defined in finalStateParticles.py" %finalState)
            if len(bsmParticle) != 1:
                raise SModelSError("Ambiguos defintion of label %s in finalStates" %bsmParticle[0].label)          
            self.BSMparticles = [bsmParticle[0]]
        else:
            self.BSMparticles = [ParticleWildcard()]
        self.vertnumb = ValueWildcard()
        self.vertparts = ListWildcard()
        
    def __cmp__(self,other):
        """
        Always returns true. The only exception is if a final state particle has been
        defined. In this case, will include the final state in the comparison.        
        The comparison is made based on vertnumb, vertparts, particles, masses of BSM particles
        and the last BSM particle appearing in the cascade decay.
        OBS: The particles inside each vertex MUST BE sorted (see branch.sortParticles())         
        :param other:  branch to be compared (Branch object)
        :return: -1 if self < other, 0 if self == other, +1, if self > other.
        """


        if not isinstance(other,(Branch,BranchWildcard)):
            return -1
        
        #If BSM particles are identical, avoid further checks                
        if self.BSMparticles and other.BSMparticles:
            #Compare final BSM state by Z2parity and quantum numbers:
            comp = self.BSMparticles[-1].cmpProperties(other.BSMparticles[-1],
                                                       properties=['Z2parity','colordim','eCharge'])
            if comp:
                return comp 
    
        return 0  #Branches are equal    
        
        
       
    def __str__(self):
        return '[*]'
    
    
    def getMasses(self):
        return ListWildcard()
    
    def setInfo(self):
        """
        Defines the number of vertices (vertnumb) and number of
        particles in each vertex (vertpats) properties, if they have not
        been defined yet.
        """

        self.vertnumb = ValueWildcard()
        self.vertparts = ListWildcard()
        
    def decayDaughter(self):
        """
        Always return False.
        """
        
        return False
        

def decayBranches(branchList, sigcut=0.*fb):
    """
    Decay all branches from branchList until all unstable intermediate states have decayed.
    
    :parameter branchList: list of Branch() objects containing the initial mothers
    :parameter sigcut: minimum sigma*BR to be generated, by default sigcut = 0.
                   (all branches are kept)
    :returns: list of branches (Branch objects)    
    """
    
        
    stableBranches,unstableBranches = [],[]
    
    for br in branchList:
        if br.maxWeight < sigcut:
            continue
        
        if br.decayDaughter():
            unstableBranches.append(br)
        else:
            stableBranches.append(br)
    
    while unstableBranches:        
        # Store branches after adding one step cascade decay
        newBranchList = []
        for inbranch in unstableBranches:
            if sigcut.asNumber() > 0. and inbranch.maxWeight < sigcut:
                # Remove the branches above sigcut and with length > topmax
                continue
            
            #If None appear amongst the decays, add the possibility for the particle not decaying prompt
            if None in inbranch.BSMparticles[-1].decays:
                stableBranches.append(inbranch)
            
            # Add all possible decays of the R-odd daughter to the original
            # branch (if any)
            newBranches = inbranch.decayDaughter()
            if newBranches:
                # New branches were generated, add them for next iteration
                newBranchList += [br for br in newBranches if br.maxWeight >= sigcut]
            elif inbranch.maxWeight >= sigcut:
                stableBranches.append(inbranch)

        # Use new unstable branches (if any) for next iteration step
        unstableBranches = newBranchList           
    
    #Sort list by initial branch pdg:        
    finalBranchList = sorted(stableBranches, key=lambda branch: branch.BSMparticles[0].pdg)      

    return finalBranchList

