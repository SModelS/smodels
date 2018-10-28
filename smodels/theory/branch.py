"""
.. module:: branch
   :synopsis: Module holding the branch class and methods.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
        
"""

from smodels.theory.auxiliaryFunctions import elementsInStr
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import fb
from smodels.theory.particle import MultiParticle, Particle
from smodels.tools.inclusiveObjects import InclusiveValue,InclusiveList
from smodels.experiment.finalStateParticles import finalStates,anyOdd
import itertools



class Branch(object):
    """
    An instance of this class represents a branch.    
    A branch-element can be constructed from a string (e.g., ('[b,b],[W]').

    :ivar evenParticles: list of even particles (Particle objects) for the final states
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
        
        self.evenParticles = []
        self.oddParticles = []
        
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
                    bsmParticle = finalStates.getParticlesWith(label='anyOdd')
                    if not bsmParticle:
                        raise SModelSError("Final state anyOdd has not been defined in finalStateParticles.py")
                    elif len(bsmParticle) != 1:
                        raise SModelSError("Ambiguos defintion of label %s in finalStates" %bsmParticle[0].label)          
                    self.oddParticles.append(bsmParticle[0])
                    particleNames = vertex[1:-1].split(',')
                    ptcs = []
                    for pname in particleNames:
                        smParticle = finalStates.getParticlesWith(label=pname)
                        if not smParticle:
                            raise SModelSError("Final state %s has not been defined in finalStateParticles.py " %pname)
                        elif len(smParticle) != 1:
                            raise SModelSError("Ambiguos defintion of label %s in finalStates" %smParticle[0].label)
                        else:
                            ptcs.append(smParticle[0])
                    self.evenParticles.append(ptcs)

            self.vertnumb = len(self.evenParticles)
            self.vertparts = [len(v) for v in self.evenParticles]
        if finalState:
            bsmParticle = finalStates.getParticlesWith(label=finalState)
            if not bsmParticle:
                raise SModelSError("Final state %s has not been defined in finalStateParticles.py" %finalState)
            elif len(bsmParticle) != 1:
                raise SModelSError("Ambiguos defintion of label %s in finalStateParticles.py" %finalState)
            else:
                bsmParticle = bsmParticle[0]
        else:
            bsmParticle = anyOdd
            
        self.oddParticles.append(bsmParticle)

    def __str__(self):
        """
        Create the branch bracket notation string, e.g. [[e+],[jet]].
        
        :returns: string representation of the branch (in bracket notation)    
        """

        st = str(self.evenParticles).replace("'", "")
        st = st.replace(" ", "")
        return st

    def __cmp__(self,other):
        """
        Compares the branch with other.        
        The comparison is made based on vertnumb, vertparts, evenParticles and oddParticles.
        The comparison allows for any ordering of the evenParticles in the vertex.
        It relies on the particle comparison, which allows for the comparison of Particles and MultiParticles.
        Only the properties which are defined for both particles are compared.
        :param other:  branch to be compared (Branch object)
        :return: -1 if self < other, 0 if self == other, +1, if self > other.
        """


        if not isinstance(other,(Branch,InclusiveBranch)):
            return -1
        
        elif isinstance(other,InclusiveBranch):
            return -1*other.__cmp__(self)
        
        
        if self.vertnumb != other.vertnumb:
            comp = self.vertnumb > other.vertnumb
            if comp: return 1
            else: return -1
        elif not self.vertparts == other.vertparts:
            comp = self.vertparts > other.vertparts
            if comp: return 1
            else: return -1
            
        #Compare even final states irrespective of ordering:
        for iv,particlesA in enumerate(self.evenParticles):
            particlesB = other.evenParticles[iv]
            equalVertex = False
            for pListA in itertools.permutations(particlesA):
                pListA = list(pListA)
                if pListA == particlesB:
                    equalVertex = True
                    break
            if not equalVertex:
                comp = (particlesA > particlesB) - (particlesA < particlesB)
                return comp
                
        #Compare BSM states:
        for iv,partA in enumerate(self.oddParticles):
            partB = other.oddParticles[iv]
            if partA != partB:
                comp = (partA > partB) - (partA < partB)
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
        bsmMasses = []    
        for bsm in self.oddParticles:
            if not isinstance(bsm.mass, list): bsmMasses.append(bsm.mass)
            else: bsmMasses.append(bsm.mass[0])
        #bsmMasses = [bsm.mass for bsm in self.oddParticles]
        
        return bsmMasses             

    def sortParticles(self):
        """
        Sort the particles inside each vertex
        """
        
        for iv,vertex in enumerate(self.evenParticles):
            self.evenParticles[iv] = sorted(vertex)

    def setInfo(self):
        """
        Defines the number of vertices (vertnumb) and number of
        even particles in each vertex (vertpats) properties, if they have not
        been defined yet.
        """

        bInfo = self.getInfo()
        self.vertnumb = bInfo['vertnumb']
        self.vertparts = bInfo['vertparts']
        
    def getInfo(self):
        """
        Get branch topology info from evenParticles.
        
        :returns: dictionary containing vertices and number of final states information  
        """

        vertnumb = len(self.evenParticles)
        vertparts = [len(v) for v in self.evenParticles]
        
        return {"vertnumb" : vertnumb, "vertparts" : vertparts}
        
    def removeVertex(self,iv):
        """
        Remove vertex iv.
        The "vertex-mother" in BSMparticles and (SM) particles in the vertex
        are removed from the branch. The vertex index corresponds
        to the BSM decay (iv = 0 will remove the first BSM particle,...)
        
        :parameter iv: Index of vertex in branch (int)
        
        """
        
        self.oddParticles = self.oddParticles[:iv] + self.oddParticles[iv+1:]
        self.evenParticles = self.evenParticles[:iv] + self.evenParticles[iv+1:]
        self.setInfo()

    def copy(self):
        """
        Generate an independent copy of self.        
        Faster than deepcopy.
        
        :returns: Branch object
        """

        #Allows for derived classes (like inclusive classes)
        newbranch = self.__class__()
        newbranch.evenParticles = self.evenParticles[:]
        newbranch.oddParticles = self.oddParticles[:]
        self.setInfo()
        newbranch.vertnumb = self.vertnumb
        newbranch.vertparts = self.vertparts[:]
        if not self.maxWeight is None:
            newbranch.maxWeight = self.maxWeight.copy()
        return newbranch


    def getLength(self):
        """
        Returns the branch length (number of odd particles).
        
        :returns: length of branch (number of odd particles)
        """
        
        return len(self.oddParticles)
    
    
    def combineWith(self,other):
        """
        Combines itself with the other branch.
        Should only be used if both branches are considered equal.
        The BSM particles appearing in both branches are combined
        into MultiParticle objects.
        
        :parameter other: branch (Branch Object)        
        """

        if self != other:
            raise SModelSError("Asked to combine distinct branches")
        
        for iptc,bsm in enumerate(other.oddParticles):
            if bsm is self.oddParticles[iptc]:
                continue #Particles are the same (do nothing)

            #Else create a particle list with particles from both
            if not isinstance(self.oddParticles[iptc],MultiParticle):
                bsmList = MultiParticle(label = 'BSM (combined)', particles=[self.oddParticles[iptc]])
            else:
                bsmList = self.oddParticles[iptc]
                
            #Now combine even particles from both branches:
            if not isinstance(bsm,MultiParticle):
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
            
            self.oddParticles[iptc] = bsmList
                        


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
        
        newBranch.oddParticles.append(oddParticles[0])
        evenParticles = sorted(evenParticles, key=lambda x: x.label.lower())
        newBranch.evenParticles.append(evenParticles)
        
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
        
        if not self.oddParticles or not self.oddParticles[-1].decays: 
            return False
        if self.oddParticles[-1].isStable():
            return False

        newBranches = []
        for decay in self.oddParticles[-1].decays:
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
            if any(x is None for x in inbranch.oddParticles[-1].decays):            
                stableBranches.append(inbranch)
            
            # Add all possible decays of the R-odd daughter to the original
            # branch (if any)
            newBranches = inbranch.decayDaughter()
            if newBranches:
                # New branches were generated, add them for next iteration
                newBranchList += [br for br in newBranches if br.maxWeight > sigcut]
            elif inbranch.maxWeight > sigcut:
                stableBranches.append(inbranch)

        # Use new unstable branches (if any) for next iteration step
        unstableBranches = newBranchList           
    
    #Sort list by initial branch pdg:        
    finalBranchList = sorted(stableBranches, key=lambda branch: branch.oddParticles[0].pdg)      

    return finalBranchList


class InclusiveBranch(Branch):
    """
    An inclusive branch class. It will return True when compared to any other branch object
    with the same final state.
    """
    
    def __init__(self,finalState=None):
        Branch.__init__(self)
        self.masses = InclusiveList()
        self.evenParticles =  InclusiveList()
        if finalState:
            bsmParticle = finalStates.getParticlesWith(label=finalState)
            if not bsmParticle:
                raise SModelSError("Final state %s has not been defined in finalStateParticles.py" %finalState)
            if len(bsmParticle) != 1:
                raise SModelSError("Ambiguos defintion of label %s in finalStates" %bsmParticle[0].label)          
            self.oddParticles = [bsmParticle[0]]
        else:
            self.oddParticles = [Particle(Z2parity='odd')]
        self.vertnumb = InclusiveValue()
        self.vertparts = InclusiveList()
        
    def __cmp__(self,other):
        """
        Always returns true. The only exception is if a final state particle has been
        defined. In this case, will include the final state in the comparison.        
        The comparison is made based on vertnumb, vertparts, evenParticles, masses of BSM particles
        and the last BSM particle appearing in the cascade decay.
        OBS: The particles inside each vertex MUST BE sorted (see branch.sortParticles())         
        :param other:  branch to be compared (Branch object)
        :return: -1 if self < other, 0 if self == other, +1, if self > other.
        """


        if not isinstance(other,(Branch,InclusiveBranch)):
            return -1
        
        #If BSM particles are identical, avoid further checks                
        if self.oddParticles and other.oddParticles:
            #Compare final BSM state by Z2parity and quantum numbers:
            comp = self.oddParticles[-1].cmpProperties(other.oddParticles[-1],
                                                       properties=['Z2parity','colordim','eCharge'])
            if comp:
                return comp 
    
        return 0  #Branches are equal    
        
        
       
    def __str__(self):
        return '[*]'
    
    
    def getMasses(self):
        return InclusiveList()
        
    def getInfo(self):
        """
        Get branch topology info (inclusive list and int).
        
        :returns: dictionary containing vertices and number of final states information  
        """

        vertnumb = InclusiveValue()
        vertparts = InclusiveList()
        
        return {"vertnumb" : vertnumb, "vertparts" : vertparts}
        
    def decayDaughter(self):
        """
        Always return False.
        """
        
        return False
        


