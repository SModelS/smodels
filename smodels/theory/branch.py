"""
.. module:: branch
   :synopsis: Module holding the branch class and methods.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>

"""

from smodels.theory.auxiliaryFunctions import elementsInStr
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger
from smodels.theory.particle import MultiParticle,ParticleList
from smodels.tools.inclusiveObjects import InclusiveValue,InclusiveList



class Branch(object):
    """
    An instance of this class represents a branch.
    A branch-element can be constructed from a string (e.g., ('[b,b],[W]').
    """

    def __init__(self, info=None, finalState=None, intermediateState=None, model=None):
        """
        Initializes the branch. If info is defined, tries to generate
        the branch using it.

        :parameter info: string describing the branch in bracket notation
                         (e.g. [[e+],[jet]])

        :parameter finalState: final state label string for the branch
                         (e.g. 'MET' or 'HSCP')
        :parameter intermediateState: list containing intermediate state labels
                                      (e.g. ['gluino','C1+'])
        :parameter model: The model (Model object) to be used when converting particle labels to
                          particle objects (only used if info, finalState or intermediateState != None).

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
            if not model or not hasattr(model,'getParticlesWith'):
                raise SModelSError("A Model object has to be defined when creating branches from strings.")
            else:
                branch = branch[0]
                vertices = elementsInStr(branch[1:-1])
                for vertex in vertices:
                    particleNames = vertex[1:-1].split(',')
                    ptcs = []
                    for pname in particleNames:
                        smParticle = model.getParticlesWith(label=pname)
                        if not smParticle:
                            raise SModelSError("Final state SM particle ``%s'' has not been defined in %s" %(pname,model))
                        elif len(smParticle) != 1:
                            raise SModelSError("Ambiguous definition of label ``%s'' in %s" %(smParticle[0].label,model))
                        else:
                            ptcs.append(smParticle[0])
                    vertexParticles = ParticleList(ptcs)
                    self.evenParticles.append(vertexParticles)

            self.vertnumb = len(self.evenParticles)
            self.vertparts = [len(v) for v in self.evenParticles]

            #Get labels of intermediate states (default is [anyOdd,anyOdd,...,MET])
            if intermediateState:
                if not isinstance(intermediateState,list):
                    raise SModelSError("Intermediate state (``%s'') should be a list)" %intermediateState)
                bsmLabels = intermediateState[:]
            else:
                bsmLabels = ['anyOdd']*self.vertnumb
            if finalState:
                bsmLabels.append(finalState)
            else:
                bsmLabels.append('MET')
            if len(bsmLabels) != self.vertnumb+1:
                raise SModelSError("Number of intermediate states (``%s'') is not consistent)" %intermediateState)
            for bsmLabel in bsmLabels:
                bsmParticle = model.getParticlesWith(label=bsmLabel)
                if not bsmParticle:
                    raise SModelSError("BSM particle ``%s'' has not been defined in databaseParticles.py" %bsmLabel)
                elif len(bsmParticle) != 1:
                    raise SModelSError("Ambiguous definition of label ``%s'' in databaseParticles.py" %bsmLabel)
                else:
                    self.oddParticles.append(bsmParticle[0])

    def __str__(self):
        """
        Create the branch bracket notation string, e.g. [[e+],[jet]].

        :returns: string representation of the branch (in bracket notation)
        """

        sortedParticles = [sorted(vertex, key = lambda ptc: str(ptc))
                           for vertex in self.evenParticles]
        st = str(sortedParticles).replace("'", "")
        st = st.replace(" ", "")
        return st

    def __repr__(self):

        return self.__str__()

    def __cmp__(self,other):
        """
        Compares the branch with other.
        The comparison is made based on vertnumb, vertparts, oddParticles and evenParticles.
        The comparison allows for any ordering of the evenParticles in the vertex.
        It relies on the particle comparison, which allows for the comparison of Particles and MultiParticles.
        Only the properties which are defined for both particles are compared.
        :param other:  branch to be compared (Branch object)
        :return: -1 if self < other, 0 if self == other, +1, if self > other.
        """

        if isinstance(other,InclusiveBranch):
            return -1*other.__cmp__(self)

        if self.vertnumb != other.vertnumb:
            comp = self.vertnumb > other.vertnumb
            if comp: return 1
            else: return -1
        elif self.vertparts != other.vertparts:
            comp = self.vertparts > other.vertparts
            if comp: return 1
            else: return -1

        #Compare BSM states:
        if self.oddParticles != other.oddParticles:
            comp = self.oddParticles > other.oddParticles
            if comp:
                return 1
            else:
                return -1

        #Compare even final states (ParticleLists)
        #The comparison of ParticleList objects is made irrespective of ordering
        for iv,vertex in enumerate(self.evenParticles):
            if vertex != other.evenParticles[iv]:
                comp = vertex > other.evenParticles[iv]
                if comp:
                    return 1
                else:
                    return -1

        return 0  #Branches are equal

    def __lt__( self, b2 ):
        return self.__cmp__(b2) == -1

    def __eq__( self, b2 ):
        return self.__cmp__(b2) == 0

    def __ne__( self, b2 ):
        return not self.__cmp__(b2) == 0

    def __getattr__(self, attr):
        """
        If the attribute has not been defined for the element
        try to fetch it from its branches.
        :param attr: Attribute name

        :return: Attribute value
        """

        #If calling another special method, return default (required for pickling)
        if attr.startswith('__') and attr.endswith('__'):
            return object.__getattr__(attr)

        try:
            val = [getattr(ptc,attr) for ptc in self.oddParticles]
            return val
        except AttributeError:
            raise AttributeError("Element nor branch has attribute %s" %attr)

    def __add__(self,other):
        """
        Adds two branches. Should only be used if the branches
        have the same topologies. The odd and even particles are combined.
        """

        if self.getInfo() != other.getInfo():
            raise SModelSError("Can not add branches with distinct topologies")

        newBranch = self.__class__()
        #Combine odd particles
        for iptc,ptc in enumerate(self.oddParticles):
            newBranch.oddParticles.append(ptc+other.oddParticles[iptc])

        #Combine even particles (if they are the same nothing changes)
        for iv,vertex in enumerate(self.evenParticles):
            vertexParticles = []
            for iptc,ptc in enumerate(vertex):
                vertexParticles.append(ptc + other.evenParticles[iv][iptc])

            vertexParticles = ParticleList(vertexParticles)
            newBranch.evenParticles.append(vertexParticles)

        if not self.maxWeight is None and not other.maxWeight is None:
            newBranch.maxWeight = self.maxWeight + other.maxWeight

        return newBranch

    def __radd__(self,other):
        """
        Adds two elements. Only elements with the same
        topology can be combined.
        """

        return self.__add__(other)

    def __iadd__(self,other):
        """
        Combine two branches. Should only be used if the elements
        have the same topologies. The branches
        odd and even particles are combined.
        """

        if self.getInfo() != other.getInfo():
            raise SModelSError("Can not add branches with distinct topologies")

        #Combine odd particles
        for iptc,ptc in enumerate(other.oddParticles):
            self.oddParticles[iptc] += ptc

        #Combine even particles (if they are the same nothing changes)
        for iv,vertex in enumerate(self.evenParticles):
            for iptc,ptc in enumerate(vertex):
                self.evenParticles[iv][iptc] += other.evenParticles[iv][iptc]
        if not self.maxWeight is None and not other.maxWeight is None:
            self.maxWeight += other.maxWeight

        return self

    def getAverage(self,attr):
        """
        Get the average value for a given attribute appearing in
        the odd particles of branch.
        """

        try:
            vals = []
            for ptc in self.oddParticles:
                v = getattr(ptc,attr)
                if isinstance(ptc,MultiParticle) and isinstance(v,list):
                    avg = v[0]
                    for x in v[1:]:
                        avg += x
                    avg =avg/len(v)
                else:
                    avg = v
                vals.append(avg)
        except (AttributeError,ZeroDivisionError):
            raise SModelSError("Could not compute average for %s" %attr)

        return vals

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
            newbranch.maxWeight = self.maxWeight
        return newbranch

    def getLength(self):
        """
        Returns the branch length (number of odd particles).

        :returns: length of branch (number of odd particles)
        """

        return len(self.oddParticles)

    def _addDecay(self, decay):
        """
        Generate a new branch adding a 1-step cascade decay
        This is described by the br object, with particle masses given by BSMList.

        :parameter decay: Decay object (see pyslha). Contains information about the decay.
        :returns: extended branch (Branch object). False if there was an error.
        """

        newBranch = self.copy()
        oddParticles = decay.oddParticles
        evenParticles = decay.evenParticles

        if len(oddParticles) != 1:
            logger.warning("Decay %s does not preserve Z2 and will be ignored" %str(decay))
            return False

        newBranch.oddParticles.append(oddParticles[0])
        newBranch.evenParticles.append(evenParticles)

        if not self.maxWeight is None:
            newBranch.maxWeight =  self.maxWeight*decay.br

        newBranch.setInfo()
        return newBranch


    def decayDaughter(self):
        """
        Generate a list of all new branches generated by the 1-step cascade
        decay of the current branch daughter.
        :returns: list of extended branches (Branch objects). Empty list if daughter is stable or if daughterID was not defined.
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


def decayBranches(branchList, sigcut=0.):
    """
    Decay all branches from branchList until all unstable intermediate states have decayed.

    :parameter branchList: list of Branch() objects containing the initial mothers
    :parameter sigcut: minimum sigma*BR (in fb) to be generated, by default sigcut = 0.
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
            if sigcut > 0. and inbranch.maxWeight < sigcut:
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
                newBranchList += [br for br in newBranches
                                  if br.maxWeight > sigcut]
            elif inbranch.maxWeight > sigcut:
                stableBranches.append(inbranch)

        # Use new unstable branches (if any) for next iteration step
        unstableBranches = newBranchList

    #Sort list by initial branch pdg:
    finalBranchList = sorted(stableBranches, key=lambda branch: branch.oddParticles[0].pdg)

    return finalBranchList


class InclusiveBranch(Branch):
    """
    An inclusive branch class. It will return True when compared to any other branch object with the same final state. If intermediateState is defined, it will
    use all the odd (BSM) particles for comparison while neglecting even particles.
    """

    def __init__(self,finalState=None,intermediateState=None,model=None):
        """
        :parameter info: string describing the branch in bracket notation
                         (e.g. [[e+],[jet]])

        :parameter finalState: final state label string for the branch
                         (e.g. 'MET' or 'HSCP')
        :parameter intermediateState: list containing intermediate state labels
                                      (e.g. ['gluino','C1+'])
        :parameter model: The model (Model object) to be used when converting particle labels to particle objects (only used if info, finalState or intermediateState != None).
        """
        Branch.__init__(self)
        self.mass = InclusiveList()
        self.totalwidth = InclusiveList()
        self.evenParticles =  InclusiveList()
        self.oddParticles = []
        #Get labels of intermediate states
        if intermediateState:
            if not isinstance(intermediateState,list):
                raise SModelSError("Intermediate state (``%s'') should be a list)" %intermediateState)
            bsmLabels = intermediateState[:]
        else:
            bsmLabels = []
        if finalState:
            bsmLabels.append(finalState)
        else:
            bsmLabels.append('anyOdd')
        for bsmLabel in bsmLabels:
            bsmParticle = model.getParticlesWith(label=bsmLabel)
            if not bsmParticle:
                raise SModelSError("BSM particle ``%s'' has not been defined in model %s" %(bsmLabel,model))
            elif len(bsmParticle) != 1:
                raise SModelSError("Ambiguous definition of label ``%s'' in model %s" %(bsmLabel,model))
            else:
                self.oddParticles.append(bsmParticle[0])

        #If intintermediateState is defined, the number of vertexParticles
        #(odd particles) must be used for comparison. Otherwise accept any number.
        if intermediateState:
            self.vertnumb = len(intermediateState)
        else:
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

        #First compare number of vertices. If no intermediate states were
        #defined, should always pass
        if self.vertnumb != other.vertnumb:
            comp = self.vertnumb > other.vertnumb
            if comp: return 1
            else: return -1

        #Compare BSM states in reverse order, so the final state is always the first to be compared.
        # (if no intermediate states were defined, self.oddParticles contain only the final state)
        for iptc,ptc in enumerate(self.oddParticles[::-1]):
            iother = -1-iptc #reverse index
            if ptc == other.oddParticles[iother]:
                continue
            comp = ptc > other.oddParticles[iother]
            if comp:
                return 1
            else:
                return -1

        return 0  #Branches are equal

    def __str__(self):
        return '[*]'

    def getInfo(self):
        """
        Get branch topology info (inclusive list and int).

        :returns: dictionary containing vertices and number of final states information
        """

        vertnumb = self.vertnumb
        vertparts = self.vertparts

        return {"vertnumb" : vertnumb, "vertparts" : vertparts}

    def decayDaughter(self):
        """
        Always return False.
        """

        return False
