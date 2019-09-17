"""
.. module:: element
   :synopsis: Module holding the Element class and its methods.
    
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
    
"""

from smodels.theory.auxiliaryFunctions import stringToTree, getTopologyName, getNodeLevelDict, getTreeRoot
from smodels.theory import crossSection
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger
from smodels.theory.particle import Particle
import networkx as nx

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
        self.tree = nx.DiGraph()
        self.weight = crossSection.XSectionList() # gives the weight for all decays promptly
        self.decayLabels = []
        self.motherElements = [("original", self)]
        self.elID = 0
        self.covered = False
        self.tested = False
                
        if info:
            if isinstance(info,str):
                self.tree = stringToTree(info,finalState)
            elif isinstance(info,nx.DiGraph):
                self.tree = info.copy() #Makes a shallow copy of the original tree
            else:
                raise SModelSError("Can not create element from input type %s" %type(info))
        
        self.setEinfo()
        
        
    def setEinfo(self):
        """
        Compute and store the canonical name for the tree
        topology. The canonical name can be used to compare and sort topologies.
        The name is stored in self.tree.topologyName
        """
        
        canonName = getTopologyName(self.tree)
        self.tree.graph['topologyName'] = canonName
    
    def __cmp__(self,other):
        """
        Compares the element with other.
        Uses the topology name (Tree canonincal name) to identify isomorphic topologies (trees).
        If trees are not isomorphic, compare topologies using the topology name.
        Else  check if nodes (particles) are equal. If particles match return 0,
        else compare lists of particles in each tree level. 

        :param other:  element to be compared (Element object)
        
        :return: -1 if self < other, 0 if self == other, +1, if self > other.
        """

        if not isinstance(other,Element):
            return -1
        
        #make sure the topology names have been computed:
        if not 'topologyName' in self.tree.graph:
            self.setEInfo()
        if not 'topologyName' in other.tree.graph:
            other.setEInfo()
            
        tnameA = self.tree.graph['topologyName']
        tnameB = other.tree.graph['topologyName']
        if tnameA != tnameB:
            return (tnameA > tnameB) - (tnameA < tnameB)
        
        #Check if the elements are isomorphic including particle comparison
        if nx.isomorphism.is_isomorphic(self.tree,other.tree,node_match=self.equalNodes):
            return 0
        else:
            #Compare particles in each level of the tree:
            #First build dictionary with nodes in each level:
            levelNodesA = getNodeLevelDict(self.tree)                                      
            levelNodesB = getNodeLevelDict(other.tree)
            for nLevel in levelNodesA:
                pListA = [self.tree.nodes[n]['particle'] for n in levelNodesA[nLevel]]
                pListB = [other.tree.nodes[n]['particle'] for n in levelNodesB[nLevel]]
                pListA = sorted(pListA)
                pListB = sorted(pListB)
                pComp = (pListA > pListB) - (pListA < pListB)
                if pComp != 0:
                    return pComp
            
            return 0 #Should never reach this point!
    
    def equalNodes(self,n1,n2):
        """
        Return True if particles in nodes are equal, False
        otherwise.

        :param n1: Dictionary with the node properties
        :param n2: Dictionary with the node properties
        
        :return: True if n1 == n2, False otherwise.
        """
        
        pA = n1['particle']
        pB = n2['particle']
        
        return pA == pB

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
        
        T = self.tree
        elStr = ""
        root = getTreeRoot(T)
        for mom,daughters in nx.bfs_successors(T,root):
            elStr += '%s->(%s),' %(T.nodes[mom]['particle'].label,
                                   ','.join(T.nodes[n]['particle'].label for n in daughters))
        return elStr[:-1]
    
    def __repr__(self):
        
        return self.__str__()

    def toStr(self):
        """
        Returns a string with the element represented in bracket notation,
        including the final states, e.g. [[[jet]],[[jet]] (MET,MET)
        """
        
        elStr = str(self)
        
        return elStr
    
    def drawTree(self,outputFile=None,show=True,
                 oddColor = 'lightcoral',evenColor = 'skyblue',
                 pvColor = 'darkgray',genericColor= 'violet',
                 nodeScale = 4):
        """
        Draws element Tree using matplotlib.
        If outputFile is defined, it will save plot to this file.
        """
        
        import matplotlib.pyplot as plt
        
        T = self.tree
        
        labels = dict([[n,str(T.nodes[n]['particle'])] 
                       for n in T.nodes()])
        for key in labels:
            if labels[key] == 'anyOdd':
                labels[key] = 'BSM'
        node_size = []
        node_color = []
        for n in T.nodes():
            node_size.append(nodeScale*100*len(labels[n]))
            if 'pv' == labels[n].lower():
                node_color.append(pvColor)
            elif hasattr(T.nodes[n]['particle'],'Z2parity'):
                if T.nodes[n]['particle'].Z2parity == 'odd':
                    node_color.append(oddColor)
                else:
                    node_color.append(evenColor)
            else:
                node_color.append(genericColor)
        pos = nx.drawing.nx_agraph.graphviz_layout(T, prog='dot')
        nx.draw(T,pos,
                with_labels=True,
                arrows=True,
                labels=labels,
                node_size=node_size,
                node_color=node_color)
        
        if outputFile:            
            plt.savefig(outputFile)
        if show:
            plt.show()


    def copy(self):
        """
        Create a copy of self.        
        Faster than deepcopy.     
        
        :returns: copy of element (Element object)   
        """

        #Allows for derived classes (like inclusive classes)
        newel = self.__class__()
        newel.tree = self.tree.copy()
        newel.weight = self.weight.copy()
        newel.motherElements = self.motherElements[:]
        newel.elID = self.elID
        return newel

    def getParticles(self):
        """
        Get the array of even particle objects in the element.
        
        :returns: list of Particle objects                
        """

        g = self.tree
        particles = []
        
        for iv,v in enumerate(g.vertices()):
            if v.out_degree():
                continue
            if g.vp.particle[iv].Z2parity != 'even':
                continue        
            particles.append(g.vp.particle[iv])

        return particles
    
    def getFinalStates(self):
        """
        Get the array of final state (last BSM particle) particle objects in the element.
        
        :returns: list of Particle objects
        """
        
        g = self.tree
        particles = []
        
        for iv,v in enumerate(g.vertices()):
            if v.out_degree():
                continue
            if g.vp.particle[iv].Z2parity != 'odd':
                continue        
            particles.append(g.vp.particle[iv])

        return particles


    def getMasses(self):
        """
        Get the array of BSM masses in the element.    
        
        :returns: list of masses (mass array)            
        """        

        massarray = [ptc.mass for ptc in self.getBSMparticles()]

        return massarray
        
    def getBSMparticles(self):
        """
        Get the list of BSM particles appearing the cascade decay,
        including the last (stable) one.

        :returns: list of Particle or ParticleList objects
        """    
        

        g = self.tree
        BSMparticles = [g.vp.particle[v] for v in g.vertices() if g.vp.particle[v].Z2parity == 'odd']

        return BSMparticles

    def getPIDs(self):
        """
        Get the list of IDs (PDGs of the intermediate states appearing the cascade decay), i.e.
        [  [[pdg1,pdg2,...],[pdg3,pdg4,...]] ].
        
        :returns: list of PDG ids
        """

        BSMpids = [particle.pdg for particle in self.getBSMparticles()]

        return BSMpids

    def getDaughters(self):
        """
        Get the list of daughter (last/stable BSM particle in the decay) PDGs. 
        Can be a nested list, if the element combines several daughters:
        [ [pdgDAUG1,pdgDAUG2],  [pdgDAUG1',pdgDAUG2']] 

        
        :returns: list of PDG ids
        """
        
        daughterPIDs = [particle.pdg for particle in self.getFinalStates()]
        
        return daughterPIDs


    def getEinfo(self):
        """
        Get element topology info from branch topology info.
        
        :returns: dictionary containing vertices and number of final states information  
        """
                
        return {"vertices" : self.tree.num_vertices(), "edgesout" : self.tree.degree_property_map('out').get_array()}


    def _getLength(self):
        """
        Get the maximum of the two branch lengths.    
        
        :returns: maximum length of the element branches (int)    
        """
        return len(self.getBSMparticles())

    
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
     
            removeVertices = []
            for i,mom in enumerate(branch.oddParticles[:-1]):
                massDiff = mom.mass - branch.oddParticles[i+1].mass
                #Get vertices which have deltaM < minmassgap and the mother is prompt:
                if massDiff < minmassgap and mom.isPrompt():
                    removeVertices.append(i)
            #Remove vertices till they exist:
            while removeVertices:
                newelement.removeVertex(ibr,removeVertices[0])
                branch = newelement.branches[ibr]
                removeVertices = []
                for i,mom in enumerate(branch.oddParticles[:-1]):
                    massDiff = mom.mass - branch.oddParticles[i+1].mass
                    #Get vertices which have deltaM < minmassgap and the mother is prompt:
                    if massDiff < minmassgap and mom.isPrompt():
                        removeVertices.append(i)
                
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
            neutralDecay = neutralSM and branch.oddParticles[-1].isMET()
            #Check if the mother can be considered MET:
            neutralBSM = branch.oddParticles[-2].isMET()
            if neutralBSM and neutralDecay:
                removeLastVertex = True
            else:
                removeLastVertex = False
                
            while len(branch.oddParticles) > 1 and removeLastVertex:
                bsmMom = branch.oddParticles[-2]
                effectiveDaughter = Particle(label='inv', mass = bsmMom.mass,
                                             eCharge = 0, colordim = 1,
                                             totalwidth = branch.oddParticles[-1].totalwidth,
                                             Z2parity = bsmMom.Z2parity)
                branch.removeVertex(len(branch.oddParticles)-2)
                #For invisible compression, keep an effective mother which corresponds to the invisible
                #daughter, but with the mass of the parent.
                branch.oddParticles[-1] = effectiveDaughter
                #Re-check if the last decay should be removed:
                if not branch.evenParticles:
                    continue
                neutralSM = all(ptc.isMET() for ptc in branch.evenParticles[-1])
                neutralDecay = neutralSM and branch.oddParticles[-1].isMET()
                neutralBSM = branch.oddParticles[-2].isMET()
                if neutralBSM and neutralDecay:
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

