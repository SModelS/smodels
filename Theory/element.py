#!/usr/bin/env python

"""
.. module:: Theory.Element
   :synopsis: missing
    
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
    
"""
from ParticleNames import PtcDic, Reven, simParticles
from Tools.PhysicsUnits import rmvunit


class Element(object):
    """
    Abstract class. DO NOT INSTANTIATE!
    """


    def __init__(self, particleString):
        """
        Constructor
        """
        self.ParticleStr = particleString
        

    def getParticleList(self):
        """
        Converts the particle string in self.ParticleStr to a list of particles.
        """
        if self.ParticleStr == "": return None
        particles = [[],[]]
        S = self.ParticleStr
        st = S.replace(" ","").replace("'","")
        st = st[st.find("[[["):st.find("]]]")+3]
        branches=[st[2:st.find("]],[[")+1],st[st.find("]],[[")+4:st.find("]]]")+1]]
        for ib,branch in enumerate(branches):
            st = branch
            while "[" in st or "]" in st:
                ptcs = st[st.find("[")+1:st.find("],[")].split(",")
                for ptc in ptcs:
                    if not ptc in Reven.values() and not PtcDic.has_key(ptc):
                        print "[AElement] unknown particle:",ptc
                        return False
                spcts=str(ptcs).replace("'","").replace(" ","")
                st=st.replace(spcts,"",1)
                particles[ib].append(ptcs)
    
        return particles
    
    
    def getEinfo( self ):
        """
        Get global topology info from particle string.
        """
        vertnumb = []
        vertparts = []
        ptcs = self.getParticleList()
        for branch in ptcs:
            vertnumb.append(len(branch))
            vertparts.append([len(v) for v in branch])
            if len(vertparts[len(vertparts)-1]) == vertnumb[len(vertnumb)-1]-1:
                vertparts[len(vertparts)-1].append(0)  #Append 0 for stable LSP
        return {"vertnumb" : vertnumb, "vertparts" : vertparts}
      
    
class AElement(Element):
    """
    An Analysis Element, contains a string with the particle list, a dictionary with mass arrays and the respective weights and the\
    analysis-dependent weight format
    """
    
    def __init__(self, particleString=""):
        """
        Constructor
        """
        self.MassWeightList = []
        super(particleString)
        
        
class CElement(Element):
    """
    A Cluster Element, contains a simple string with the particle list and its weight
    """
    
    def __init__(self, particleString="", weight={}):
        self.Weight = weight
        super(particleString)
        
        
class EElement(object):
    """
    An Event Element, contains of several branches and weight information.
    """
    
    def __init__(self, S=None):
        """
        If S != None, an Element is created from a string description.
        """
        self.B = []
        self.weight = []
        if S:
            st = S.replace(" ","").replace("'","")
            st = st[st.find("[[["):st.find("]]]")+3]
            b1=st[2:st.find("]],[[")+1]
            b2=st[st.find("]],[[")+4:st.find("]]]")+1]
            import copy
            self.B.append ( copy.deepcopy( BElement ( b1 ) ) )
            self.B.append ( copy.deepcopy ( BElement ( b2 ) ) )
            
            
    def getEinfo( self ):
        """
        Get global topology info from element structure.
        """
        vertnumb = []
        vertparts = []
        for el in self.B:
            vertnumb.append(len(el.masses))
            vertparts.append([len(x) for x in el.particles])
            if len(vertparts[len(vertparts)-1]) == vertnumb[len(vertnumb)-1]-1:
                vertparts[len(vertparts)-1].append(0)  #Append 0 for stable LSP
            return {"vertnumb" : vertnumb, "vertparts" : vertparts}
    
    def allParticles ( self ):
        """
        Returns all particles from all branches.
        """
        ret=[]
        for b in self.B:
            ret.append ( b.particles )
            return ret
    
    def isSimilar ( ElA, ElB,order=True,igmass=False ):
        """
        Compare two EElements
        If particles are similar and all masses equal, returns True,
        otherwise returns False
        If order = False, test both branch orderings (for an element doublet only)
        If igmass = True, only compare particles
        """
        if type(ElA) != type(ElB): return False
        El1 = ElA.B
        El2 = ElB.B
        if len(El1) == 2 and not order:
            ptcsA = [El2[0].particles,El2[1].particles]
            massA = [El2[0].masses,El2[1].masses]
            ptcs = [El1[0].particles,El1[1].particles]
            mass = [El1[0].masses,El1[1].masses]
            ptcs_b = [El1[1].particles,El1[0].particles]
            mass_b = [El1[1].masses,El1[0].masses]
            if igmass:
                mass = massA
                mass_b = massA
            if simParticles(ptcsA,ptcs) and mass == massA:
                return True
            elif simParticles(ptcsA,ptcs_b) and mass_b == massA:
                return True
            else:
                return False
        else:
            for i in range(len(El1)):
                if not simParticles(El1[i].particles,El2[i].particles): return False
                if not igmass and El1[i].masses != El2[i].masses: return False
        return True
    
    def isEqual ( ElA, ElB,order=True):
        """
        Compare two EElements
        If all masses and particles are equal, returns True,
        otherwise returns False
        If order = False, test both branch orderings (for an element doublet only)
        """
        if type(ElA) != type(ElB): return False
        El1 = ElA.B
        El2 = ElB.B
        if len(El1) == 2 and not order:
            ptcsA = [El2[0].particles,El2[1].particles]
            massA = [El2[0].masses,El2[1].masses]
            ptcs = [El1[0].particles,El1[1].particles]
            mass = [El1[0].masses,El1[1].masses]
            ptcs_b = [El1[1].particles,El1[0].particles]
            mass_b = [El1[1].masses,El1[0].masses]
            if simParticles(ptcsA,ptcs,useDict=False) and mass == massA:
                return True
            elif simParticles(ptcsA,ptcs_b,useDict=False) and mass_b == massA:
                return True
            else:
                return False
        else:
            for i in range(len(El1)):
                if not simParticles(El1[i].particles,El2[i].particles,useDict=False): return False
                if El1[i].masses != El2[i].masses: return False
                
        return True
    
    def __eq__ ( self, other ):
        return self.isEqual ( other )
    
    def __str__ ( self ):
        """
        Returns the canonical name of the element, e.g. [[jet],[jet]].
        """
        ret="["
        for i in self.B:
            ret+=str(i)+","
        if len(ret)>1:
            ret=ret[:-1]
        ret+="]"
        return ret
    
    def describe ( self ):
        """
        Returns a lengthy description of the event element.
        """
        ret="Branch #1={{"+str(self.B[0])+"}}, Branch #2={{"+str(self.B[1])+"}}"
        return ret


class BElement(object):
    """
    A branch-element.
    """
    
    def __init__( self, S=None ):
        """
        Constructor.
        A branch-element can be constructed from a string S (e.g. ('[b,b],[W]').
        """
        self.masses = []
        self.particles = []
        self.momID = 0
        if type(S)==type(""):
            st = S.replace(" ","")
            while "[" in st or "]" in st:
                ptcs = st[st.find("[")+1:st.find("],[")].split(",")
                for ptc in ptcs:
                    if not ptc in Reven.values() and not PtcDic.has_key(ptc):
                        print "[BElement] unknown particle:",ptc
                        return
                spcts=str(ptcs).replace("'","").replace(" ","")
                st=st.replace(spcts,"",1)
                self.particles.append ( ptcs )


    def isEqual ( ElA, ElB, order=True ):
        if not simParticles(ElA.particles,ElB.particles,useDict=False): return False
        if ElA.masses != ElB.masses: return False
        return True


    def __eq__ ( self, other ):
        return self.isEqual ( other )


    def isSimilar ( self, elB, order=True, igmass=False ):
        """
        Compare elB with self.
        If particles are similar and all masses equal, returns True,
        otherwise returns False.
        If order = False, test both branch orderings (for an element doublet only)
        If igmass = True, only compare particles
        """
        if type (elB) != type(self): return False
        if not simParticles(self.particles,elB.particles): return False
        if not igmass and self.masses != elB.masses: return False
        return True

    def __str__ ( self ):
        """
        The canonical SModels description of the BElement.
        """
        st = str(self.particles).replace("'","")
        st = st.replace(" ","")
        return st


    def describe ( self ):
        """
        A lengthy description of the BElement.
        """
        ret="particles=%s masses=%s" % \
                ( self.particles, [ rmvunit(x,"GeV") for x in self.masses ] )
        return ret

        
class MassWeight(object):
    """
    A simple object holding a mass array and its respective weight.
    """
    
    def __init__(self, mass=[], weight={}):
        """
        Constructor
        """
        self.mass = mass
        self.weight = weight
        

