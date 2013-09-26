#!/usr/bin/env python

"""
.. module:: Theory.Element
     :synopsis: missing
        
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""
from ParticleNames import PtcDic, Reven, simParticles
import copy

class Branch(object):
    """ A branch """

    def __init__( self, S=None ):
        """ A branch-element can be constructed from a string S (e.g. ('[b,b],[W]')"""
        self.masses = []
        self.particles = []
        self.momID = None
        self.daughterID = None
        self.maxWeight = None
        if type(S)==type(""):
            st = S.replace(" ","")
            while "[" in st or "]" in st:
                ptcs = st[st.find("[")+1:st.find("],[")].split(",")
                for ptc in ptcs:
                    if not ptc in Reven.values() and not PtcDic.has_key(ptc):
                        print "[Branch] unknown particle:",ptc
                        return
                spcts=str(ptcs).replace("'","").replace(" ","")
                st=st.replace(spcts,"",1)
                self.particles.append ( ptcs )
                
    def __str__ ( self ):
        """ the canonical SModels description of the Branch. """
        st = str(self.particles).replace("'","")
        st = st.replace(" ","")
        return st
                
    def __eq__ ( self, other ):
        return self.isEqual ( other )
    
    def __ne__ ( self, other ):
        return not self.isEqual ( other )


    def isEqual ( self, other, useDict=True):
        if type (other) != type(self): return False
        if not simParticles(self.particles,other.particles,useDict): return False
        if self.masses != other.masses: return False
        return True

    
    def addDecay(self,BR,Massdic):
        """ Generates a new branch adding a 1-step cascade decay described by the BR object, with particle masses given by Massdic"""
        import copy
        import ParticleNames
        newBranch = copy.deepcopy(self)
        newparticles = []
        newmass = []

        for partID in BR.ids:
            if partID in ParticleNames.Reven:
                newparticles.append(ParticleNames.Reven[partID])
            elif partID in ParticleNames.Rodd:
                newmass.append(Massdic[partID])
                newBranch.daughterID = partID
            else:
                print '[addDecay] unknown particle:',partID
                return False

        if len(newmass) != 1:
            print '[addDecay] R-parity violating decay:',BR
            return False

        newBranch.particles.append(newparticles)
        newBranch.masses.append(newmass[0])
        newBranch.maxWeight = self.maxWeight*BR.br

        return newBranch

    def addDecays(self,BRs,Massdic):
        """ Generates a list of branches adding all 1-step cascade decays described by the BRs object, with particle masses given by Massdic"""

        if len(BRs) == 0: return [copy.deepcopy(self)]
        newBranches = []
        for br in BRs: newBranches.append(self.addDecay(br,Massdic))
        return newBranches
