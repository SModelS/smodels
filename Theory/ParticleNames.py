#!/usr/bin/env python

"""
.. module:: ParticleNames
        :synopsis: methods for getting particle names from pdg ids, and other helpers

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import auxiliaryFunctions
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


#Converts pdg number to particle name according to the dictionaries Rodd
# and Reven
def getName(pdg):
    """ Converts pdg number to particle name 
        according to the dictionaries Rodd and Reven 

        :type pdg: int
        :returns: particle name (e.g. gluino, mu-, ...)
    """
    p=int(pdg)
    if p in Rodd: return Rodd[p]
    if p in Reven:
        return Reven[p]
    else:
        return False


def getPdg(name):
    """ Converts a name to the pdg number
        according to the dictionaries Rodd and Reven

        :type name: string
        :returns: particle pdg, None if name couldnt be resolved
    """
    for (pdg,pname) in Rodd.items():
        if name==pname: return abs(pdg)
    for (pdg,pname) in Reven.items():
        if name==pname: return abs(pdg)
    return None

#class Part(object):
#  
#    def __init__(self,ptc=""):
#        self.label = ptc
#    
#    def __str__(self):
#        return self.label
#  
#    def __eq__(self,other):
#        if self.label == other.label: return True
#        else:
#            if PtcDic.has_key(self.label): ptcs = set(PtcDic[self.label])
#            else: ptcs = set([self.label])
#            if PtcDic.has_key(other.label): ptcs2 = set(PtcDic[other.label])
#            else: ptcs2 = set([other.label])
#            if ptcs.intersection(ptcs2): return True
#            else: return False


def elementsInStr(instring):
    """ Parses instring (or a list of strings) and return a list of elements (in string format) appearing in instring"""
  
    if type(instring) == type('st'): outstr = instring
    elif type(instring) == type([]):
        outstr = ""
        for st in instring:
            if type(st) != type('st'):
                logger.error("Input must be a string or a list of strings")
                return False        
            outstr += st    #Combines list of strings in a single string
  
    elements = []
    outstr = outstr.replace(" ","").replace("'","")
    elStr = ""
    nc = 0
#Parses the string and looks for matching ['s and ]'s, when the matching is complete, store element    
    for c in outstr:
        delta = 0
        if c == '[':  delta = -1
        elif c == ']': delta = 1
        nc += delta
        if nc != 0: elStr += c
        if nc == 0 and delta != 0:
            elements.append(elStr+c)
            elStr = ""
            #Syntax checks:
            ptclist = elements[-1].replace(']',',').replace('[',',').split(',')
            for ptc in ptclist:
                if not ptc: continue
                if not ptc in Reven.values() and not PtcDic.has_key(ptc):
                    logger.error("Unknown particle "+ptc)
                    return False

#Check if there are not unmatched ['s and/or ]'s in the string:        
    if nc != 0:
        logger.error("Wrong input (incomplete elements?) "+instring)
        return False     
      
    return elements


def vertInStr(instring):
    """ Parses instring (or a list of strings) and returns the list of particle vertices appearing in instring"""
  
    if type(instring) == type('st'): outstr = instring
    elif type(instring) == type([]):
        outstr = ""
        for st in instring:
            if type(st) != type('st'):
                logger.error("Input must be a string or a list of strings")
                return False        
            outstr += st    #Combines list of strings in a single string
  
    vertices = []
    outstr = outstr.replace(" ","").replace("'","")
    vertStr = ""
    nc = 0
#Parses the string and looks for matching ['s and ]'s, when the matching is complete, store element    
    for c in outstr:
        delta = 0
        if c == '[':  delta = -1
        elif c == ']': delta = 1
        nc += delta
        if c == '[': vertStr = ""
        if nc != 0 and c != '[' and c != ']': vertStr += c
        if delta > 0 and vertStr:
            vertices.append(vertStr.split(','))            
            #Syntax checks:
            for ptc in vertices[-1]:
                if not ptc: continue
                if not ptc in Reven.values() and not PtcDic.has_key(ptc):
                    logger.error("Unknown particle "+ptc)
                    return False
            vertStr = ""

#Check if there are not unmatched ['s and/or ]'s in the string:        
    if nc != 0:
        logger.error("Wrong input (incomplete elements?) "+instring)
        return False     
      
    return vertices
 
#def simParticles(ptype1,ptype2,useDict=True):
#    """ Compares 2 particle names or 2 nested name arrays. \
#        Allows for dictionary labels
#        (Ex: L = l, l+ = l, l = l-,...) 
#        For the last nested level ignore particle ordering 
#        :param ptype1: first (nested) list of particle names, e.g. ['l','jet']
#        :param ptype2: second (nested) list of particle names
# 
#        :param useDict: use the translation dictionary, i.e. allow e to stand for e+ or e-, l+ to stand for e+ or mu+, etc
# 
#        :returns: boolean
#    """
#        
#    if ptype1 == ptype2: return True
#    if type(ptype1) != type(ptype2): return False
#    if type(ptype1) == type(str()):
#        ptype1v = [ptype1]
#        ptype2v = [ptype2]
#  
#    dims1, dims2 = [],[]
#    ptcs1 = auxiliaryFunctions.flattenList(ptype1v,dims1)
#    ptcs2 = auxiliaryFunctions.flattenList(ptype2v,dims2)
#    if dims1 != dims2: return False
#    if ptcs1 == ptcs2: return True
#    elif not useDict: return False
##Check dictionaries to see if particles are similar:
#    vertices1 = vertInStr(str(ptcs1))
#    vertices2 = vertInStr(str(ptcs2))
#    for iv,v1 in enumerate(vertices1):
#        v2 = vertices2[iv]
#        if v1 == v2: continue
#        
#     
#    return True


def simParticles(ptype1,ptype2,useDict=True):
     """ Compares 2 particle names or 2 nested name arrays. \
             Allows for dictionary labels
             (Ex: L = l, l+ = l, l = l-,...) 
             For the last nested level ignore particle ordering 
             FIXME nesting? 
 
         :param ptype1: first (nested) list of particle names, e.g. ['l','jet']
         :param ptype2: second (nested) list of particle names
 
         :param useDict: use the translation dictionary, i.e. allow e to stand for e+ or e-, l+ to stand for e+ or mu+, etc
 
         :returns: boolean
     """
     import copy
     
     wrongFormat = False
     if type(ptype1) != type(ptype2): return False
 #Check for nested arrays (should be in the standard notation [[[]],[[]]]):    
     if type(ptype1) == type([]):
         if len(ptype1) != len(ptype2): return False
         for ib,br in enumerate(ptype1):
             if type(br) != type(ptype2[ib]) or type(br) != type([]): wrongFormat = True
             if len(ptype1[ib]) != len(ptype2[ib]): return False    #Check number of vertices in branch
             for iv,vt in enumerate(br):
                 if type(vt) != type(ptype2[ib][iv]) or type(vt) != type([]): wrongFormat = True
                 if len(ptype1[ib][iv]) != len(ptype2[ib][iv]): return False #Check number of particles in vertex
                 for ptc in ptype1[ib][iv]+ptype2[ib][iv]:
                     if not ptc in PtcDic.keys()+Reven.values(): wrongFormat = True
     if wrongFormat:     
         print "[ParticleNames.simParticles]: Wrong input format!",ptype1,ptype2
         return False
     
 #Put input in standard notation
     if type(ptype1) == type("str"):
         ptype1v = [[[ptype1]],[[ptype1]]]
         ptype2v = [[[ptype2]],[[ptype2]]]
     else:
         ptype1v = ptype1[:]
         ptype2v = ptype2[:]
 
     for ibr,br in enumerate(ptype1v): #Loop through branches
         for iv,vt in enumerate(br): #Loop over vertices    
 #Check    if lists match, ignoring possible dictionary entries
             pmatch = True
             for ptc in ptype1v[ibr][iv]:
                 if ptype1v[ibr][iv].count(ptc) != ptype2v[ibr][iv].count(ptc): pmatch = False
             if pmatch: continue
             elif not useDict: return False        
 #If they do not match and useDict=True, generate all possible lists from dictionary entries:
             allptcs = [[ptype1v[ibr][iv]],[ptype2v[ibr][iv]]]
             for allpt in allptcs:
                 ptc0 = copy.deepcopy(allpt[0])
                 for ipt in range(len(ptc0)):
                     if PtcDic.has_key(ptc0[ipt]):
                         for jpt in range(len(allpt)):
                             if allpt[jpt] == []: continue
                             newptc = copy.deepcopy(allpt[jpt])
                             for ptc in PtcDic[ptc0[ipt]]:
                                 newptc[ipt] = ptc
                                 allpt.append(copy.deepcopy(newptc))
                             allpt[jpt] = []
                 while allpt.count([]) > 0: allpt.remove([])
                 
 #Now compare all possibilities:
             match = False
             iA = 0
             while not match and iA < len(allptcs[0]):
                 ptcA = allptcs[0][iA]
                 for ptcB in allptcs[1]:
                     if len(ptcA) != len(ptcB): return False
                     pmatch = True
                     for ptc in ptcA:
                         if ptcA.count(ptc) != ptcB.count(ptc): pmatch = False
                     if pmatch:
                         match = True
                         break
                 iA += 1
             if not match: return False
 
 #if it reached here, entries are similar:
     return True


Rodd={
1000021 : "gluino", 1000022: "N1", 1000023 : "N2", 1000025 : "N3", 1000035 : "N4", 
1000024 : "C1", 1000037 : "C2", 1000039 : "gravitino", 1000001 : "squark", 1000002 : "squark", 
1000003 : "squark", 1000004 : "squark", 2000001 : "squark", 2000002 : "squark", 2000003 : "squark", 
2000004 : "squark", 1000005 : "sbottom", 2000005 : "sbottom", 1000006 : "stop", 2000006 : "stop", 
1000011 : "slepton", 1000013 : "slepton", 1000015 : "stau", 2000011 : "slepton", 2000013 : "slepton", 
2000015 : "stau", 1000012 : "sneutrino", 1000014 : "sneutrino", 1000016 : "sneutrino", 
2000012 : "sneutrino", 2000014 : "sneutrino", 2000016 : "sneutrino", -1000021 : "gluino", 
-1000022: "N1", -1000023 : "N2", -1000025 : "N3", -1000035 : "N4", 
-1000024 : "C1", -1000037 : "C2", -1000039 : "gravitino", 
-1000001 : "squark", -1000002 : "squark", -1000003 : "squark", -1000004 : "squark", 
-2000001 : "squark", -2000002 : "squark", -2000003 : "squark", -2000004 : "squark", 
-1000005 : "sbottom", -2000005 : "sbottom", -1000006 : "stop", -2000006 : "stop", 
-1000011 : "slepton", -1000013 : "slepton", -1000015 : "stau", -2000011 : "slepton", 
-2000013 : "slepton", -2000015 : "stau", -1000012 : "sneutrino", -1000014 : "sneutrino", 
-1000016 : "sneutrino", -2000012 : "sneutrino", -2000014 : "sneutrino", -2000016 : "sneutrino"
}

Reven={
 25 : "higgs", -25: "higgs", 35 : "H0", -35 : "H0", 36 : "A0", -36 : "A0", 37 : "H+", -37 : "H-", 
 23 : "Z", -23 : "Z", 22 : "photon", -22 : "photon", 24 : "W+", -24 : "W-", 16 : "nu", -16 : "nu", 
 15 : "ta-", -15 : "ta+", 14 : "nu", -14 : "nu", 13 : "mu-", -13 : "mu+", 12 : "nu", -12 : "nu", 
 11 : "e-", -11 : "e+", 5 : "b", -5 : "b", 6 : "t+", -6 : "t-", 1 : "jet", 2 : "jet", 3 : "jet", 
 4 : "jet", 21 : "jet", -1 : "jet", -2 : "jet", -3 : "jet", -4 : "jet", -21 : "jet" 
 }

PtcDic={
"e" : ["e+","e-"], "mu" : ["mu+", "mu-"], "ta" : ["ta+","ta-"], "l+" : ["e+","mu+"],"l-" : ["e-","mu-"],
"l" : ["e-","mu-","e+","mu+"], "W" : ["W+","W-"], "t" : ["t+","t-"], "L+" : ["e+","mu+","ta+"], 
"L-" : ["e-","mu-","ta-"], "L" : ["e+","mu+","ta+","e-","mu-","ta-"]
}

