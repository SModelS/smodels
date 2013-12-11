#!/usr/bin/env python

"""
.. module:: ParticleNames
    :synopsis: methods for getting particle names from pdg ids, and other helpers

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

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
  import copy,sys
  
  wrongFormat = False
  if type(ptype1) != type(ptype2): return False
#Check for nested arrays (should be in the standard notation [[[]],[[]]]):  
  if type(ptype1) == type([]):
    if len(ptype1) != len(ptype2): return False
    for ib,br in enumerate(ptype1):
      if type(br) != type(ptype2[ib]) or type(br) != type([]): wrongFormat = True
      if len(ptype1[ib]) != len(ptype2[ib]): return False  #Check number of vertices in branch
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
    ptype1v = copy.deepcopy([[[ptype1]],[[ptype1]]])
    ptype2v = copy.deepcopy([[[ptype2]],[[ptype2]]])
  else:
    ptype1v = copy.deepcopy(ptype1)
    ptype2v = copy.deepcopy(ptype2)

  for ibr,br in enumerate(ptype1v): #Loop through branches
    for iv,vt in enumerate(br): #Loop over vertices  
#Check  if lists match, ignoring possible dictionary entries
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
1000021 : "gluino", 1000022: "N1", 1000023 : "N2", 1000025 : "N3", 1000035 : "N4", 1000024 : "C1", 1000037 : "C2", 1000039 : "gravitino",  1000001 : "squark", 1000002 : "squark", 1000003 : "squark", 1000004 : "squark", 2000001 : "squark", 2000002 : "squark", 2000003 : "squark", 2000004 : "squark", 1000005 : "sbottom", 2000005 : "sbottom", 1000006 : "stop", 2000006 : "stop", 1000011 : "slepton", 1000013 : "slepton", 1000015 : "stau", 2000011 : "slepton", 2000013 : "slepton", 2000015 : "stau", 1000012 : "sneutrino", 1000014 : "sneutrino", 1000016 : "sneutrino", 2000012 : "sneutrino", 2000014 : "sneutrino", 2000016 : "sneutrino", -1000021 : "gluino", -1000022: "N1", -1000023 : "N2", -1000025 : "N3", -1000035 : "N4", -1000024 : "C1", -1000037 : "C2", -1000039 : "gravitino",  -1000001 : "squark", -1000002 : "squark", -1000003 : "squark", -1000004 : "squark", -2000001 : "squark", -2000002 : "squark", -2000003 : "squark", -2000004 : "squark", -1000005 : "sbottom", -2000005 : "sbottom", -1000006 : "stop", -2000006 : "stop", -1000011 : "slepton", -1000013 : "slepton", -1000015 : "stau", -2000011 : "slepton", -2000013 : "slepton", -2000015 : "stau", -1000012 : "sneutrino", -1000014 : "sneutrino", -1000016 : "sneutrino", -2000012 : "sneutrino", -2000014 : "sneutrino", -2000016 : "sneutrino"

}

Reven={
 25 : "higgs", -25: "higgs", 35 : "H0", -35 : "H0", 36 : "A0", -36 : "A0", 37 : "H+", -37 : "H-", 23 : "Z", -23 : "Z", 22 : "photon", -22 : "photon", 24 : "W+", -24 : "W-", 16 : "nu", -16 : "nu", 15 : "ta-", -15 : "ta+", 14 : "nu", -14 : "nu", 13 : "mu-", -13 : "mu+", 12 : "nu", -12 : "nu", 11 : "e-", -11 : "e+", 5 : "b", -5 : "b", 6 : "t+", -6 : "t-", 1 : "jet", 2 : "jet", 3 : "jet", 4 : "jet", 21 : "jet", -1 : "jet", -2 : "jet", -3 : "jet", -4 : "jet", -21 : "jet", 5000022 : "X", -5000022 : "X"
 }

PtcDic={
"e" : ["e+","e-"], "mu" : ["mu+", "mu-"], "ta" : ["ta+","ta-"], "l+" : ["e+","mu+"],"l-" : ["e-","mu-"],"l" : ["e-","mu-","e+","mu+"], "W" : ["W+","W-"], "t" : ["t+","t-"], "L+" : ["e+","mu+","ta+"], "L-" : ["e-","mu-","ta-"], "L" : ["e+","mu+","ta+","e-","mu-","ta-"]
}


