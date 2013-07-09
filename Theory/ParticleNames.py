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
  """
  p=int(pdg)
  if p in Rodd: return Rodd[p]
  if p in Reven:
    return Reven[p]
  else:
    return False

def simParticles(ptype1,ptype2,useDict=True):
  """ Compares 2 particle names or 2 nested name arrays. Allows for dictionary labels
  (Ex: L = l, l+ = l, l = l-,...) 
  For the last nested level ignore particle ordering """
  if len(ptype1) != len(ptype2): return False
  import copy

  ptype1v = [[ptype1]]
  ptype2v = [[ptype2]]


#First flatten nested arrays until next-to-last level:
  isNested = True
  while isNested:
    newptype1v = []
    newptype2v = []
    if len(ptype1v) != len(ptype2v): return False
    for i in range(len(ptype1v)):
      if type(ptype1v[i]) == type(list()):
        if len(ptype1v[i]) != len(ptype2v[i]): return False
        for j in range(len(ptype1v[i])):
          newptype1v.append(ptype1v[i][j])
          newptype2v.append(ptype2v[i][j])
      else:
        newptype1v.append(ptype1v[i])
        newptype2v.append(ptype2v[i])

    ptype1v = newptype1v
    ptype2v = newptype2v
    isNested = False
    for i in range(len(ptype1v)):
      if len(ptype1v[i]) != len(ptype2v[i]): return False
      if len(ptype1v[i]) == 0: continue   #Empty list
      if type(ptype1v[i]) == type(list()) and type(ptype1v[i][0]) == type(list()): isNested = True
      if type(ptype2v[i]) == type(list()) and type(ptype2v[i][0]) == type(list()): isNested = True

  if len(ptype1v) != len(ptype2v): return False

#Compare last level lists one by one, ignoring the order:
  for i in range(len(ptype1v)):
    if len(ptype1v[i]) != len(ptype2v[i]): return False


#Check  if lists match, ignoring possible dictionary entries
    pmatch = True
    for ptc in ptype1v[i]:
      if ptype1v[i].count(ptc) != ptype2v[i].count(ptc): pmatch = False
    if pmatch: continue
    elif not useDict: return False

#If they do not match and useDict=True, generate all possible lists from dictionary entries:
    allptcs = [[ptype1v[i]],[ptype2v[i]]]
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
 25 : "higgs", -25: "higgs", 35 : "H0", -35 : "H0", 36 : "A0", -36 : "A0", 37 : "H+", -37 : "H-", 23 : "Z", -23 : "Z", 22 : "photon", -22 : "photon", 24 : "W+", -24 : "W-", 16 : "nu", -16 : "nu", 15 : "ta-", -15 : "ta+", 14 : "nu", -14 : "nu", 13 : "mu-", -13 : "mu+", 12 : "nu", -12 : "nu", 11 : "e-", -11 : "e+", 5 : "b", -5 : "b", 6 : "t+", -6 : "t-", 1 : "jet", 2 : "jet", 3 : "jet", 4 : "jet", 21 : "jet", -1 : "jet", -2 : "jet", -3 : "jet", -4 : "jet", -21 : "jet"
 }

PtcDic={
"e" : ["e+","e-"], "mu" : ["mu+", "mu-"], "ta" : ["ta+","ta-"], "l+" : ["e+","mu+"],"l-" : ["e-","mu-"],"l" : ["e-","mu-","e+","mu+"], "W" : ["W+","W-"], "t" : ["t+","t-"], "L+" : ["e+","mu+","ta+"], "L-" : ["e-","mu-","ta-"], "L" : ["e+","mu+","ta+","e-","mu-","ta-"]
}


