#!/usr/bin/python

# (c) Wolfgang Waltenberger, 2013

""" This module goes thru an LHE particle list and counts 
    the number of susy particles we have in both legs """

def searchParticleByMother ( particles, mother=0 ):
  ret=[]
  for (idx,i) in enumerate(particles):
    if i.mother==mother: ret.append(idx)
  return ret

def reconstructChains ( particles, lsp=1000022 ):
  mothers=searchParticleByMother ( particles, 0 )
  chains=[]
  for mother in mothers:
    currentmother=mother
    onechain=[particles[currentmother].pdgid]
    while abs(particles[currentmother].pdgid)!=lsp:
      daughters=searchParticleByMother ( particles, currentmother )
      for daughter in daughters:
        daughterpdg=particles[daughter].pdgid
        if abs(daughterpdg)>1000000:
          currentmother=daughter
          onechain.append( daughterpdg)
    chains.append(onechain)
  return chains
