#!/usr/bin/env python

import SMSDecomposition

""" stuff to simplify an lhe file """


def newMothers ( after, daughters, mother ):
  for (idx,o) in enumerate(after):
    m=o.mother
    if m in daughters: o.mother=mother
     
def removeDaughters ( after, daughters ):
  """ remove daughters from after """
  ret=[]
  for (idx,o) in enumerate(after):
    if not idx in daughters: ret.append ( o )
  return ret

def charge ( pid ):
  """ get the charge from the particle id """
  apid=abs(pid)
  if apid in [ 1000021, 1000022, 1000023, 1000025, 1000035, 1000039 ]: return 0
  if apid in [ 1000012, 1000014, 1000016, 2000012, 2000014, 2000016 ]: return 0
  if apid in [ 25, 35, 23, 22, 16, 14, 12 ]: return 0
  if apid in [ 1, 3, 5, 1000001, 1000003, 1000005, \
               2000001, 2000003, 2000005 ]: return -1./3.*( pid / apid )
  if apid in [ 2, 4, 6, 1000002, 1000004, 1000006, \
               2000002, 2000004, 2000006 ]: return 2./3.*( pid / apid )
  if apid in [ 11, 13, 15, 1000011, 1000013, 1000015, 2000011, 2000013, 2000015 ]: return -1
  return int ( pid / apid )

def leptonNUmber ( pid ):
  apid=abs(pid)
  if apid in [ 11, 13, 15, 1000011, 1000013, 1000015, 2000011, 2000013, 2000015 ]:
    return int ( pid / apid )

def getLSP ( particles ):
    """ obtain the pdgid of the lsp """
    mylsp=1 # look for the lsp
    m=99999999
    for p in particles:
      if p.pdgid > 1000000 and abs(p.M) < abs(m):
        mylsp=p.pdgid
        m=abs(p.M)
    return mylsp

def createTree ( particles ):
  """ this method creates trees of the particle indices """
  ret={}
  for (index,p) in enumerate(particles):
    mother=p.mother
    if mother<1: continue
    if not ret.has_key ( mother ):
      ret[mother]=[]
    ret[mother].append ( index )
  return ret

def printP ( i, p ):
  print "[lheSimplifier] [%d] %s m=%f mother=%d px=%f" % ( i, SMSDecomposition.particle( p.pdgid ) , p.M, p.mother, p.px )

def simplify ( before, verbose=False ):
  """ this code is intended to perform the simplification """
  if verbose:
    print "[lheSimplifier] [simplifyEvent]"
    print "[lheSimplifier] [before]"
    for (ctr,p) in enumerate(before): printP ( ctr, p )
  tree=createTree ( before )
  after=before
  if verbose:
    print "[lheSimplifier] tree=",tree
  for ( mother, daughters) in tree.items():
    if verbose:
      print 
      print "[lheSimplifier] %d: %s" % ( mother, daughters )
      print "[lheSimplifier] mother mass:",abs(before[mother].M)
    heaviestdaughter=-1
    heaviestdaughtermass=-1.
    charges=[]
    for d in daughters: 
      charges.append ( charge(before[d].pdgid))
      if abs(before[d].M) > heaviestdaughtermass:
        heaviestdaughter=d
        heaviestdaughtermass=abs(before[d].M)
    deltaM=abs(before[mother].M)-heaviestdaughtermass
    if verbose: 
      print "[lheSimplifier] charges of daughters:",charges,"sum:",sum(charges)
      print "[lheSimplifier] heaviest daughter:",heaviestdaughtermass
      print "[lheSimplifier] deltaM=",deltaM
    if deltaM<5. and sum(charges)==0: ## rule number one
      if verbose: print "[lheSimplifier] we simplify!!!!"
      grandmother=after[mother].mother ## we need to copy the granny
      after[mother]=after[heaviestdaughter]
      after[mother].mother=grandmother
      after=removeDaughters ( after, daughters ) ## now remove the daughters from after
      ## if a daughter is mother to some other particle,
      ## the old grannie becomes the new mom
      newMothers ( after, daughters, mother )

  if verbose:
    print
    print "[lheSimplifier] [after]"
    for (ctr,p) in enumerate(after): printP ( ctr, p )

  return after

