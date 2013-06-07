#!/usr/bin/env python

"""
.. module:: TopologyBuilder
    :synopsis: ...
    
.. moduleauthor:: someone <email@example.com>
    
"""
    
""" A builder to create GTop objects from e.g. SMSEvents """

def fromEvent( Event, weight = {}, DoCompress=False, DoInvisible=False, \
                       minmassgap=None ):
  """ Create a topology (GTop object) from Event.
      If the DoCompress and/or DoInvisible flags are on, also generate
  compressed topologies with small mass gaps and/or neutrinos emitted
  in the last step of the cascade ("effective LSP"). """
  from SMSDataObjects import GTop, EElement, BElement
  from ParticleNames import Rodd, Reven, ptype
  from Tools.PhysicsUnits import addunit, rmvunit
  PList=Event.particles
  if len(PList)==0: return None
    
  ETopList = []
  ETop = GTop()
  
  momspdg = [0,0]
  mompos = [0,0]
  imom = 0
  nptcs = 0    #Particle counter just for sanity checks
  
  minmassgap=rmvunit(minmassgap,"GeV")
  if DoCompress and minmassgap==None: 
    print "[TopologyBuilder.py] Please, set minmassgap"
    return None
  if DoCompress and minmassgap==-1: 
    print "[TopologyBuilder.py] Please, set minmassgap"
    return None
  
#First get Mothers:  
  for i in range(len(PList)):
    if PList[i].moms[0] == 1 or PList[i].moms[1] == 1:
      momspdg[imom] = PList[i].pdg
      mompos[imom] = i
      imom +=1

  NewEl = EElement()
#Each mother = different branch. Loop over branchs: 
  for ib in range(2):
    mother = mompos[ib]
    newmom = mompos[ib]
    El = BElement()
        
    nptcs += 1
    
    ndaugs = 2
    while ndaugs > 0:
      ndaugs = 0
      El.particles.append([])
      for i in range(len(PList)):
        if PList[i].moms[0] != mother+1 and PList[i].moms[1] != mother+1: continue
        if abs(PList[i].pdg) in Rodd:
          newmom = i
          ndaugs += 1
        elif abs(PList[i].pdg) in Reven:
          pname = ptype(PList[i].pdg)
          El.particles[len(El.particles)-1].append(pname)
          ndaugs += 1
        else:  # dunno the particle so it must be a BSM thingie
          newmom = i
          ndaugs += 1

      El.masses.append(addunit(PList[mother].mass,'GeV'))
      mother = newmom

    if El.particles[len(El.particles)-1] == []:
      El.particles.pop()   #Remove last empty insertion if LSP is stable
    El.momID = momspdg[ib]
    NewEl.B.append(El)

  import copy
  NewEl.weight = copy.deepcopy(weight)
  Einfo = NewEl.getEinfo()
  ETop.ElList.append(NewEl)
  ETop.vertnumb = Einfo["vertnumb"]
  ETop.vertparts = Einfo["vertparts"]


  ETopList.append(ETop)    
  
  #Do compression:
  if DoCompress or DoInvisible:
    FinalTopList = compressTopology(ETopList,DoCompress,DoInvisible,minmassgap)
  else:
    FinalTopList = ETopList

  return FinalTopList
    

def compressTopology(ETopList,DoCompress,DoInvisible,minmassgap):
  """ Given a topology list, keep compressing the element it 
      can be compressed no more.
      Returns a list with the old toplogies and the compressed ones.
      To avoid double counting the input list should have a single element. """
  from Tools.PhysicsUnits import rmvunit
  minmassgap=rmvunit ( minmassgap, "GeV" )

#Keep compressing the topologies generated so far until no new compressions can happen:
  if len(ETopList) > 0:
    added = True
  else:
    added = False
  while added:    
    added = False
#Check for mass compressed topologies    
    if DoCompress:
      for Top in ETopList:
        ETopComp = Top.massCompressedTopology(minmassgap)
        if ETopComp:
          exists = False
          for Topp in ETopList:
            if Topp.isEqual(ETopComp): exists = True 
          if not exists:   #Avoid double counting (conservative)
            ETopList.append(ETopComp)
            added = True
      
#Check for invisible compressed topologies 
#(look for effective LSP, such as LSP + neutrino = LSP')      
    if DoInvisible:
      for Top in ETopList:
        ETopInComp = Top.invisibleCompressedTopology()
        if ETopInComp:
          exists = False
          for Topp in ETopList:
            if Topp.isEqual(ETopInComp): exists = True
          if not exists:   #Avoid double counting (conservative)
            ETopList.append(ETopInComp)
            added = True
  return ETopList


