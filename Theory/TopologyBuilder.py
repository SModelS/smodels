#!/usr/bin/env python

"""
.. module:: TopologyBuilder
    :synopsis: A builder to create GTop objects from e.g. SMSEvents
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
    
"""

def compressTopology(ETopList,DoCompress,DoInvisible,minmassgap):
  """ Given a topology list, keep compressing the element till it 
      can be compressed no more.
      Returns a list with the old topologies and the compressed ones.
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


