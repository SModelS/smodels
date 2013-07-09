#!/usr/bin/env python

"""
.. module:: SLHADecomposer
    :synopsis: SLHA-based SMS decomposition
    
.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
    
"""
    
def decompose(slhafile,Xsec=None,sigcut=None,DoCompress=False,DoInvisible=False,minmassgap=-1):
  """ Do SLHA-based decomposition.
      FIXME currently using pyslha2 because we need this hack to handle SLHA files with XSECTION blocks.

    :param slhafile: file with mass spectrum and branching ratios and optionally with cross-sections
    :param Xsec: optionally a dictionary with cross-sections for pair production, by default reading the cross sections from the SLHA file.
    :param sigcut: minimum sigma*BR to be generated, by default sigcut = 0.1 fb
    :param DoCompress: FIXME
    :param DoInvisible: FIXME
    :param minmassgap: FIXME
    :returns: a TopologyList. 
  """
  import sys, os, copy

  workdir = os.getcwd() 
  pyslhadir = workdir + "/pyslha-1.4.3"
  sys.path.append(pyslhadir)
  import TopologyBuilder, SMSDataObjects, ParticleNames
  from Tools.PhysicsUnits import addunit, rmvunit
  import pyslha2 as pyslha

  if DoCompress and rmvunit(minmassgap,'GeV') == -1: 
    print "[SLHAdecomposer] Please, set minmassgap"
    return False

#sigcut by default 0.1 fb
  if not rmvunit(sigcut, 'fb'):
    sigcut = addunit(0.1, 'fb')

#creates Xsec dictionary
  if not Xsec:
    Xsec = {'7 TeV (LO)':{}, '8 TeV (LO)':{}, '7 TeV (NLL)':{}, '8 TeV (NLL)':{}}
    slha = open(slhafile, 'r')
    lines = slha.readlines()
    currentblock = None
    for l in lines:
      if l.startswith("#"):
         continue
      if 'XSECTION' in l:
         currentblock = 'XSECTION'
         sqrt = l.split()[1]
         pids = [int(l.split()[2]), int(l.split()[3])]
         pids.sort()
      elif ('BLOCK' in l) or ('DECAY' in l):
         currentblock = 'BLOCK'
      else:
         if currentblock == 'BLOCK':          #ignores BLOCK and DECAY entries
           continue
         elif currentblock == 'XSECTION':
           if l.split()[1] == '0':
             key = (pids[0],pids[1])
             Xsec['%s TeV (LO)'%sqrt][key] = addunit(float(l.split()[6]),'fb')
           elif l.split()[1] =='2':
             key = (pids[0],pids[1])
             Xsec['%s TeV (NLL)'%sqrt][key] = addunit(float(l.split()[6]),'fb')
           else:
             print '[SLHADecomposer] unknown sqrt in XSECTION line', l
         else:
           print '[SLHADecomposer] unknown Block in XSECTION line', l

  print Xsec
#Read SLHA file
  res = pyslha.readSLHAFile(slhafile)
  
#Get list of particles with maximum production cross-section above sigcut and maximum cross-sections
  Pdic = {}
  for k in Xsec.keys():
    for k2 in Xsec[k]:
      if Xsec[k][k2] > sigcut:
        for ip in range(2):
          if not Pdic.has_key(k2[ip]) or Pdic[k2[ip]] < Xsec[k][k2]: 
            Pdic.update({k2[ip] : Xsec[k][k2]})

#Get list of branching ratios for all particles:
  BRdic = {}
  for k in res[1].keys():
    brs = copy.deepcopy(res[1][abs(k)].decays)
    for br in brs:
      br.ids = [-x for x in br.ids]
    BRdic.update({k : res[1][abs(k)].decays, -k : brs})

#Get mass list for all particles
  Massdic = {}
  for k in res[1].keys():
    if k and res[1][k].mass != None:
      Massdic.update({k : addunit(abs(res[1][k].mass),'GeV'), -k : addunit(abs(res[1][k].mass),'GeV')})
      
#Loop over all particles and generate all possible 1branch-elements with sigmamax*BR > sigcut
  ElList = []
  WList = []
  for ptc in Pdic.keys():
    NewEl = SMSDataObjects.BElement()
    NewEl.momID = [ptc,ptc]
    NewEl.masses.append(Massdic[ptc])
    weight = Pdic[ptc]
    ElList.append(NewEl)
    WList.append(weight)
    
  FinalList = []
  WFinal = []
  newel = True  
  while newel:
    newel = False
    NewList = []
    NewWeight = []
    for iel in range(len(ElList)):
      BaseEl = ElList[iel]
      ptc = BaseEl.momID.pop()
      weight = WList[iel]
      
      if len(BRdic[ptc]) == 0:     # Stable final state (LSP)
        BaseEl.momID = BaseEl.momID[0]
        FinalList.append(BaseEl)
        WFinal.append(weight)
        continue
      
      for BR in BRdic[ptc]:
        NewEl = copy.deepcopy(BaseEl)
        vertparts = []
        mass = []
        for x in BR.ids:
          if x in ParticleNames.Reven:
            vertparts.append(ParticleNames.Reven[x])
          elif x in ParticleNames.Rodd:
            mass.append(Massdic[x])
            NewEl.momID.append(x)
          else:
            print '[SLHAdecomposer] unknown particle:',x
            return False
          
        NewEl.particles.append(vertparts)
        if len(mass) == 1:
          NewEl.masses.append(mass[0])
        else:
          print '[SLHAdecomp] unknown decay (R-parity violation?)'
          return False
        
        if weight*BR.br > sigcut:
          NewList.append(NewEl)
          NewWeight.append(weight*BR.br)
          newel = True
                
    if newel:
      ElList = copy.deepcopy(NewList)
      WList = copy.deepcopy(NewWeight)
            

  SMSTopList = SMSDataObjects.TopologyList()
#Combine 1branch elements according to production cross-section:
  for ptcs in Xsec[Xsec.keys()[0]].keys():
    for iel in range(len(FinalList)):
      for jel in range(len(FinalList)):
        if ptcs[0] == ptcs[1] and jel < iel: continue     #Avoid double counting

        if FinalList[iel].momID == ptcs[0] and FinalList[jel].momID == ptcs[1]:
          Els = SMSDataObjects.EElement()    
          Els.B = [copy.deepcopy(FinalList[iel]),copy.deepcopy(FinalList[jel])]
          weight = {}
          for w in Xsec.keys():
            weight.update({w : Xsec[w][ptcs]*WFinal[iel]*WFinal[jel]/(Pdic[ptcs[0]]*Pdic[ptcs[1]])})
          Els.weight = weight          
        
          if max(weight.values()) < sigcut: continue
        
          Einfo = Els.getEinfo()
          Top = SMSDataObjects.GTop()
          Top.vertnumb = Einfo["vertnumb"]
          Top.vertparts = Einfo["vertparts"]
          Top.ElList.append(Els)
          for ib in range(len(Top.vertnumb)):
            if len(Top.vertparts[ib]) != Top.vertnumb[ib]:
              print '[SLHAdecomposer] error creating topology'
              return False

          TopList = []
          TopList.append(Top) 
          # Do compression:          
          if DoCompress or DoInvisible:
            FinalTopList = TopologyBuilder.compressTopology(TopList,DoCompress,DoInvisible,minmassgap)
          else:
            FinalTopList = TopList
          # Add topologies to topology list:
          SMSTopList.addList ( FinalTopList )

      
  return SMSTopList    

