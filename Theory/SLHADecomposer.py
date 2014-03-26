#!/usr/bin/env python

"""
.. module:: SLHADecomposer
    :synopsis: SLHA-based SMS decomposition
    
.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
    
"""

import logging
    
def decompose(slhafile,Xsec=None,sigcut=None,DoCompress=False,DoInvisible=False,minmassgap=-1,XsecsInfo=None):
  """Do SLHA-based decomposition.
      FIXME currently using pyslha2 because we need this hack to handle SLHA files with XSECTION blocks.

    :param slhafile: file with mass spectrum and branching ratios and optionally with cross-sections
    :param Xsec: optionally a dictionary with cross-sections for pair production, by default reading the cross sections from the SLHA file.
    :param XsecsInfo: information about the cross-sections (sqrts, order and label). Only relevant for Xsec=None (reading from slha file).\
       If defined as input or in CrossSection.XSectionInfo restricts the cross-sections values in the SLHA file to the ones in XsecsInfo. \
       If not defined, it will be generated from the SLHA file and stored in CrossSection.XSectionInfo.
Only generated if cross-sections are read from SLHA file and not previously created
    :param sigcut: minimum sigma*BR to be generated, by default sigcut = 0.1 fb
    :param DoCompress: turn mass compressed topologies on/off
    :param DoInvisible: turn invisibly compressed topologies on/off
    :param minmassgap: maximum value for considering two R-odd particles degenerate (only revelant for DoCompress=True)
    :returns: a TopologyList. 
  """
  import sys, os, copy

  workdir = os.getcwd() 
  pyslhadir = workdir + "/pyslha-1.4.3"
  sys.path.append(pyslhadir)
  import TopologyBuilder, SMSDataObjects, ParticleNames
  from Tools.PhysicsUnits import addunit, rmvunit
  import pyslha2 as pyslha
  import CrossSection
  
  log = logging.getLogger(__name__)

  if DoCompress and rmvunit(minmassgap,'GeV') == -1: 
    print "[SLHAdecomposer] Please, set minmassgap"
    return False

#sigcut by default 0.1 fb
  if not rmvunit(sigcut, 'fb'):
    sigcut = addunit(0.1, 'fb')

 #creates Xsec dictionary if Xsec=None and store cross-section information if not previously defined
  if not Xsec:
    XsecsInfoFile = CrossSection.XSecInfoList('')  #To store information about all cross-sections in the SLHA file
    Xsec = {}
    slha = open(slhafile, 'r')
    lines = slha.readlines()
    xsecblock = False
    for l in lines:
      if l.startswith("#") or len(l.replace('\n',''))<2: continue
      if 'XSECTION' in l:
        xsecblock = True
        sqrtS =  eval(l.split()[1])/1000.    #Values in the SLHA file are in GeV
        pids = (eval(l.split()[5]),eval(l.split()[6]))
        continue
      if not xsecblock: continue  #ignores other entries
      cs_order = eval(l.split()[1])
      cs = addunit(eval(l.split()[6]),'fb')
      wlabel = str(int(sqrtS))+' TeV'
      if cs_order == 0: wlabel += ' (LO)'
      elif cs_order == 1: wlabel += ' (NLO)'
      elif cs_order == 2: wlabel += ' (NLL)'
      else:
        print '[SLHADecomposer] unknown QCD order in XSECTION line', l
        return False
      xsInfo = CrossSection.SingleXSecInfo()
      xsInfo.sqrts = addunit(sqrtS,'TeV')
      xsInfo.order = cs_order
      xsInfo.label = wlabel
      if not wlabel in Xsec.keys():
        Xsec[wlabel] = {}
        XsecsInfoFile.xsecs.append(xsInfo)
      Xsec[wlabel][pids] = cs
   

    if XsecsInfo is None:
      try:
        XsecsInfo = CrossSection.XSectionInfo  #Check if cross-section information has been defined
      except:
        pass
    if not XsecsInfo:
      XsecsInfo = XsecsInfoFile              #Store information from file
      CrossSection.XSectionInfo = XsecsInfo
      log.warning ( "Cross-section information not found. Using values from SLHA file" )
    else:
      for xsec in XsecsInfoFile.xsecs:
        if not xsec in XsecsInfo.xsecs: Xsec.pop(xsec.label)    #Remove entries which do not match the previously defined cross-sections


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
      if not BaseEl.momID: continue
      ptc = BaseEl.momID.pop()
      weight = WList[iel]
      
      if len(BRdic[ptc]) == 0:     # Stable final state (LSP)
        BaseEl.momID = BaseEl.momID[0]
        FinalList.append(BaseEl)
        WFinal.append(weight)
        continue
      
      for BR in BRdic[ptc]:
        if BR.br < 0.:
          log.warning ("Ignoring negative BRs" )
          continue
        NewEl = copy.deepcopy(BaseEl)
        vertparts = []
        mass = []
        for x in BR.ids:
          if x in ParticleNames.Reven:
            vertparts.append(ParticleNames.Reven[x])
          else:
            mass.append(Massdic[x])
            NewEl.momID.append(x)
          
        NewEl.particles.append(vertparts)
        if len(mass) == 1: NewEl.masses.append(mass[0])
        
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
    for iel,El1 in enumerate(FinalList):
      for jel,El2 in enumerate(FinalList):

        if El1.momID == ptcs[0] and El2.momID == ptcs[1]:
          Els = SMSDataObjects.EElement()    
          Els.B = [copy.deepcopy(El1),copy.deepcopy(El2)]
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

