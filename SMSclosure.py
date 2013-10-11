#!/usr/bin/python

import sys, glob, os
from Theory import SLHADecomposer, CrossSection, ClusterTools
from Tools.PhysicsUnits import addunit, rmvunit
from Experiment import SMSAnalysisFactory

XSectionInfo = CrossSection.XSecInfoList('8 TeV (NLL)')
CrossSection.XSectionInfo = XSectionInfo
#Basic input for decomposition:
if os.path.isfile(sys.argv[1]):
  slhafileList = [sys.argv[1]]  
else:
  slhafileList = glob.glob(sys.argv[1]+'*.slha')
  outfile = open(sys.argv[1][sys.argv[1].rfind('/')+1:]+'.dat','w')


for slhafile in slhafileList:
  res = 0.
  DoCompress = False
  DoInvisible = False
  minmassgap = addunit(10.,'GeV')
  sigmacut = addunit(0.000001,'fb')
  SMSTopList = SLHADecomposer.decompose(slhafile,None,sigmacut,DoCompress,DoInvisible,minmassgap)
#Get total cross-section
  totweight = []
  for topo in SMSTopList:
    weightlist = [el.weight for el in topo.ElList]
    sumw = ClusterTools.sumweights(weightlist)
    totweight.append(sumw)
  totxsec = ClusterTools.sumweights(totweight)
  
#Load analysis
  ListOfAnalyses = SMSAnalysisFactory.load(sys.argv[2],sys.argv[3])
  for Analysis in ListOfAnalyses: Analysis.add(SMSTopList)
#Compute theoretical predictions to analyses results:
  for Analysis in ListOfAnalyses: Analysis.computeTheoryPredictions()
  

#Get the result:
  Analysis = ListOfAnalyses[0]
  Mmass = Analysis.Top.ElList[0].MassWeightList[0].mass[0][0].asNumber()
  LSPmass = Analysis.Top.ElList[0].MassWeightList[0].mass[0][-1].asNumber()
  if len(Analysis.ResultList) > 1:
    print 'more than one cluster in',slhafile
    sys.exit()
  elif len(Analysis.ResultList) == 0:
    res = -1.
    limit = 1.  #dummy
  else:      
    cluster = Analysis.ResultList[0]
    res = rmvunit(cluster.result_dic.values()[0].values()[0],'fb')
    limit = rmvunit(cluster.explimit,'fb')
    if res/rmvunit(totxsec.values()[0],'fb') < 0.03: print slhafile
#Print result:
  try:
    outfile.write(str(Mmass)+" "+str(LSPmass)+" "+str(res)+" "+str(limit)+" "+str(rmvunit(totxsec.values()[0],'fb'))+"\n")
  except:    
    print str(Mmass)+" "+str(LSPmass)+" "+str(res)+" "+str(limit)+" "+str(rmvunit(totxsec.values()[0],'fb'))+"\n"

try: outfile.close()
except: pass  
sys.exit()
