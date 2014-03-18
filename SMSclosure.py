#!/usr/bin/python


import sys, glob, os
smsdir = os.getcwd() + '/SMSfrozen/'
sys.path.append(smsdir)

from theory import SLHADecomposer, CrossSection, ClusterTools
from tools.physicsUnits import addunit, rmvunit
from experiment import SMSAnalysisFactory

XSectionInfo = CrossSection.XSecInfoList('8 TeV (NLL)')
CrossSection.XSectionInfo = XSectionInfo
#Basic input for decomposition:
if os.path.isfile(sys.argv[1]):
  slhafileList = [sys.argv[1]]  
else:
  slhafileList = glob.glob(sys.argv[1]+'*.slha')
  outfile = open(sys.argv[4],'w')


for slhafile in slhafileList:
  print slhafile
  res = 0.
  DoCompress = False
  DoInvisible = False
  minmassgap = addunit(10.,'GeV')
  sigmacut = addunit(0.000001,'fb')
  SMSTopList = SLHADecomposer.decompose(slhafile,None,sigmacut,DoCompress,DoInvisible,minmassgap)
#Get total cross-section
  totweight = []
  for topo in SMSTopList:
    weightlist = [el.weight for el in topo.elList]
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
  Mmass = Analysis.Top.elList[0].MassWeightList[0].mass[0][0].asNumber()
  LSPmass = Analysis.Top.elList[0].MassWeightList[0].mass[0][-1].asNumber()
  if len(Analysis.ResultList) > 1:
    print 'more than one cluster in',slhafile
    sys.exit()
  elif len(Analysis.ResultList) == 0:
    res = -1.
    limit = 1.  #dummy
    cond = -1. #dummy
  else:      
    cluster = Analysis.ResultList[0]
    res = rmvunit(cluster.result_dic.values()[0].values()[0],'fb')
    limit = rmvunit(cluster.explimit,'fb')
    cond = cluster.getMaxCondition()
#Print result:
  outfile.write(str(Mmass)+" "+str(LSPmass)+" "+str(res)+" "+str(limit)+" "+str(cond)+" "+str(rmvunit(totxsec.values()[0],'fb'))+"\n")

outfile.close()
sys.exit()
