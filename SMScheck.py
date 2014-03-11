#!/usr/bin/python

import sys, glob, os
from theory import SLHADecomposer, CrossSection, SLHATools
from tools.PhysicsUnits import addunit
from experiment import SMSAnalysisFactory, LimitGetter

XSectionInfo = CrossSection.XSecInfoList('8 TeV (NLL)')
CrossSection.XSectionInfo = XSectionInfo
#Basic input for decomposition:
if os.path.isfile(sys.argv[1]):
  slhafileList = [sys.argv[1]]
else:
  slhafileList = glob.glob(sys.argv[1]+'*')
  excl = open("excluded.dat","w")
  allow = open("allowed.dat","w")



for slhafile in slhafileList:
  print slhafile
 

#Generate events and compute cross-sections:
  nevts = 10000
  SLHATools.writeXSecToSLHAFile(slhafile,nevts,printLHE=False)
  continue

  DoCompress = False
  DoInvisible = False
  minmassgap = addunit(10.,'GeV')
  sigmacut = addunit(0.01,'fb')
  SMSTopList = SLHADecomposer.decompose(slhafile,None,sigmacut,DoCompress,DoInvisible,minmassgap)
#Load analysis
  ListOfAnalyses = SMSAnalysisFactory.load(sys.argv[2],sys.argv[3])
  for Analysis in ListOfAnalyses: Analysis.add(SMSTopList)
#Compute theoretical predictions to analyses results:
  for Analysis in ListOfAnalyses: Analysis.computeTheoryPredictions()
  

#Simple output example:  
  for Analysis in ListOfAnalyses:
    if len(slhafileList) == 1:
      print Analysis.label
#      for Top in SMSTopList:
#        for el in Top.ElList:
#          print el.B[0].masses,el,el.weight
      for cluster in Analysis.ResultList:
        print cluster.mass,cluster.result_dic.values(),cluster.explimit
      continue
    
    if len(Analysis.ResultList) > 1:
      print 'more than one cluster in',slhafile
      sys.exit()
    elif len(Analysis.ResultList) == 0:
      print 'skip',slhafile
      continue  
    cluster = Analysis.ResultList[0]
    if cluster.result_dic.values()[0].values()[0]/cluster.explimit > 1.:
      excl.write(str(cluster.mass[0][0].asNumber())+" "+str(cluster.mass[0][-1].asNumber())+"\n")
    else:
      allow.write(str(cluster.mass[0][0].asNumber())+" "+str(cluster.mass[0][-1].asNumber())+"\n")


sys.exit()
