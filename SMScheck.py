#!/usr/bin/python

import sys, glob, os
from Theory import SLHADecomposer, CrossSection, SLHATools
from Tools.PhysicsUnits import addunit
from Experiment import SMSAnalysisFactory, LimitGetter

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
  
#  infile = open(slhafile,'r')
#  lines = infile.readlines()
#  for line in lines:
#    if '1000021' in line and len(line.split()) == 2: mM = eval(line.split()[1])  #Gluino
#    if '1000022' in line and len(line.split()) == 2: m = eval(line.split()[1])   #N1
#    if '1000024' in line and len(line.split()) == 2: mI = eval(line.split()[1])  #C1
#  
#  x = (mI-m)/(mM-m)
#  if mI-m < 81.:
#    print slhafile
#    os.remove(slhafile)
#  elif abs(x-0.7) < 0.01 and not '070' in slhafile:
#    print slhafile,x
#  elif abs(x-0.5) < 0.01 and '070' in slhafile:
#    print slhafile,x
#  elif abs(x-0.7) > 0.01 and abs(x-0.5) > 0.01:
#    print slhafile,x    
#  continue
    
  
#  fname = slhafile[slhafile.rfind('/'):]
#  mM = eval(fname[fname.find('_')+1:fname.rfind('_')])
#  m = eval(fname[fname.rfind('_')+1:fname.rfind('.')])
#  if mM-m < 116:
#    print slhafile
#    os.remove(slhafile)
#  continue

  

#Generate events and compute cross-sections:
  nevts = 10000
  SLHATools.writeXSecToSLHAFile(slhafile,nevts,printLHE=False)
  sys.exit()

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
