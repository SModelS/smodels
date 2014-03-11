#!/usr/bin/python

import sys
from theory import SLHADecomposer
from tools.PhysicsUnits import addunit
from experiment import SMSAnalysisFactory

#Basic input for decomposition:
slhafile = "slha/andrePT4.slha"
DoCompress = True
DoInvisible = True
minmassgap = addunit(10.,'GeV')
sigmacut = addunit(0.01,'fb')
SMSTopList = SLHADecomposer.decompose(slhafile,None,sigmacut,DoCompress,DoInvisible,minmassgap)

#Load analysis
ListOfAnalyses = SMSAnalysisFactory.load()
for Analysis in ListOfAnalyses: Analysis.add(SMSTopList)
#Compute theoretical predictions to analyses results:
for Analysis in ListOfAnalyses: Analysis.computeTheoryPredictions()
  

#Simple output example:  
for Analysis in ListOfAnalyses:
  for cluster in Analysis.ResultList:
    print cluster.mass,Analysis.label,cluster.explimit

sys.exit()
