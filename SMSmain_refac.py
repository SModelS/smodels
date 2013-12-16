#!/usr/bin/python

import sys
from prettytable import PrettyTable
from Theory import SLHADecomposer, SLHATools, CrossSection
from Tools.PhysicsUnits import addunit, rmvunit
from Tools import SMSPrettyPrinter, VariousHelpers
from Tools.SMSPrettyPrinter import wrap
from Tools.VariousHelpers import logging
from Experiment import SMSAnalysisFactory
from Theory.theoryPrediction import PredictionForAnalysis

# useXsec = CrossSection.XSectionInfo()
# useXsec.sqrts = addunit(8,'TeV')
# useXsec.order = 2
# useXsec.label = 'tev8'
# CrossSection.UseXSecs = [useXsec]

listOfAnalyses = SMSAnalysisFactory.load()

printer=SMSPrettyPrinter.SMSPrettyPrinter()

slhafile = "slha/andrePT4.slha"
nevts = 10000

Xsec=None
SLHATools.writeXSecToSLHAFile(slhafile,nevts,printLHE=False)
DoCompress = True
DoInvisible = True
minmassgap = addunit(5.,'GeV')
sigmacut = addunit(0.1,'fb')
SMSTopList = SLHADecomposer.decompose(slhafile,sigmacut,DoCompress,DoInvisible,minmassgap)

EvTop_table = PrettyTable(["Topology","#Vertices", "#Insertions", "#Elements", "Sum of weights"])
EvElement_table = PrettyTable(["Topology","Element","Particles B[0]","Particles B[1]", "Masses B[0]","Masses B[1]","Element Weight"])

eltot = 0
totweight = []
#Print Results:
# for i in range(len(SMSTopList)):
for (i,topo) in enumerate(SMSTopList):
  sumw = topo.getTotalWeight().getDictionary()
  EvTop_table.add_row([i,topo.vertnumb,topo.vertparts,len(topo.ElList),wrap(printer.pformat(sumw),width=30)])
  eltot += len(topo.ElList)

 
      
#Print element list for Topology[i]:  
  if i == 0:       
    for j,el in enumerate(topo.ElList):
      EvElement_table.add_row([i,j,el.getParticles()[0],el.getParticles()[1],wrap(printer.pformat(el.getMasses()[0]),width=25),wrap(printer.pformat(el.getMasses()[1]),width=25),wrap(printer.pformat(el.weight.getDictionary()),width=30)])
    EvElement_table.add_row(["---","---","---","---","---","---","---"])  

     
print "Number of Global topologies = ",len(SMSTopList)      
print(EvTop_table)
print "Total Number of Elements = ",eltot
print "Total weight = ",SMSTopList.getTotalWeight()
# print(EvElement_table)

print '\n \n \n'

 
for iana,ana in enumerate(listOfAnalyses):    
    print ana.label,ana.listOfPlots[0].label,ana.sqrts    
    pred = PredictionForAnalysis(ana)
    pred.getTheoryXSecs(SMSTopList)
    pred.clusterElements()
     
    for el in pred.analysisTopo.ElList:
        print el,el.getMasses()
        print el.weight


sys.exit()

