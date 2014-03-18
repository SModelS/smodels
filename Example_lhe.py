#!/usr/bin/python

import sys
from prettytable import PrettyTable
from theory import lheDecomposer, SLHADecomposer, XSecComputer, ClusterTools, CrossSection, SLHATools
from tools.physicsUnits import addunit, rmvunit
from tools import SMSPrettyPrinter, VariousHelpers
from tools.SMSPrettyPrinter import wrap
from Experiment import SMSAnalysisList, SMSAnalysisFactory, LimitGetter, SMSHelpers


printer=SMSPrettyPrinter.SMSPrettyPrinter()

#Example of how to define cross-sections. If not defined, default values
#will be generated and stored in CrossSection.XSectionInfo the first
#time this information is needed 
XSectionInfo = CrossSection.XSecInfoList('8 TeV (NLL)')
CrossSection.XSectionInfo = XSectionInfo

#Input file and number of events to be used (if =None, use all)
lhefile = "lhe/ued_1.lhe"
nevts=None
#Creat analyses list:
ListOfAnalyses = SMSAnalysisFactory.load()
#Turn mass compression and invisible compression on and off
DoCompress = True
DoInvisible = True
minmassgap = addunit(5.,'GeV')  #Minimum mass gap for mass compression
sigmacut = addunit(0.1,'fb')  #Minimum cross-section*BR to be kept during decomposition
#Decompose and get the topologies/elements:
SMSTopList = LHEDecomposer.decompose(lhefile,None,nevts,DoCompress,DoInvisible,minmassgap)


#Print decomposition info:
EvTop_table = PrettyTable(["Topology","#Vertices", "#Insertions", "#Elements", "Sum of weights"])
EvElement_table = PrettyTable(["Topology","Element","Particles B[0]","Particles B[1]", "Masses B[0]","Masses B[1]","Element Weight"])
eltot = 0
totweight = []
for (i,topo) in enumerate(SMSTopList):
  weightlist = [el.weight for el in topo.ElList]
  sumw = ClusterTools.sumweights(weightlist)
  totweight.append(sumw)
  EvTop_table.add_row([i,topo.vertnumb,topo.vertparts,len(topo.ElList),wrap(printer.pformat(sumw),width=30)])
  eltot += len(topo.ElList)
  #Print element list for Topology[i]:  
  if i == 0:       
    for j in range(len(SMSTopList[i].ElList)):
      EvElement_table.add_row([i,j,SMSTopList[i].ElList[j].B[0].particles,SMSTopList[i].ElList[j].B[1].particles,wrap(printer.pformat(SMSTopList[i].ElList[j].B[0].masses),width=25),wrap(printer.pformat(SMSTopList[i].ElList[j].B[1].masses),width=25),wrap(printer.pformat(SMSTopList[i].ElList[j].weight),width=30)])
    EvElement_table.add_row(["---","---","---","---","---","---","---"])  
#Print Summary:      
print "Number of Global topologies = ",len(SMSTopList)      
print(EvTop_table)
print "Total Number of Elements = ",eltot
print "Total weight = ",ClusterTools.sumweights(totweight)
print(EvElement_table)

print '\n \n \n'


#Add event topologies to analyses:
for Analysis in ListOfAnalyses: Analysis.add(SMSTopList)
   
#Loop over analyses and compute theoretical predictions to each of them:
for Analysis in ListOfAnalyses:
  print "---------------Analysis Label = ",Analysis.label
  Results_table = PrettyTable(["Result","Conditions","Mass","Theoretical Value","Experimental Limits"])
  if max(Analysis.Top.vertnumb) > 3:
    print "Skipping analysis",Analysis.label,"with more than 3 vertices"
    continue

  th = Analysis.computeTheoryPredictions()  #Here we match the decomposition and the analysis
  if not th: continue
  elif th == 'Cluster Failed':
    print th
    continue
  
  res = Analysis.results.keys()[0]  

  if len(Analysis.ResultList) == 0: continue

#Print results for analysis:
  for cluster in Analysis.ResultList:
    mass = cluster.mass
    tvalue = cluster.result_dic.values()[0]
    conds = cluster.conditions_dic.values()[0]
    explimit = cluster.explimit
    sigmalimit = [Analysis.label,cluster.explimit]
    Results_table.add_row([wrap(printer.pformat(res),width=30),wrap(printer.pformat(conds),width=30),wrap(printer.pformat(mass),width=30),wrap(printer.pformat(tvalue),width=30),wrap(printer.pformat(sigmalimit),width=30)])
    Results_table.add_row(["---------","---------","---------","---------","---------"])

  print(Results_table)


sys.exit()
