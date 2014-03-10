#!/usr/bin/python

import sys
from prettytable import PrettyTable
from theory import LHEDecomposer, SLHADecomposer, XSecComputer, ClusterTools, CrossSection, SLHATools, TheoryPrediction
from tools.PhysicsUnits import addunit, rmvunit
from tools import SMSPrettyPrinter, VariousHelpers
from tools.SMSPrettyPrinter import wrap
from tools.VariousHelpers import logging
from experiment import SMSAnalysisList, SMSAnalysisFactory, LimitGetter, SMSHelpers


printer=SMSPrettyPrinter.SMSPrettyPrinter()

#Example of how to define cross-sections. If not defined, default values
#will be generated and stored in CrossSection.XSectionInfo the first
#time this information is needed 
#XSectionInfo = CrossSection.XSecInfoList('8 TeV (NLL), 8 TeV (nlo)')
#CrossSection.XSectionInfo = XSectionInfo

#Generate events and compute cross-sections:
nevts = 10000
slhafile = "slha/andrePT4.slha"
#slhafile = "slha/DESY_stop.slha"

WriteToFile = True
if not WriteToFile:
  Wv = XSecComputer.compute(nevts,slhafile,rpythia = True)
  W = Wv["Wdic"]
  Xsec = Wv["Xsecdic"]
  lhefile = Wv["lhefile"]
else:
  Xsec=None
  SLHATools.writeXSecToSLHAFile(slhafile,nevts,printLHE=False)

#PYTHIA must have MSTP(42)=0 ! no mass smearing (narrow width approximation)
#Creat analyses list:
DoFactory = True
if DoFactory:
  ListOfAnalyses = SMSAnalysisFactory.load()
else:  
  ListOfAnalyses = SMSAnalysisList.load()




DoCompress = True
DoInvisible = True
minmassgap = addunit(10.,'GeV')

DoSLHAdec = True
if DoSLHAdec:
  sigmacut = addunit(0.1,'fb')
  if DoCompress or DoInvisible: sigmacut = sigmacut/10.  #When compression is turned on, relax sigmacut
  SMSTopList = SLHADecomposer.decompose(slhafile,sigmacut,DoCompress,DoInvisible,minmassgap)
else:
  SMSTopList = LHEDecomposer.decompose(lhefile,W,nevts,DoCompress,DoInvisible,minmassgap)


EvTop_table = PrettyTable(["Topology","#Vertices", "#Insertions", "#Elements", "Sum of weights"])
EvElement_table = PrettyTable(["Topology","Element","Particles B[0]","Particles B[1]", "Masses B[0]","Masses B[1]","Element Weight"])

eltot = 0
totweight = []
#Print Results:
# for i in range(len(SMSTopList)):
for (i,topo) in enumerate(SMSTopList):
  weightlist = [el.weight for el in topo.elList]
  sumw = ClusterTools.sumweights(weightlist)
  totweight.append(sumw)
  EvTop_table.add_row([i,topo.vertnumb,topo.vertparts,len(topo.elList),wrap(printer.pformat(sumw),width=30)])
  eltot += len(topo.elList)

 
      
#Print element list for Topology[i]:  
  if i >= 0:       
    for j in range(len(SMSTopList[i].elList)):
      EvElement_table.add_row([i,j,SMSTopList[i].elList[j].B[0].particles,SMSTopList[i].elList[j].B[1].particles,wrap(printer.pformat(SMSTopList[i].elList[j].B[0].masses),width=25),wrap(printer.pformat(SMSTopList[i].elList[j].B[1].masses),width=25),wrap(printer.pformat(SMSTopList[i].elList[j].weight),width=30)])
    EvElement_table.add_row(["---","---","---","---","---","---","---"])  
      
      
print "Number of Global topologies = ",len(SMSTopList)      
print(EvTop_table)
print "Total Number of Elements = ",eltot
print "Total weight = ",ClusterTools.sumweights(totweight)
#print(EvElement_table)

print '\n \n \n'


#sys.exit()


#print '\n \n \n'
#sys.exit()

    
#Compute theoretical predictions to analyses results:
for Analysis in ListOfAnalyses:
  print "---------------Analysis Label = ",Analysis.label
  Results_table = PrettyTable(["Result","Conditions","Mass","Theoretical Value","Experimental Limits"])
  
  
  if max(Analysis.Top.vertnumb) > 3:
    print "Skipping analysis",Analysis.label,"with more than 3 vertices"
    continue

  result = TheoryPrediction.TheoryPrediction()
  if not result.computeTheoryPredictions(Analysis,SMSTopList): continue #Theoretical values for result and conditions

  if len(result.clusterResults) == 0: continue

  for cluster in result.clusterResults:
    mass = cluster.mass
    res = cluster.result_dic.keys()[0]
    tvalue = cluster.result_dic.values()
    conds = cluster.conditions_dic.values()
    sigmalimit = [result.analysis.label,cluster.explimit]

    Results_table.add_row([wrap(printer.pformat(res),width=30),wrap(printer.pformat(conds),width=30),wrap(printer.pformat(mass),width=30),wrap(printer.pformat(tvalue),width=30),wrap(printer.pformat(sigmalimit),width=30)])

    Results_table.add_row(["---------","---------","---------","---------","---------"])

  print(Results_table)


sys.exit()
