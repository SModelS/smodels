#!/usr/bin/python

import sys
from prettytable import PrettyTable
from Theory import LHEDecomposer, SLHADecomposer, XSecComputer, ClusterTools, CrossSection
from Tools.PhysicsUnits import addunit, rmvunit
from Tools import SMSPrettyPrinter, VariousHelpers
from Tools.SMSPrettyPrinter import wrap
from Tools.VariousHelpers import logging
from Experiment import SMSAnalysisList, SMSAnalysisFactory, LimitGetter


printer=SMSPrettyPrinter.SMSPrettyPrinter()

#Example of how to define cross-sections. If not defined, default values
#will be generated and stored in CrossSection.XSectionInfo the first
#time this information is needed 
#XSectionInfo = CrossSection.XSecInfoList('8 TeV (NLL), 8 TeV (nlo)')
#CrossSection.XSectionInfo = XSectionInfo

#Generate events and compute cross-sections:
nevts = 10000
#slhafile = "slha/andrePT4.slha"
slhafile = "slha/DESY_stop.slha"
Wv = XSecComputer.compute(nevts,slhafile,rpythia = True)
W = Wv["Wdic"]
Xsec = Wv["Xsecdic"]
lhefile = Wv["lhefile"]

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
  maxlum = VariousHelpers.getMaxLum(ListOfAnalyses) # Maximum cross-section*BR to be included
  if rmvunit(maxlum,'fb-1'):    
    sigmacut = addunit(1./rmvunit(maxlum,'fb-1'),'fb')
  else:
    sigmacut = addunit(0.1,'fb')
  if DoCompress or DoInvisible: sigmacut = sigmacut/10.  #When compression is turned on, relax sigmacut
  SMSTopList = SLHADecomposer.decompose(slhafile,Xsec,sigmacut,DoCompress,DoInvisible,minmassgap)
else:
  SMSTopList = LHEDecomposer.decompose(lhefile,W,nevts,DoCompress,DoInvisible,minmassgap)


EvTop_table = PrettyTable(["Topology","#Vertices", "#Insertions", "#Elements", "Sum of weights"])
EvElement_table = PrettyTable(["Topology","Element","Particles B[0]","Particles B[1]", "Masses B[0]","Masses B[1]","Element Weight"])

eltot = 0
totweight = []
#Print Results:
# for i in range(len(SMSTopList)):
for (i,topo) in enumerate(SMSTopList):
  weightlist = [el.weight for el in topo.ElList]
  sumw = ClusterTools.sumweights(weightlist)
  totweight.append(sumw)
  EvTop_table.add_row([i,topo.vertnumb,topo.vertparts,len(topo.ElList),wrap(printer.pformat(sumw),width=30)])
  eltot += len(topo.ElList)

 
      
#Print element list for Topology[i]:  
  if i >= 0:       
    for j in range(len(SMSTopList[i].ElList)):
      EvElement_table.add_row([i,j,SMSTopList[i].ElList[j].B[0].particles,SMSTopList[i].ElList[j].B[1].particles,wrap(printer.pformat(SMSTopList[i].ElList[j].B[0].masses),width=25),wrap(printer.pformat(SMSTopList[i].ElList[j].B[1].masses),width=25),wrap(printer.pformat(SMSTopList[i].ElList[j].weight),width=30)])
    EvElement_table.add_row(["---","---","---","---","---","---","---"])  
      
      
print "Number of Global topologies = ",len(SMSTopList)      
print(EvTop_table)
print "Total Number of Elements = ",eltot
print "Total weight = ",ClusterTools.sumweights(totweight)
#print(EvElement_table)

print '\n \n \n'


#sys.exit()
#Add event topologies to analyses:
for Analysis in ListOfAnalyses:
  Analysis.add(SMSTopList)


#Print analyses output:
AnElement_table = PrettyTable(["Analyses","Element","Masses","Element Weight"])  



for Ana in ListOfAnalyses:
  label = Ana.label
  ifirst = True
  for iel,El in enumerate(Ana.Top.ElList):
    ptcs = El.ParticleStr
    for im,massweight in enumerate(El.MassWeightList):
      mass = massweight.mass
      if not ifirst: label = ""
      if im != 0: ptcs = ""
      ifirst = False
      AnElement_table.add_row([label,ptcs,wrap(printer.pformat(mass),width=100),wrap(printer.pformat(massweight.weight),width=30)])
  AnElement_table.add_row(["---","---","---","---"])  
    

#print(AnElement_table)


#print '\n \n \n'
#sys.exit()

    
#Compute theoretical predictions to analyses results:
for Analysis in ListOfAnalyses:
  print "---------------Analysis Label = ",Analysis.label
  Results_table = PrettyTable(["Result","Conditions","Mass","Theoretical Value","Experimental Limits"])
  
  
  if max(Analysis.Top.vertnumb) > 3:
    print "Skipping analysis",Analysis.label,"with more than 3 vertices"
    continue


  if not Analysis.computeTheoryPredictions(): continue #Theoretical values for result and conditions

  res = Analysis.results.keys()[0]  

  if len(Analysis.ResultList) == 0: continue

  for cluster in Analysis.ResultList:
    theoRes = cluster.oldformat()   #Convert to old format (except for small change in conditions output)
    mass = theoRes['mass']
    tvalue = theoRes['result']
    conds = theoRes['conditions']
    sigmalimit = [Analysis.plots.values()[0][1][0],cluster.explimit]

    Results_table.add_row([wrap(printer.pformat(res),width=30),wrap(printer.pformat(conds),width=30),wrap(printer.pformat(mass),width=30),wrap(printer.pformat(tvalue),width=30),wrap(printer.pformat(sigmalimit),width=30)])

    Results_table.add_row(["---------","---------","---------","---------","---------"])

  print(Results_table)


sys.exit()
