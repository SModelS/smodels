#!/usr/bin/python

""" another sandbox to try things out """

import set_path, argparse
from Theory import LHEDecomposer, SMSDataObjects, XSecComputer
from Experiment import TxNames
import Experiment.SMSAnalysisFactory as Analyses
from Tools import SMSFeynmanGraphs

print "[run.py] loading analyses ...."
ListOfAnalyses = Analyses.load( topos="T2" )
print "[run.py] done loading %d analyses" % len(ListOfAnalyses)

nevts=1000
# slhafile="AndreSLHA/andrePT4.slha"
slhafile="../slha/T2.slha"
Wv=XSecComputer.compute(nevts,slhafile,rpythia = True, donlo = True)
print "Wdic=",Wv["Wdic"]
print Wv["Xsecdic"]
lhefile="../lhe/T2_1.lhe"
# lhefile="T2_100.lhe"
topos=LHEDecomposer.decompose ( lhefile, Wv["Wdic"], None )

for topo in topos:
    topo.checkConsistency(verbose=True)
    print "topo.leadingElement=",topo.leadingElement()
    for element in topo.elements():
      print "======================="
      tx=TxNames.getTx(element)
      print "element=",element
      print "B0=",element.B[0]
      print "tx=",tx

print
print
for topo in topos.topos:
  print "topos SMS topology",topo
  element=topo.leadingElement() 
  print "topos leading element",element
  tx=TxNames.getTx(element)
  print "topos Tx=",tx
  
for Analysis in ListOfAnalyses:
  Analysis.evaluateResults()
  Analysis.add ( topos )
  print "------------- Analysis Label = ",Analysis.label, Analysis.sqrts, Analysis.Top
  print "    `-- element ",Analysis.Top.leadingElement()
  print "    `-- ",Analysis.Top.leadingElement().B[0]
  d="[[[jet]],[[jet]]]"
  #d="[[['jet','jet'],['jet']],[['jet']]]"
  # S=SMSmethods.strtoel ( d )
  e=SMSDataObjects.EElement ( d )
  SMSFeynmanGraphs.asciidraw ( e )
  # print "S=",S,"e=",e,Analysis.Top.leadingElement()==e
  #print "              S=``%s'' -> ``%s''" % ( S, Analysis.Top.leadingElement().B[0] )

