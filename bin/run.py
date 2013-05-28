#!/usr/bin/python

""" another sandbox to try things out """

import set_path, argparse
from Theory import LHEReader, SMSDecomp, SMSmethods, TopologyBuilder, SMSDataObjects
from Experiment import TxNames
import Experiment.SMSAnalysisFactory as Analyses
## import SMSanalyses as Analyses
## import SMSanalysesTest as Analyses


print "[run.py] loading analyses ...."
ListOfAnalyses = Analyses.load( topos="T2" )
print "[run.py] done loading %d analyses" % len(ListOfAnalyses)

# n=1000
# slhafile="AndreSLHA/andrePT4.slha"
lhefile="../lhe/T2_1.lhe"
# lhefile="T2_100.lhe"
reader=LHEReader.LHEReader ( lhefile, 1 )
SMSTopList=[]
for Event in reader:
  print Event
  SMSTops=TopologyBuilder.fromEvent ( Event )
  print "SMSTops[0].leadingElement=",SMSTops[0].leadingElement()

  for SMSTop in SMSTops:
    SMSTop.checkConsistency(verbose=True)
    print "SMSTop.leadingElement=",SMSTop.leadingElement()
    SMSmethods.AddToList ( SMSTop, SMSTopList )
    #  element=SMSTop.leadingElement()
    for element in SMSTop.elements():
      print "======================="
      tx=TxNames.getTx(element)
      print "element=",element
      print "B0=",element.B[0]
      print "tx=",tx

print
print
for SMSTop in SMSTopList: 
  print "SMSTopList SMS topology",SMSTop
  element=SMSTop.leadingElement() 
  print "SMSTopList leading element",element
  tx=TxNames.getTx(element)
  print "SMSTopList Tx=",tx
  
for Analysis in ListOfAnalyses:
  Analysis.evaluateResults()
  Analysis.add ( SMSTopList )
  print "------------- Analysis Label = ",Analysis.label, Analysis.sqrts, Analysis.Top
  print "    `-- element ",Analysis.Top.leadingElement()
  print "    `-- ",Analysis.Top.leadingElement().B[0]
  d="[[[jet]],[[jet]]]"
  #d="[[['jet','jet'],['jet']],[['jet']]]"
  # S=SMSmethods.strtoel ( d )
  e=SMSDataObjects.EElement ( d )
  # print "S=",S,"e=",e,Analysis.Top.leadingElement()==e
  #print "              S=``%s'' -> ``%s''" % ( S, Analysis.Top.leadingElement().B[0] )

