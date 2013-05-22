#!/usr/bin/python

import sys, argparse
sys.path.append ( "../" )
import SMSglobals, SMSmethods, TxNames, LHEReader
import LHEReader
## import SMSAnalysisFactory as Analyses
## import SMSanalyses as Analyses
import SMSanalysesTest as Analyses
import SMSDecomp

SMSglobals.initglob()
print "[run.py] loading analyses ...."
SMSglobals.ListOfAnalyses = Analyses.load()
print "[run.py] done loading %d analyses" % len(SMSglobals.ListOfAnalyses)

# n=1000
# slhafile="AndreSLHA/andrePT4.slha"
lhefile="../regression/T2_1.lhe"
reader=LHEReader.LHEReader ( lhefile )
for Event in reader:
  print Event
  SMSTop=SMSmethods.getEventTop ( Event.particles )
  element=SMSTop[0].leadingElement()
  tx=TxNames.getTx(element)
  print "element=",element,
  print "tx=",tx

for Analysis in SMSglobals.ListOfAnalyses:
  print "---------------Analysis Label = ",Analysis.label   
  for res in Analysis.results.keys():
    theoRes = SMSmethods.EvalRes(res,Analysis)
    print "[SMSmain.py] res=",res,"theoRes len=",theoRes
