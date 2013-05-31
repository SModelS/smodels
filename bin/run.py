#!/usr/bin/python

""" another sandbox to try things out """

import set_path, tempfile
from Theory import LHEDecomposer, XSecComputer
from Experiment import TxNames, SMSAnalysisFactory, SMSResults

print "[run.py] loading analyses ...."
analyses = SMSAnalysisFactory.load( anas="alphaT8TeV", topos="T2bb" )
print "[run.py] done loading %d analyses" % len(analyses)

nevts=1000
slhafile="../slha/T2bb.slha"

Tmp=tempfile.mkdtemp()
print "[run.py] now run pythia in",Tmp
Wv=XSecComputer.compute(nevts,slhafile,rpythia = True, donlo = True, datadir=Tmp)
print "[run.py] done running pythia"

lhefile=Wv.lhefile ( 8 )
topos=LHEDecomposer.decompose ( lhefile, Wv.weights(), nevts=nevts )
XSecComputer.clean ( Tmp )

print
for Analysis in analyses:
  Analysis.add ( topos )
  print "========================================="
  print "[run.py] analysis ",Analysis
  lead=Analysis.Top.leadingElement()
  print "[run.py]  element ",lead
  print "[run.py]   masses ",lead.B[0].masses,lead.B[1].masses
  print "[run.py]    plots ",Analysis.plots
  for (constraint,condition) in Analysis.results.items():
    print "[run.py] -- constraint=",constraint,"condition=",condition
    theoRes=Analysis.evaluateResult( constraint )
    print "[run.py] -- plot=",Analysis.plots[constraint]
    print "[run.py] -- len theoRes=",len(theoRes)
    print "[run.py] --    result 0=",theoRes[0]
    Tx0=Analysis.plots[constraint][0]
    ana0=Analysis.plots[constraint][1][0]
    masses=lead.B[0].masses[0] ## ,lead.B[1].masses[0]
    print "[run.py] --    plot 0=",Tx0,ana0,masses
    ul=SMSResults.getSmartUpperLimit(ana0,Tx0,masses)
    refxsec=theoRes[0]['result']['8 TeV (NLL)']
    excluded=refxsec>ul
    print "[run.py] upper limit=",ul,"theory prediction=",refxsec,excluded
