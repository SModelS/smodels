#!/usr/bin/python

""" another sandbox to try things out """

import set_path, tempfile
from Theory import LHEDecomposer, XSecComputer
from Experiment import TxNames, SMSAnalysisFactory, SMSResults, LimitGetter

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
  lims=LimitGetter.limit ( Analysis )
  print "lims=",lims
  for (constraint,condition) in Analysis.results.items():
    Tx0=Analysis.plots[constraint][0]
    ana0=Analysis.plots[constraint][1][0]
    theoRes=Analysis.evaluateResult( constraint )
    refxsec=theoRes[0]['result']['8 TeV (NLL)']
    excluded=refxsec>lims[ana0+Tx0]['ul']
    print "[run.py] upper limit=",lims[ana0+Tx0],"theory prediction=",refxsec,excluded
