#!/usr/bin/python

""" another sandbox to try things out """

import set_path, tempfile
from Theory import LHEDecomposer, XSecComputer
from Experiment import TxNames, SMSAnalysisFactory, SMSResults, LimitGetter

print "[run.py] loading analyses ...."
## analyses = SMSAnalysisFactory.load( anas="alphaT8TeV", topos="T2bb" )
analyses = SMSAnalysisFactory.load( )
print "[run.py] done loading %d analyses" % len(analyses)

nevts=1000
slhafile="../slha/TChiNuSlep.slha"

Tmp=tempfile.mkdtemp()
print "[run.py] now run pythia in",Tmp
Wv=XSecComputer.compute(nevts,slhafile,rpythia = True, donlo = True, datadir=Tmp)
print "[run.py] done running pythia"

lhefile=Wv.lhefile ( 8 )

## print "[run.py] weights=",Wv.weights()
topos=LHEDecomposer.decompose ( lhefile, Wv.weights(), nevts=nevts )
XSecComputer.clean ( Tmp )

print
for Analysis in analyses:
  Analysis.add ( topos )
  lims=LimitGetter.limit ( Analysis )
  if len(lims)==0: continue
  print "[run.py] -------------------------------------------------"
  print "[run.py] analysis=",Analysis,"lims=",lims
  Analysis.evaluateResult()
  if len(Analysis.ResultList)==0: continue
  for cluster in Analysis.ResultList:
    theoRes=cluster.result_dic
    print "[run.py] ****"
    print "[run.py] theoRes=",theoRes
