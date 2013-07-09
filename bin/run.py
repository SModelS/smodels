#!/usr/bin/python

""" another sandbox to try things out """

import set_path, tempfile
from Theory import LHEDecomposer, XSecComputer
from Experiment import TxNames, SMSAnalysisFactory, SMSResults, LimitGetter

print "[run.py] loading analyses ...."
# analyses = SMSAnalysisFactory.load( )
analyses = SMSAnalysisFactory.load( anas="alphaT8TeV", topos="T2bb" )
print "[run.py] done loading %d analyses" % len(analyses)

nevts=1000
# slhafile="../slha/TChiNuSlep.slha"
slhafile="../slha/T2bb.slha"

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
  Ret=Analysis.computeTheoryPredictions()
  # lims=LimitGetter.limit ( Analysis )
  #lims=[]
  #if len(lims)==0: continue
  print "[run.py] -------------------------------------------------"
  print "[run.py] analysis=",Analysis ## ,"lims=",lims
  print "[run.py] plots=",Analysis.plots
  print "[run.py] v=",Analysis.plots.values()[0][1][0]
  #for lim in lims:
  #  print "[run.py] limit=",lim
  if len(Analysis.ResultList)==0: continue
  for cluster in Analysis.ResultList:
    print "[run.py] cluster",cluster
    theoRes=cluster.oldformat()
    print "[run.py] cluster",theoRes
