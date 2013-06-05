#!/usr/bin/python

""" another sandbox to try things out """

import set_path, tempfile
from Theory import LHEDecomposer, XSecComputer
from Experiment import TxNames, SMSAnalysisFactory, SMSResults, LimitGetter
from Tools.PhysicsUnits import fb

print "[run.py] loading analyses ...."
analyses = SMSAnalysisFactory.load( anas="alphaT8TeV", topos="T2" )
print "[run.py] done loading %d analyses" % len(analyses)

nevts=1
#slhafile="../slha/T2bb.slha"

#Tmp=tempfile.mkdtemp()
#print "[run.py] now run pythia in",Tmp
#Wv=XSecComputer.compute(nevts,slhafile,rpythia = True, donlo = True, datadir=Tmp)
#print "[run.py] done running pythia"

#lhefile=Wv.lhefile ( 8 )
lhefile="../lhe/unweighted_events_sms.lhe"

#weights= {'7 TeV (LO)': {(-1000005, 1000005): 1.15E-01 [fb], (1000022, 1000022): 4.63E-01 [fb]}, '8 TeV (LO)': {(-1000005, 1000005): 1.90E-01 [fb], (1000022, 1000022): 1.90E-01 [fb]}, '7 TeV (NLL)': {(1000022, 1000022): 4.63E-01 [fb], (-1000005, 1000005): 2.05E-01 [fb]}, '8 TeV (NLL)': {(1000022, 1000022): 1.90E-01 [fb], (-1000005, 1000005): 3.27E-01 [fb]}}
weights= { '8 TeV (NLL)': { (-6100002, 6100002): .463 * fb } }

topos=LHEDecomposer.decompose ( lhefile, weights, nevts=nevts )
#XSecComputer.clean ( Tmp )

print
for Analysis in analyses:
  Analysis.add ( topos )
  lims=LimitGetter.limit ( Analysis )
  print "[run.py] lims=",lims
