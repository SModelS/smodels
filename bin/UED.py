#!/usr/bin/python

""" another sandbox to try things out """

import set_path, tempfile
from Theory import LHEDecomposer, XSecComputer
from Experiment import TxNames, SMSAnalysisFactory, SMSResults, LimitGetter
from Tools.PhysicsUnits import fb

print "[UED.py] loading analyses ...."
analyses = SMSAnalysisFactory.load( anas="alphaT8TeV", topos="T2" )
print "[UED.py] done loading %d analyses" % len(analyses)

nevts=1
lhefile="../lhe/ued_1.lhe"

weights= { '8 TeV (NLL)': { (-6100002, 6100002): .463 * fb } }

topos=LHEDecomposer.decompose ( lhefile, weights, nevts=nevts )

print
for Analysis in analyses:
  Analysis.add ( topos )
  lims=LimitGetter.limit ( Analysis )
  print "[run.py] lims=",lims
