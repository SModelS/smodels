#!/usr/bin/python

""" A first UED example """

import set_path, tempfile
from Theory import LHEDecomposer, XSecComputer
from Experiment import TxNames, SMSAnalysisFactory, SMSResults, LimitGetter
from Tools.PhysicsUnits import fb

print "[UED.py] loading analyses ...."
analyses = SMSAnalysisFactory.load( anas="alphaT8TeV", topos="T2" )
print "[UED.py] done loading %d analyses" % len(analyses)

nevts=1
lhefile="../lhe/ued_1.lhe" ## thats the lhe file we're using

## we supply the weight manually
weights= { '8 TeV (NLL)': { (-6100002, 6100002): .463 * fb } }

## now create the list of topologies
topos=LHEDecomposer.decompose ( lhefile, weights, nevts=nevts )

print
for Analysis in analyses:
  Analysis.add ( topos ) ## add all the theory information to the analysis objcts
  lims=LimitGetter.limit ( Analysis ) ## obtain the relevant limits
  print "[UED.py] lims=",lims
