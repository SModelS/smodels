#!/usr/bin/python

""" A first UED example """

import set_path, tempfile
from Theory import LHEDecomposer, XSecComputer
from Experiment import TxNames, SMSAnalysisFactory, SMSResults, LimitGetter, SMSHelpers
from Tools.PhysicsUnits import fb
from Tools import FeynmanGraphs
from Tools.VariousHelpers import logging

## SMSHelpers.Base="/home/lisa/daten/smodels"

print "[UED.py] loading analyses ...."
## analyses = SMSAnalysisFactory.load( anas="alphaT8TeV", topos="T2" )
analyses = SMSAnalysisFactory.load()
#analyses = SMSAnalysisFactory.load( topos="T2" )

print "[UED.py] done loading %d analyses" % len(analyses)

nevts=10
lhefile="../lhe/ued_2.lhe" ## thats the lhe file we're using
#lhefile="../lhe/firstry.lhe" ## thats the lhe file we're using


## we create the correct object manually
Wv=XSecComputer.loFromLHE ( lhefile, totalxsec = None ) 
weights={ '8 TeV (LO)': Wv[0] }
print "[UED.py] weights=",Wv[0]

## now create the list of topologies
topos=LHEDecomposer.decompose ( lhefile, weights, nevts=nevts )
for topo in topos: FeynmanGraphs.asciidraw ( topo.leadingElement() )

for topo in topos:
  for element in topo.elements():
    FeynmanGraphs.asciidraw( element )
    ## FeynmanGraphs.draw( element )

print
for Analysis in analyses:
  Analysis.add ( topos ) ## add all the theory information to the analysis objcts
  lims=LimitGetter.limit ( Analysis, addTheoryPredictions=[ '8 TeV (LO)' ] ) ## obtain the relevant limits
  if len(lims): print "[UED.py] lims=",lims

