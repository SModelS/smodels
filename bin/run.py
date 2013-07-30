#!/usr/bin/python

""" another sandbox to try things out """

import set_path, tempfile
from Theory import LHEDecomposer, XSecComputer, ClusterTools
from Experiment import TxNames, SMSAnalysisFactory, SMSResults, LimitGetter
from Tools.VariousHelpers import logging

log = logging.getLogger(__name__)

nevts=1000
# slhafile="../slha/TChiNuSlep.slha"
slhafile="../slha/T2bb.slha"

Tmp=tempfile.mkdtemp()
log.info ( "now run pythia in "+Tmp )
Wv=XSecComputer.compute(nevts,slhafile,rpythia = True, datadir=Tmp)
log.info ( "done running pythia" )
# print "cmdic=",ClusterTools.CMdic

lhefile=Wv.lhefile ( sqrts=8 )

## print "[run.py] weights=",Wv.weights()
topos=LHEDecomposer.decompose ( lhefile, Wv.weights(), nevts=nevts )
XSecComputer.clean ( Tmp )

print "topos=",topos

log.info ( "loading analyses" )
# analyses = SMSAnalysisFactory.load( )
analyses = SMSAnalysisFactory.load( anas="alphaT8TeV", topos="T2bb" )
log.info ( "done loading %d analyses" % len(analyses) )

print
for Analysis in analyses:
  Analysis.add ( topos )
  Ret=Analysis.computeTheoryPredictions()
  # lims=LimitGetter.limit ( Analysis )
  #lims=[]
  #if len(lims)==0: continue
  log.info ( "-------------------------------------------------" )
  log.info ( "analysis="+str(Analysis) ) 
  log.info ( "plots="+str(Analysis.plots) )
  log.info ( "v="+str(Analysis.plots.values()[0][1][0]) )
  #for lim in lims:
  #  print "[run.py] limit=",lim
  if len(Analysis.ResultList)==0: continue
  for cluster in Analysis.ResultList:
    log.info ( "cluster"+str(cluster) )
    log.info ( "cluster class %s" % str ( cluster.__class__ ) )
    theoRes=cluster.oldformat()
    log.info ( "cluster oldformat"+str(theoRes) )
