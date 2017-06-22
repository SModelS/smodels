#!/usr/bin/python

import sys, os, time
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.tools.smodelsLogging import setLogLevel
from smodels.tools.colors import colors
from smodels.tools.physicsUnits import pb, fb
from smodels.theory import slhaDecomposer

colors.on = True
# setLogLevel ( "debug" )

smstoplist = smstoplist = slhaDecomposer.decompose( "T2tt.slha" )
# print ( "smstoplist=",len(smstoplist ) )
dir = "corrdb/"
d=Database( dir, discard_zeroes = True )
print(d)
e=d.getExpResults()
e0=e[0]
print (e)
sig1,sig2=1.0*fb,1.0*fb
y1,y2=float(sig1*e0.globalInfo.lumi),float(sig2*e0.globalInfo.lumi)
print ( "yields",y1,y2 )
print ( "upper limit", e0.getCombinedUpperLimitFor ( [ sig1 , sig2 ] ) )
print ( "upper limit", e0.getUpperLimitFor ( "SR1: MET > 200" ) / y1 )
print ( "upper limit", e0.getUpperLimitFor ( "SR2: MET > 300" ) / y2 )

tpred = theoryPredictionsFor ( e0, smstoplist )
print ( "theorypredfor=%s" % tpred )
