#!/usr/bin/python

from __future__ import print_function
import sys, os, time
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.exceptions import SModelSExperimentError
from smodels.tools.smodelsLogging import setLogLevel
from smodels.tools.colors import colors
from smodels.tools.physicsUnits import pb, fb, GeV
from smodels.theory import slhaDecomposer

colors.on = True
# setLogLevel ( "debug" )

smstoplist = smstoplist = slhaDecomposer.decompose( "T2tt.slha" )
# print ( "smstoplist=",len(smstoplist ) )
dir = "corrdb/"
d=Database( dir, discard_zeroes = True )
print(d)
results=d.getExpResults()

massvec = [[400.*GeV,75*GeV], [400*GeV,75*GeV]]

for e in results:
    print ( e.globalInfo.id )
    e0 = e
    sig1,sig2=1.0*fb,1.0*fb
    y1,y2=float(sig1*e0.globalInfo.lumi),float(sig2*e0.globalInfo.lumi)
    eff1 = e0.getEfficiencyFor ( "T2tt", massvec, "SR1: MET > 200" )
    eff2 = e0.getEfficiencyFor ( "T2tt", massvec, "SR2: MET > 300" )
    sig1,sig2=1.,1.
    print ( "yields",y1,y2 )
    print ( "signals",sig1, sig2 )
    print ( "signals",eff1, eff2 )
    try:
        print ( "upper limit combined", e0.getCombinedUpperLimitFor ( [ eff1 , eff2 ] ) )
    except SModelSExperimentError as e:
        print ( "exception: %s" % e )
        
    print ( "upper limit SR1", e0.getUpperLimitFor ( "SR1: MET > 200" )  )
    print ( "eff1", e0.getEfficiencyFor ( "T2tt", massvec, "SR1: MET > 200" ) )
    print ( "upper limit SR2", e0.getUpperLimitFor ( "SR2: MET > 300" )  )
    print ( "eff2", e0.getEfficiencyFor ( "T2tt", massvec, "SR2: MET > 300" ) )

    tpred = theoryPredictionsFor ( e0, smstoplist, useBestDataset=False )
    print ( "theorypredfor=%s" % tpred )
