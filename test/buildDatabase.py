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

for e in results[:]:
    print ( e.globalInfo.id )
    dsets = [ "SR1: MET > 200", "SR2: MET > 300" ]
    topo = "T2tt"
    effs = []
    for ds in dsets:
        eff = e.getEfficiencyFor ( topo, massvec, ds )
        if not eff: continue
        effs.append ( eff )
        ul = e.getUpperLimitFor ( ds )
        eul = e.getUpperLimitFor ( ds, expected=True )
        print ( "observed upper limit",ds," (on sigma*eff)", ul )
        print ( "expected upper limit",ds," (on sigma*eff)", eul )
        print ( "observed upper limit",ds," (on sigma)", ul/eff )
        print ( "expected upper limit",ds," (on sigma)", eul/eff )
        print ( "eff", eff )
    try:
        print ( "observed upper limit combined (on sigma)", e.getCombinedUpperLimitFor ( effs ) )
        print ( "expected upper limit combined (on sigma)", e.getCombinedUpperLimitFor ( effs, expected=True ) )
    except SModelSExperimentError as e:
        print ( "exception: %s" % e )
        
    tpred = theoryPredictionsFor ( e, smstoplist, useBestDataset=False )
    print ( "theorypredfor=%s" % tpred )
