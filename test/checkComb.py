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
setLogLevel ( "debug" )

smstoplist = smstoplist = slhaDecomposer.decompose( "T2tt.slha" )
print ( "smstoplist=",len(smstoplist ) )
dir = "covdb2/"
d=Database( dir, discard_zeroes = True )
print(d)

massvec = [[400.*GeV,75*GeV], [400*GeV,75*GeV]]

for e in d.getExpResults():
    print ( e.globalInfo.id )
    #for ds in e.datasets:
    #    print ( ds.dataInfo.dataId )
    predictions = theoryPredictionsFor ( e, smstoplist, useBestDataset=False, 
                                         combinedResults=True )
    for pred in predictions:
        print ( pred.describe() )

