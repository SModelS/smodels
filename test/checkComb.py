#!/usr/bin/python3

from __future__ import print_function
import sys, os, time
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.exceptions import SModelSExperimentError
from smodels.tools.smodelsLogging import setLogLevel
from smodels.tools.colors import colors
from smodels.tools.physicsUnits import pb, fb, GeV
from smodels.tools.likelihoods import LikelihoodComputer
from smodels.tools.statistics import UpperLimitComputer
from smodels.tools import ioObjects
from smodels.tools import printer
from smodels.theory import slhaDecomposer

LikelihoodComputer.deltas_default = 1e-12
LikelihoodComputer.debug_mode = True
UpperLimitComputer.debug_mode = True

prt = printer.PyPrinter( output="file", filename = "./out.py" )

colors.on = True
setLogLevel ( "debug" )
# setLogLevel ( "info" )

file = "T2tt_533_40_533_40.slha"
file = "T2tt.slha"

smstoplist = smstoplist = slhaDecomposer.decompose( file )
print ( "smstoplist=",len(smstoplist ) )
dir = "covdb/"

if len ( sys.argv) > 1:
    dir = sys.argv[1]

d=Database( dir, discard_zeroes = True )
print(d)

expRes = d.getExpResults( datasetIDs=['all'] )

print ( "%d results." % len(expRes) )

for e in expRes:
    print ( e.globalInfo.id )
    #for ds in e.datasets:
    #    print ( ds.dataInfo.dataId )
    predictions = theoryPredictionsFor ( e, smstoplist, useBestDataset=False, 
                                         combinedResults=True )
    results = ioObjects.ResultList(predictions )
    prt.addObj ( results )
    for pred in predictions[:]:
        print ( pred.describe() )

prt.flush()
