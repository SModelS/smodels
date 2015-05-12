#!/usr/bin/python

import argparse, types
argparser = argparse.ArgumentParser(description="validates one pair of txname / axes")
argparser.add_argument ( '-T', '--txname', nargs='?', help='txname [default T1]',
                type=types.StringType, default='T1' )
argparser.add_argument ( '-a', '--axes', nargs='?', help='axes description [default 2*Eq(mother,x)_Eq(lsp,y)]',
                type=types.StringType, default='2*Eq(mother,x)_Eq(lsp,y)' )
args=argparser.parse_args() 


import sys,os
sys.path.insert(0,"../../../../../smodels-utils/")
sys.path.insert(0,"../../../../../smodels/")

from smodels.experiment.databaseObjects import DataBase
from validation.plotProducer import ValidationPlot, getExpIdFromPath
from smodels.tools.physicsUnits import pb

## ValidationPlot.computeWrongnessFactor = computeWrongnessFactor

filename="%s_%s.py" % ( args.txname, args.axes.replace("(","").replace(")","").replace(",","").replace("*","") )
execfile(filename)

database = DataBase("../../../../")
expRes = database.getExpResults(analysisIDs=[ getExpIdFromPath() ],datasetIDs=[None])


plot=ValidationPlot( expRes, args.txname, args.axes )
plot.data=validationData
agreement = plot.computeAgreementFactor()
print "agreement=",agreement

plot.getPlot()
plot.savePlot()

# import IPython
# IPython.embed()
