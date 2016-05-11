#!/usr/bin/python

import argparse, types
argparser = argparse.ArgumentParser(description="validates one pair of txname/axes")
argparser.add_argument ( '-T', '--txname', nargs='?', 
         help='txname [default T1]', type=types.StringType, default='T1' )
argparser.add_argument ( '-p', '--plottype', nargs='?', help=
         'plottype [plain, bestregion, upperlimits, crosssections. default plain]',
         type=types.StringType, default='plain' )
argparser.add_argument ( '-a', '--axes', nargs='?', 
         help='axes description [default 2*Eq(mother,x)_Eq(lsp,y)]',
         type=types.StringType, default='2*Eq(mother,x)_Eq(lsp,y)' )
argparser.add_argument ( '-n', '--nthpoint', nargs='?', 
         help='plot only every nth point [1]',
         type=types.IntType, default=1 )
argparser.add_argument ( '-s', '--signal_factor', nargs='?', 
         help='factor to multiply theory cross section with, before comparing [1.0]',
         type=types.FloatType, default=1.0 )

args=argparser.parse_args() 


import sys,os
home = os.path.expanduser("~")
sys.path.insert(0,os.path.join(home,"smodels-utils/"))
sys.path.insert(0,os.path.join(home,"smodels/"))

from smodels.experiment.databaseObjects import Database
from validation.plotProducer import ValidationPlot, getExpIdFromPath
from validation.plotProducer import getDatasetIdsFromPath
from smodels.tools.physicsUnits import pb, fb
NAN=float('nan')

## ValidationPlot.computeWrongnessFactor = computeWrongnessFactor

filename="%s_%s.py" % ( args.txname, args.axes.replace("(","").replace(")","").replace(",","").replace("*","") )
execfile(filename)

database = Database(os.path.join(home,"smodels-database"))
expRes = database.getExpResults(analysisIDs=[ getExpIdFromPath() ],datasetIDs=getDatasetIdsFromPath() )

for res in expRes:
    plot=ValidationPlot( res, args.txname, args.axes )
    plot.data=validationData
    agreement = plot.computeAgreementFactor( signal_factor = args.signal_factor )
    print "agreement=",agreement

    if args.plottype=="plain":
        plot.getPlot()
    else:
        plot.getSpecialPlot( what=args.plottype, nthpoint=args.nthpoint, signal_factor = args.signal_factor )
    plot.savePlot()

    # import IPython
    # IPython.embed()
