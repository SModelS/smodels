#!/usr/bin/env python

import argparse, types
argparser = argparse.ArgumentParser(description="validates one pair of txname / axes")
argparser.add_argument ( '-T', '--txname', nargs='?', help='txname [T1]',
                type=types.StringType, default='T1' )
argparser.add_argument ( '-a', '--axes', nargs='?', help='axes description [2*Eq(mother,x)_Eq(lsp,y)]',
                type=types.StringType, default='2*Eq(mother,x)_Eq(lsp,y)' )
argparser.add_argument ( '-k', '--kfactor', nargs='?', help='k factor [1.0]',
                type=types.FloatType, default=1.0 )
args=argparser.parse_args() 


import sys,os
sys.path.insert(0,"/home/lessa/smodels-utils/validation")
sys.path.insert(0,"/home/lessa/smodels-utils")
sys.path.insert(0,"/home/lessa/smodels")

from validation.plotProducer import validateTxName,validatePlot,validateExpRes, getExpIdFromPath
from smodels.experiment.databaseObjects import Database
import logging
from smodels.theory.crossSection import logger as cl
from smodels.theory.slhaDecomposer import logger as dl
from smodels.experiment.txnameObject import logger as tl
cl.setLevel(level=logging.DEBUG)
dl.setLevel(level=logging.DEBUG)
tl.setLevel(level=logging.DEBUG)

database = Database("../../../../")

#How to validate all plots for all Txnames in one ExpRes:
expRes = database.getExpResults(analysisIDs=[getExpIdFromPath()],datasetIDs=[None])

## axes="2*Eq(mother,x)_Eq(lsp,y)"
slhamain = '../../../../../smodels-utils/slha/'
## txname="T6bbWW"

print validatePlot(expRes,args.txname,args.axes,slhamain+"%s.tar" % args.txname, kfactor=args.kfactor )
