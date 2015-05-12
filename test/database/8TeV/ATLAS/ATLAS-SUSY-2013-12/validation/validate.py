#!/usr/bin/env python

import sys
sys.path.insert(0,"../../../../../smodels-utils/")
sys.path.insert(0,"../../../../../smodels/")


from validation.plotProducer import validateExpRes, getExpIdFromPath
from smodels.experiment.databaseObjects import DataBase
import logging
from smodels.theory.crossSection import logger as cl
from smodels.theory.slhaDecomposer import logger as dl
from smodels.experiment.txnameObject import logger as tl
cl.setLevel(level=logging.DEBUG) 
dl.setLevel(level=logging.DEBUG)
tl.setLevel(level=logging.DEBUG)


print "exp id=",getExpIdFromPath()

database = DataBase("../../../../")
#How to validate all plots for all Txnames in one ExpRes:
expRes = database.getExpResults(analysisIDs=[getExpIdFromPath()],datasetIDs=[None])
slhamain = '../../../../../smodels-utils/slha/'
kfactorDict = { "TChiChipmSlepL" : 1.25, "TChiWZ": 1.25, "TChiWW": 1.25, "TChiChipmStauL" : 1.25, "TChiWH" : 1.25, "TChiWZoff" : 1.25}
validateExpRes(expRes,slhamain, kfactorDict = kfactorDict )

