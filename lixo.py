#!/usr/bin/env python

import sys
from smodels.theory import slhaDecomposer
from smodels.theory import lheDecomposer
from smodels.tools.physicsUnits import fb, GeV
from smodels.tools.printer import printout
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.databaseObjects import Database

#Set the address of the database folder
database = Database("../smodels-database/")

 listOfExpRes = database.getExpResults(datasetIDs=[None])

for expRes in listOfExpRes:
      print (expRes.getValuesFor('sqrts')/TeV).asNumber()
      sys.exit()
