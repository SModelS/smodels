#!/usr/bin/env python

import sys
from smodels.theory import slhaDecomposer
from smodels.theory import lheDecomposer
from smodels.tools.physicsUnits import fb, GeV
from smodels.tools.printer import printout
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.databaseObjects import Database
from smodels.tools.physicsUnits import GeV, fb, TeV


from smodels.tools import modpyslha as pyslha

slhafile = 'inputFiles/slha/lightEWinos.slha'
res = pyslha.readSLHAFile(slhafile)
print res
