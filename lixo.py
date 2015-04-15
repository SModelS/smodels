#!/usr/bin/env python

import sys
from smodels.theory import slhaDecomposer
from smodels.theory import lheDecomposer, branch
from smodels.tools.physicsUnits import fb, GeV
from smodels.tools.printer import printout
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.databaseObjects import Database
from smodels.tools.physicsUnits import GeV, fb, TeV


from smodels.tools import modpyslha as pyslha


b = branch.Branch()
b.PIDs = [[10002]]
b2 = b.copy()
b2.PIDs[0].append(101)
print b.PIDs
print b2.PIDs
