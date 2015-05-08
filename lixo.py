#!/usr/bin/env python

"""
.. module:: Example
   :synopsis: Basic main file example for using SModelS.
   
   This file must be run under the installation folder.

"""

""" Import basic functions (this file must be executed in the installation folder) """

import sys
from smodels.theory import slhaDecomposer
from smodels.theory import lheDecomposer
from smodels.tools.physicsUnits import fb, GeV
from smodels.tools.printer import printout
from smodels.theory.theoryPrediction import TheoryPredictionList,TheoryPrediction
from smodels.experiment.databaseObjects import Database
from unum import ShouldBeUnitlessError

a1 = TheoryPrediction()
a2 = TheoryPrediction()
a3 = TheoryPrediction()
a = TheoryPredictionList()
a._theoryPredictions = [a1,a2,a3]

b1 = TheoryPrediction()
b2 = TheoryPrediction()
b = TheoryPredictionList()
b._theoryPredictions = [b1,b2]

print len(a),len(b)

c = a + b

print len(sum([a,b]))
