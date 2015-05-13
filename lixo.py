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
import random

r = 10.
for i in range(10000):  
      x = random.random()*fb
      ul = random.random()*fb
      r1 = x/ul
      if r is None or r1 < r:
            r = r1

print r
