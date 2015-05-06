#!/usr/bin/env python

from __future__ import print_function

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
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.databaseObjects import Database
from unum import ShouldBeUnitlessError

from smodels.tools.statistics import upperLimit


Nobs = 0
Nexp = 3.0
sigmaexp = 2.8
lumi = 20.1/fb

print(upperLimit(Nobs, Nexp, sigmaexp, lumi))
