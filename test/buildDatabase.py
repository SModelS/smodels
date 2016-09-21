#!/usr/bin/python

import sys
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database

d=Database("./database/")
print d
