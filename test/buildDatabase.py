#!/usr/bin/python

import sys, os, time
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
from smodels.tools.smodelsLogging import setLogLevel
from smodels.tools.colors import colors
colors.on = True
setLogLevel ( "debug" )

def build ( force=False ):
    Dir = "tinydb/"
    force_load = None
    if force: force_load = "txt"
    d=Database( Dir, force_load = force_load, discard_zeroes = True )

force=False
if len(sys.argv)>1 and sys.argv[1]=="-f":
    force=True
build(force)
