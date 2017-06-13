#!/usr/bin/python

import sys, os, time
sys.path.insert(0,"../")
from smodels.experiment.databaseObj import Database
from smodels.tools.smodelsLogging import setLogLevel
from smodels.tools.colors import colors
colors.on = True
# setLogLevel ( "debug" )

def build ( force=False, dir=None ):
    force_load = None
    if force: force_load = "txt"
    d=Database( dir, force_load = force_load, discard_zeroes = True, 
                progressbar=True )
    if force:
        d.createBinaryFile()
    print ( d )

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Simple database builder (i.e. pickelizer)")
    parser.add_argument('-f', '--force', help='force text mode', action='store_true')
    parser.add_argument('-d', '--dir', help='database dir', default="./tinydb", type=str )
    args=parser.parse_args()
    build( args.force, args.dir )
