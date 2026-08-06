#!/usr/bin/env python3

import glob, shutil, os

def rmOldPickles():
    files = glob.glob ( "**/*.pcl", recursive=True )
    files += glob.glob ( "**/.*.pcl", recursive=True )
    for f in files:
        if not "notebookTests" in f:
            print ( f"[buildPickles] rm {f}" )
            os.unlink ( f )

def buildDatabases():
    from smodels.experiment.databaseObj import Database
    versionfiles = glob.glob ( "*/version" )
    for versionfile in versionfiles:
        dbpath = versionfile.replace("/version","")
        print ( dbpath )
        try:
            db = Database ( dbpath )
        except Exception as e:
            print ( f"error: {e}" )

if __name__ == "__main__":
    rmOldPickles()
    buildDatabases()
