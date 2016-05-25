#!/usr/bin/python

""" A script that splits up the results into official and fastlim.
    (maybe later even more packages ) """

import os
import commands

fastlimdir = "../../smodels-fastlim"

def run ( cmd, dryrun=False ):
    if dryrun:
        print "Dry-run: skipping %s." % cmd
    else:
        print "Executing: %s." % cmd
        commands.getoutput ( cmd )

def backupScript():
    ## first we copy ourself to /tmp
    cmd="cp ./moveFastlimResults.py /tmp/"
    commands.getoutput ( cmd )

def rmDirs():
    if os.path.exists ( fastlimdir ):
        cmd = "rm -r %s" % fastlimdir
        run ( cmd )

def mkDirs():
    cmd="mkdir %s" % fastlimdir
    run ( cmd )

def isFastlim ( path, dryrun ):
    dname = os.path.dirname ( path )
    bname = os.path.basename ( path )
    print "%s is fastlim!" % bname
    cmd = "mkdir -p %s/%s" % ( fastlimdir, dname )
    run ( cmd )
    cmd = "mv %s %s/%s" % ( path, fastlimdir, dname )
    if dryrun:
        cmd = "cp -a %s %s/%s" % ( path, fastlimdir, dname )
    run ( cmd )
    cmd = "rm -r %s/%s/*/orig" % ( fastlimdir, path )
    run ( cmd )
    cmd = "rm -r %s/%s/*/convert.py" % ( fastlimdir, path )
    run ( cmd )
    cmd = "rm -r %s/%s/validation" % ( fastlimdir, path )
    run ( cmd )
    cmd = "rm -r %s/%s/sms.root" % ( fastlimdir, path )
    run ( cmd )

def createFastlimTarball():
    cmd = "cd %s; tar czvf smodels-fastlim.tgz ./" % fastlimdir
    run ( cmd )

## now traverse the *TeV dirs
def traverse( dryrun ):
    for i in os.listdir("."):
        if not os.path.isdir ( i ) or i in [ ".git" ]:
            continue
        for j in os.listdir ( i ):
            fulldir = os.path.join ( i, j )
            if not os.path.isdir ( fulldir ):
                continue
            for analysis in os.listdir ( fulldir ):
                fullpath = os.path.join ( fulldir, analysis )
                gif=open ( fullpath + "/globalInfo.txt" )
                lines=gif.readlines()
                for line in lines:
                    if "fastlim" in line:
                        isFastlim ( fullpath, dryrun )
                        break
                gif.close()

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument( '-d', '--dryrun', 
            help='Dry-run, dont actuall move or create anything',
            action='store_true')
    args = ap.parse_args()
    backupScript()
    rmDirs()
    mkDirs()
    traverse( args.dryrun )
    createFastlimTarball()
