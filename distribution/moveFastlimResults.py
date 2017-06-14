#!/usr/bin/python

""" A script that splits up the results into official and fastlim.
    (maybe later even more packages ) """

import os
import commands
import sys
sys.path.insert(0,"." )

# from createTarballs import clearGlobalInfos

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
    # clearGlobalInfos ( fastlimdir )

def createFastlimTarball():
    cmd = "cd %s; tar czvf ../smodels-fastlim.tgz ./" % fastlimdir
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
                gi = fullpath + "/globalInfo.txt"
                if not os.path.exists ( gi ):
                    continue
                gif=open ( gi )
                lines=gif.readlines()
                for line in lines:
                    if "fastlim" in line and "contact" in line:
                        isFastlim ( fullpath, dryrun )
                        break
                gif.close()

def error ( text ):
    print "ERROR: %s" % text

def moveBibFile ( dryrun ):
    """ move fastlim-specific bibliography file """
    fastlim_bib = "references-fastlim.bib"
    if not os.path.exists ( fastlim_bib ):
        error ( "%s is missing!" % fastlim_bib )
    else:
        cmd = "mv %s %s" % ( fastlim_bib, fastlimdir )
        run ( cmd, dryrun )
        
def moveReadmeFile ( dryrun ):
    """ move fastlim-specific README file """
    fastlim_readme = "README_fastlim"
    if not os.path.exists ( fastlim_readme ):
        error ( "%s is missing!" % fastlim_readme )
    else:
        cmd = "mv %s %s" % ( fastlim_readme, os.path.join(fastlimdir,'README' ))
        run ( cmd, dryrun )        

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
    moveBibFile ( args.dryrun )
    moveReadmeFile ( args.dryrun )
    createFastlimTarball()

