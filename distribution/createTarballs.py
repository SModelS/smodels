#!/usr/bin/env python

"""
.. module:: createTarballs
   :synopsis: Script that is meant to create the distribution tarballs

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
import commands
import os
import time

def getVersion():
    """
    Obtain the smodels version """
    sys.path.insert(0,"../")
    from smodels import installation
    # return "1.1.0"
    return installation.version()

dummyRun=False ## True
version = getVersion()
dirname = "smodels-v%s" % version
# fastlimdir = "smodels-v%s-fastlim-1.0" % version[:3]

RED = "\033[31;11m"
GREEN = "\033[32;11m"
YELLOW = "\033[33;11m"
RESET = "\033[7;0m"


def comment( text, urgency="info" ):
    col=YELLOW
    pre=""
    if "err" in urgency.lower():
        pre="ERROR: "
        col=RED
    print "%s[%s] %s%s %s" % ( col, time.asctime(), pre, text, RESET )
    f=open("create.log","a")
    f.write (  "[%s] %s\n" % ( time.asctime(),text ) )
    f.close()
    if col == RED:
        sys.exit(-1)

def isDummy ( ):
    if dummyRun:
        comment ( "DUMMY RUN!!!!" )
    return dummyRun

def run ( cmd ):
    cmd=cmd.strip()
    print ( "%scmd: %s%s" % (GREEN,cmd,RESET) )
    f=open("create.log","a")
    f.write ( "cmd: %s\n" % (cmd) )
    o=commands.getoutput ( cmd )
    if len(o)>0:
        print (o)
        f.write ( o + "\n" )
    f.close()
    return o

def removeNonValidated():
    """ remove all non-validated analyses from
        database """
    comment ( "Now remove non-validated results." )
    from smodels.experiment.databaseObj import Database
    # d = Database ( "%s/smodels-database" % dirname )
    d = Database ( "%s/smodels-database" % dirname, force_load = "txt" )
    ers = d.expResultList
    comment ( "Loaded the database with %d results." % ( len(ers) ) )
    for er in ers:
        if hasattr ( er.globalInfo, "private" ) and er.globalInfo.private:
            comment ( "%s is private. delete!" % ( er.globalInfo.id ) )
            cmd = "rm -r %s" % ( er.path )
            run ( cmd )
        else:
            hasDataSets=False
            for dataset in er.datasets:
                hasTxNames=False
                for txn in dataset.txnameList:
#                    if txn.validated in [ None, False ]:
                    if txn.validated in [ False ]:
                        comment ( "%s/%s/%s is not validated. Delete it." % \
                                  ( er, dataset, txn ) )
                        cmd="rm %s" % txn.path
                        run ( cmd )
                    else:
                        hasTxNames=True
                if not hasTxNames:
                        comment ( "%s/%s has no validated txnames. remove folder." %\
                                  (er, dataset ) )
                        cmd = "rm -rf %s" % dataset.path
                        run ( cmd )
                if hasTxNames:
                    hasDataSets=True
            if not hasDataSets:
                comment ( "%s has no validated datasets. remove folder." % \
                          (er) )
                cmd = "rm -rf %s" % er.path
                run ( cmd )
    # comment ( "base=%s" % d.base )
    for tev in os.listdir ( d.base ):
        fullpath = os.path.join ( d.base, tev )
        if not os.path.isdir ( fullpath ):
            continue
        tevHasResults=False
        for experiment in os.listdir ( fullpath ):
            exppath = os.path.join ( fullpath, experiment )
            if not os.path.isdir ( exppath ):
                continue
            if os.listdir ( exppath ) == []:
                comment ( "%s/%s is empty. Delete it!" % ( tev, experiment ) )
                cmd = "rm -rf %s" % exppath
                run ( cmd )
            else:
                tevHasResults=True
        if not tevHasResults:
            comment ( "%s is empty. Delete it!" % ( tev ) )
            cmd = "rm -rf %s" % fullpath
            run ( cmd )

def rmlog():
    """ clear the log file """
    cmd="rm -f create.log"
    commands.getoutput ( cmd )
    cmd="rm -f %s/smodels-database/create.log" % dirname
    commands.getoutput ( cmd )


def mkdir():
    """
    Create a temporary directory for creating the tarball.
    """
    ## for i in ( dirname, fastlimdir ):
    for i in ( [ dirname ] ): ## , fastlimdir ):
        comment ("Creating temporary directory %s" % i )
        run ( "mkdir -p %s" % i )

def rmdir():
    """
    Remove the temporary directories
    """
    for i in ( dirname, ): ## fastlimdir ):
        if os.path.exists(i):
            comment ( "Removing temporary directory %s" % i )
            run ("rm -rf %s" % i )

def clone():
    """
    Git clone smodels itself into dirname, then remove .git, .gitignore,
    distribution, and test.
    """
    comment ( "Git-cloning smodels into %s (this might take a while)" % dirname )
    cmd = "git clone -b v%s git@smodels.hephy.at:smodels %s" % (version, dirname)
#     cmd = "git clone git@smodels.hephy.at:smodels %s" % (dirname)
    if dummyRun:
        cmd = "cp -a ../../smodels-v%s/* %s" % ( version, dirname )
    run ( cmd )
    for i in os.listdir( dirname ):
        if i in [".git", ".gitignore", "distribution", "test" ]:
            run ( "rm -rf %s/%s" % (dirname,i) )

def rmpyc ():
    """
    Remove .pyc files.
    """
    comment ( "Removing all pyc files ... " )
    run ("cd %s; rm -f *.pyc */*.pyc */*/*.pyc" % dirname )


def makeClean ():
    """
    Execute 'make clean' in host directory.
    """
    comment ( "Make clean ...." )
    run ("cd ../lib/ ; make clean")

def fetchDatabase():
    """
    Execute 'git clone' to retrieve the database.
    """
    dbversion = version
    if dbversion.find("-")>0:
        dbversion=dbversion[:dbversion.find("-")]
    comment ( "git clone the database (this might take a while)" )
    cmd = "cd %s; git clone -b v%s git@smodels.hephy.at:smodels-database"  % \
            (dirname, dbversion)
            
    if dummyRun:
        cmd = "cd %s; cp -a ../../../smodels-database-v%s smodels-database" % \
              ( dirname, version )
    run ( cmd )
    rmcmd = "cd %s/smodels-database; " \
            "rm -rf .git .gitignore *.py *.sh *.tar *.pyc" % \
             ( dirname )
    run ( rmcmd )


def clearGlobalInfo ( filename ):
    f=open(filename)
    lines=f.readlines()
    f.close()
    g=open("/tmp/tmp.txt","w")
    skip = [ "publishedData", "comment", "private", "checked", \
             "prettyName", "susyProcess", "dataUrl" ]
    #skip.append ( "validated" )
    #skip.append ( "axes" )
    for line in lines:
        to_skip = False
        for s in skip:
            if line.find(s)==0:
                to_skip = True
        if not to_skip:
            g.write ( line )
    g.close()
    cmd = "cp /tmp/tmp.txt %s" % filename
    run ( cmd )

def cleanDatabase():
    """
    Clean up the database, e.g. remove orig and validation folders
    """
    comment ( "Now cleaning up database in %s/smodels-database" % dirname )
    fullpath = "%s/smodels-database" % dirname
    walker = os.walk ( fullpath )
    for record in walker:
        File=record[0]
        # comment ( "Now in %s" % File )
        # comment ( "Now in %s: %s" % (File, record[1] ) )
        removals = [ "orig", ".git", "validation" ]
        rmFiles = [ "run_convert.sh", "checkFastlimValidation.py",  \
                    "checkFastlimValidation.ipynb", "convert.py","convertCMS.py", "sms.root", "general.comment" ]
        for r in removals:
            if r in File:
                cmd = "rm -rf %s" % File
                run ( cmd )
        for rf in rmFiles:
            fullpath = os.path.join ( File, rf )
            if os.path.exists ( fullpath):
                os.unlink ( fullpath )

def clearGlobalInfos( path ):
    walker = os.walk ( path )
    for record in walker:
        File=record[0]
        # print ( "record=",record )
        for i in record[2]:
            if i[0]=="T" and i[-4:]==".txt":        
                fullpath = os.path.join ( File, i )
                clearGlobalInfo ( fullpath )
        gIpath = os.path.join ( File, "globalInfo.txt" )
        if os.path.exists ( gIpath ):
            clearGlobalInfo ( gIpath )

def splitDatabase():
    """
    Split up between the official database and the optional database
    """
    comment ( "Now move all the non-official entries in the database." )
    cwd=os.getcwd()
    comment ( "debug cwd: %s" % cwd )
    comment ( "debug dirname: %s" % dirname )
    cmd = "cd %s/smodels-database/; %s/moveFastlimResults.py" % \
          ( dirname, cwd )
    run ( cmd )
    cmd = "mv ./smodels-fastlim.tgz %s/smodels-v%s-fastlim-1.0.tgz" % \
          ( cwd, version[:3] )
    run ( cmd )

def clearTxtFiles():
    clearGlobalInfos ( "%s/smodels-database/" % dirname )
    clearGlobalInfos ( "./smodels-fastlim/" )

def createTarball():
    """
    Create the tarball.
    """
    comment ( "Create tarball smodels-v%s.tgz" % version )
    run ("tar czvf smodels-v%s.tgz %s" % (version, dirname))

def rmExtraFiles():
    """
    Remove additional files.
    """
    comment ( "Remove a few unneeded files" )
    extras = [ "inputFiles/slha/nobdecay.slha", "inputFiles/slha/broken.slha", 
               "docs/documentation/smodels.log", "inputFiles/slha/complicated.slha" ]
    for i in extras:
        cmd = "rm -rf %s/%s" % ( dirname, i )
        run ( cmd )

def convertRecipes():
    """
    Compile recipes from .ipynb to .py and .html.
    """
    comment ( "Converting the recipes" )
    cmd = "cd %s/docs/manual/source/recipes/; make convert remove_ipynb" % dirname
    run (cmd)

def makeDocumentation():
    """
    create the documentation via sphinx """
    comment ( "Creating the documentation" )
    cmd = "cd %s/docs/manual/; make clean html; rm -r make.bat Makefile source " % dirname
    run (cmd)
    cmd = "cd %s/docs/documentation/; make clean html; rm -r make.bat  Makefile source update" % dirname
    run (cmd)

def explode ():
    """
    Explode the tarball.
    """
    comment ( "Explode the tarball ..." )
    cmd = "tar xzvf smodels-v%s.tgz" % version
    run (cmd)

def make ():
    """
    Execute 'make' in dirname/lib.
    """
    comment ( "Now run make in dirname/lib ..." )
    cmd = "cd %s/lib; make" % dirname
    run (cmd)

def runExample ():
    """
    Execute Example.py.
    """
    comment ( "Now run Example.py ..." )
    cmd = "cd %s/; ./Example.py 2>&1 | tee out.log" % dirname
    run (cmd)
    comment ( "Now check diff" )
    cmd = "diff -u default.log %s/out.log" % dirname
    d = run ( cmd )
    if len ( d ) > 0:
        comment ( "Example test failed!!", "error" )
    else:
        comment ( "Looking good." )


def test ():
    """
    Test the tarball, explode it, execute 'make', and 'runSModelS.py'.
    """
    comment ( "--------------------------" )
    comment ( "    Test the setup ...    " )
    comment ( "--------------------------" )
    rmdir ()
    explode ()
    make ()
    runExample()

def testDocumentation():
    """ Test the documentation """
    comment ( "Test the documentation" )
    cmd="ls %s/docs/manual/build/html/index.html" % dirname
    run (cmd)

def create():
    """
    Create a tarball for distribution.
    """
    isDummy()
    rmlog() ## first remove the log file
    comment ( "Creating tarball for distribution, version %s" % version )
    makeClean()
    rmdir()
    mkdir() ## .. then create the temp dir
    clone() ## ... clone smodels into it ...
    fetchDatabase() ## git clone the database
    cleanDatabase() ## clean up database, remove orig, validated
    splitDatabase() ## split database into official and optional
    removeNonValidated() ## remove all non-validated analyses
    clearTxtFiles() ## now clear up all txt files
    convertRecipes()
    makeDocumentation()
    rmExtraFiles() ## ... remove unneeded files ...
    rmpyc() ## ...  remove the pyc files created by makeDocumentation ...
    rmlog()  ##Make sure to remove log files
    createTarball() ## here we go! create!
    test ()    
    # rmdir(dirname)
    testDocumentation()
    isDummy()

if __name__ == "__main__":
    create()
    # clearGlobalInfos ( "%s/smodels-database/" % dirname )
    # fetchDatabase()
    # cleanDatabase()
    # splitDatabase()
    # cleanDatabase()
    # clearGlobalInfo ( "./globalInfo.txt" )
    # runExample()
