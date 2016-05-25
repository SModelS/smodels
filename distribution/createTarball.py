#!/usr/bin/env python

"""
.. module:: createTarball
   :synopsis: Script that is meant to create the distribution tarball

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
   
"""

import sys
import commands
import os
import time

RED = "\033[31;11m"
GREEN = "\033[32;11m"
RESET = "\033[7;0m"

dummyRun=True

def comment( text ):
    print "%s[%s] %s %s" % ( RED, time.asctime(),text,RESET )
    f=open("create.log","a")
    f.write (  "[%s] %s\n" % ( time.asctime(),text ) )
    f.close()

def isDummy ( ):
    if dummyRun:
        comment ( "DUMMY RUN!!!!" )

def run ( cmd ):
    print "%scmd: %s%s" % (GREEN,cmd,RESET)
    f=open("create.log","a")
    f.write ( "cmd: %s\n" % (cmd) )
    o=commands.getoutput ( cmd )
    print o
    f.write ( o + "\n" )
    f.close()

def getVersion():
    """
    Obtain the smodels version """
    sys.path.insert(0,"../")
    from smodels import installation
    return installation.version()

version = getVersion()
dirname = "smodels-v%s" % version
fastlimdir = "smodels-fastlim-v%s" % version

def rmlog():
    """ clear the log file """
    cmd="rm -f create.log"
    commands.getoutput ( cmd )

def mkdir():
    """
    Create a temporary directory for creating the tarball.
    """
    for i in ( dirname, fastlimdir ):
        comment ("Creating temporary directory %s" % i )
        run ( "mkdir -p %s" % dirname )

def rmdir():
    """
    Remove the temporary directories
    """
    for i in ( dirname, fastlimdir ):
        if os.path.exists(i):
            comment ( "Removing temporary directory %s" % i )
            run ("rm -rf %s" % i )

def cp():
    """
    Copy files to temporary directory.
    """
    comment ( "Copying the files to %s" % dirname )
    for i in os.listdir("../"):
        if i not in [".git", ".gitignore", "distribution", "test", "__pycache__" ]:
            run ("cp -r ../%s %s/" % (i, dirname))

def clone():
    """
    Git clone smodels itself into dirname, then remove .git, .gitignore,
    distribution, and test.
    """
    comment ( "Git-cloning smodels into %s (this might take a while)" % dirname )
    cmd = "git clone git@smodels.hephy.at:smodels %s" % (dirname)
    if isDummy:
        cmd = "cp -a ../../smodels-v%s/* %s" % ( version, dirname )
    run ( cmd )
    for i in os.listdir( dirname ):
        if i in [".git", ".gitignore", "distribution", "test"]:
            run ( "rm -rf %s/%s" % (dirname,i) )

def rmpyc ():
    """
    Remove .pyc files.
    """
    comment ( "Removing all pyc files ... " )
    run ("cd %s; rm -f *.pyc */*.pyc */*/*.pyc" % dirname)


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
    comment ( "git clone the database (this might take a while)" )
    cmd = "cd %s; git clone -b v%s git@smodels.hephy.at:smodels-database"  % \
            (dirname, version)
    if isDummy:
        cmd = "cd %s; cp -a ../../../smodels-database-v%s smodels-database" % \
              ( dirname, version )
    run ( cmd )
    rmcmd = "cd %s/smodels-database; " \
            "rm -rf .git .gitignore *.py *.sh *.tar *.pyc" % \
             ( dirname )
    run ( rmcmd )

def splitDatabase():
    """
    Split up between the official database and the optional database
    """
    comment ( "Now move all the non-official entries in the database." )
    cwd=os.getcwd()
    comment ( "debug cwd: %s" % cwd )
    comment ( "debug dirname: %s" % dirname )
    dflag=""
    #if isDummy():
    #    dflag="-d"

    cmd = "cd %s/smodels-database/; %s/moveFastlimResults.py %s" % \
          ( dirname, cwd, dflag )
    run ( cmd )

    cmd = "cp ../smodels-fastlim/smodels-fastlim.tar.gz %s/smodels-fastlim-v%s.tar.gz" % ( cwd, version )
    run ( cmd )

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
    extras = [ "inputFiles/slha/nobdecay.slha", "inputFiles/slha/lightEWinos.slha", "docs/documentation/smodels.log" ]
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
    cmd = "cd %s/docs/manual/; make html; rm -r make.bat Makefile source " % dirname
    run (cmd)
    cmd = "cd %s/docs/documentation/; make html; rm -r make.bat  Makefile source update" % dirname
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
    cmd = "cd %s/; ./Example.py" % dirname
    run (cmd)

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
    ## cp()
    clone() ## ... clone smodels into it ...
    fetchDatabase() 
    splitDatabase() ## split database into official and optional
    convertRecipes()
    makeDocumentation()
    rmExtraFiles() ## ... remove unneeded files ...
    rmpyc() ## ...  remove the pyc files created by makeDocumentation ...
    createTarball() ## here we go! create!
    test ()
    # rmdir(dirname)
    testDocumentation()
    isDummy()


if __name__ == "__main__":
    create()
