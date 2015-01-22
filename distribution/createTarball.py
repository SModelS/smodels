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
RESET = "\033[7;0m"

def comment( text ):
    print "%s[%s] %s %s" % ( RED, time.asctime(),text,RESET )

def getVersion():
    """
    Obtain the smodels version """
    sys.path.append("../")
    from smodels import installation
    return installation.version()

version = getVersion()
dirname = "smodels-v%s" % version


def mkdir():
    """
    Create a temporary directory for creating the tarball.
    """
    comment ("Creating temporary directory %s" %  dirname )
    o = commands.getoutput("mkdir -p %s" % dirname)
    print o


def rmdir():
    """
    Remove the temporary directory.
    """
    if os.path.exists(dirname):
        comment ( "Removing temporary directory %s" % dirname )
        o = commands.getoutput("rm -rf %s" % dirname)
        print o


def cp():
    """
    Copy files to temporary directory.
    """
    comment ( "Copying the files to %s" % dirname )
    for i in os.listdir("../"):
        if i not in [".git", ".gitignore", "distribution", "test"]:
            # print i
            o = commands.getoutput("cp -r ../%s %s/" % (i, dirname))
    # print o

def clone():
    """
    Git clone smodels itself into dirname, then remove .git, .gitignore, distribution, and test.
    """
    print RED, "Git-cloning smodels in", dirname, RESET
    o = commands.getoutput("cd %s; git clone git@smodels.hephy.at:smodels" % (dirname) )
    print o
    for i in os.listdir( dirname ):
        if i in [".git", ".gitignore", "distribution", "test"]:
            o = commands.getoutput ( "rm -rf %s/%s" % (dirname,i) )

def rmpyc ():
    """
    Remove .pyc files.
    """
    comment ( "Removing all pyc files ... " )
    o = commands.getoutput("cd %s; rm -f *.pyc */*.pyc */*/*.pyc" % dirname)
    print o


def makeClean ():
    """
    Execute 'make clean' in host directory.
    """
    comment ( "Make clean ...." )
    o = commands.getoutput("cd ../lib/ ; make clean")
    print o


def fetchDatabase():
    """
    Execute 'git clone' to retrieve the database.
    """
    comment ( "git clone the database ... " )
    cmd = "cd %s; git clone -b v%s git@smodels.hephy.at:smodels-database ;" \
        " rm -rf smodels-database/.git smodels-database/.gitignore " % \
            (dirname, version)
    o = commands.getoutput(cmd)
    print o


def createTarball():
    """
    Create the tarball.
    """
    comment ( "Create tarball smodels-v%s.tar.gz" % version )
    o = commands.getoutput("tar czvf smodels-v%s.tar.gz %s" % (version, dirname))
    print o


def rmExtraFiles():
    """
    Remove additional files.
    """
    pass


def convertRecipes():
    """
    Compile recipes from .ipynb to .py and .html.
    """
    comment ( "Converting the recipes" )
    cmd = "cd %s/docs/manual/source/recipes/; make convert remove_ipynbs" % dirname
    o = commands.getoutput (cmd)
    print o

def makeDocumentation():
    """
    create the documentation via sphinx """
    comment ( "Creating the documentation" )
    cmd = "cd %s/docs/manual/; make html; rm -r source/" % dirname
    o = commands.getoutput (cmd)
    print o
    cmd = "cd %s/docs/documentation/; make html; rm -r source/" % dirname
    o = commands.getoutput (cmd)
    print o

def explode ():
    """
    Explode the tarball.
    """
    comment ( "Explode the tarball ..." )
    cmd = "tar xzvf smodels-v%s.tar.gz" % version
    o = commands.getoutput (cmd)


def make ():
    """
    Execute 'make' in dirname/lib.
    """
    comment ( "Now run make in dirname/lib ..." )
    cmd = "cd %s/lib; make" % dirname
    o = commands.getoutput (cmd)
    print o


def runExample ():
    """
    Execute Example.py.
    """
    comment ( "Now run Example.py ..." )
    cmd = "cd %s/; ./Example.py" % dirname
    o = commands.getoutput (cmd)
    print o


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
    cmd="ls %s/docs/manual.html" % dirname 
    o = commands.getoutput (cmd)
    print o


def create():
    """
    Create a tarball for distribution.
    """
    comment ( "Creating tarball for distribution, version %s" % version )
    makeClean()
    rmdir()
    mkdir()
    cp()
    # clone()
    rmpyc()
    rmExtraFiles()
    fetchDatabase()
    makeDocumentation()
    convertRecipes()
    createTarball()
    test ()
    # rmdir(dirname)
    testDocumentation()


if __name__ == "__main__":
    create()
