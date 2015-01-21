#!/usr/bin/env python

"""
.. module:: createTarball
   :synopsis: Script that is meant to create the distribution tarball

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
   
"""

import sys
import commands
import os


RED = "\[31;11m"
RESET = "\[37;0m"


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
    print RED, "Creating temporary directory", dirname, RESET
    o = commands.getoutput("mkdir -p %s" % dirname)
    print o


def rmdir():
    """
    Remove the temporary directory.
    """
    if os.path.exists(dirname):
        print RED, "Removing temporary directory", dirname, RESET
        o = commands.getoutput("rm -rf %s" % dirname)
        print o


def cp():
    """
    Copy files to temporary directory.
    """
    print RED, "Copying the files to", dirname, RESET
    for i in os.listdir("../"):
        if i not in [".git", ".gitignore", "distribution", "test"]:
            # print i
            o = commands.getoutput("cp -r ../%s %s/" % (i, dirname))
    # print o


def rmpyc ():
    """
    Remove .pyc files.
    """
    print RED, "Removing all pyc files ... ", RESET
    o = commands.getoutput("cd %s; rm -f *.pyc */*.pyc */*/*.pyc" % dirname)
    print o


def makeClean ():
    """
    Execute 'make clean' in host directory.
    """
    print RED, "Make clean ....", RESET
    o = commands.getoutput("cd ../lib/ ; make clean")
    print o


def fetchDatabase():
    """
    Execute 'git clone' to retrieve the database.
    """
    print RED, "git clone the database ... ", RESET
    cmd = "cd %s; git clone -b v%s git@smodels.hephy.at:smodels-database ;" \
        " rm -rf smodels-database/.git smodels-database/.gitignore " % \
            (dirname, version)
    o = commands.getoutput(cmd)
    print o


def createTarball():
    """
    Create the tarball.
    """
    print RED, "Create tarball smodels-v%s.tar.gz" % version, RESET
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
    print RED, "Converting the recipes", RESET
    cmd = "cd %s/docs/Manual/recipes/; make convert remove_ipynbs" % dirname
    o = commands.getoutput (cmd)
    print o


def explode ():
    """
    Explode the tarball.
    """
    print RED, "Explode the tarball ...", RESET
    cmd = "tar xzvf smodels-v%s.tar.gz" % version
    o = commands.getoutput (cmd)


def make ():
    """
    Execute 'make' in dirname/lib.
    """
    print RED, "Now run make in dirname/lib ...", RESET
    cmd = "cd %s/lib; make" % dirname
    o = commands.getoutput (cmd)
    print o


def runExample ():
    """
    Execute Example.py.
    """
    print RED, "Now run Example.py ...", RESET
    cmd = "cd %s/; ./Example.py" % dirname
    o = commands.getoutput (cmd)
    print o


def test ():
    """
    Test the tarball, explode it, execute 'make', and 'runSModelS.py'.
    """
    print RED, "--------------------------", RESET
    print RED, "    Test the setup ...    ", RESET
    print RED, "--------------------------", RESET
    rmdir ()
    explode ()
    make ()
    runExample()


def create():
    """
    Create a tarball for distribution.
    """
    print RED, "Creating tarball for distribution, version", version, RESET
    makeClean()
    rmdir()
    mkdir()
    cp()
    rmpyc()
    rmExtraFiles()
    fetchDatabase()
    convertRecipes()
    createTarball()
    test ()
    # rmdir(dirname)


if __name__ == "__main__":
    create()
