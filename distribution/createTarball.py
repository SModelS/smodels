#!/usr/bin/env python

"""
.. module:: createTarball
   :synopsis: script that is meant to create the distribution tarball

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
   
"""

import sys
import commands
import os

RED="[31;11m"
RESET="[37;0m" 

def getVersion():
    """ obtain the smodels version """
    sys.path.append("../")
    from smodels import installation
    return installation.version()

version=getVersion()
dirname="smodels-v%s" % version

def mkdir( ):
    """ create the directory which we will create the tarball with """
    print RED,"Creating temporary directory",dirname,RESET
    o=commands.getoutput ( "mkdir -p %s" % dirname )
    print o

def rmdir( ):
    """ remove the temporary directory """
    if os.path.exists(dirname):
        print RED,"Removing temporary directory",dirname,RESET
        o=commands.getoutput ( "rm -rf %s" % dirname )
        print o

def cp( ):
    """ copy files to temp directory """
    print RED,"Copying the files to",dirname,RESET
    for i in os.listdir("../"):
        if i not in [ ".git", ".gitignore", "distribution", "test" ]:
      #      print i
            o=commands.getoutput ( "cp -r ../%s %s/" % (i,dirname) )
    #print o

def rmpyc ( ):
    """ remove all pyc files """
    print RED,"Removing all pyc files ... ",RESET
    o=commands.getoutput("cd %s; rm -f *.pyc */*.pyc */*/*.pyc" % dirname )
    print o

def makeClean ():
    """ perform 'make clean' in host directory """
    print RED,"Make clean ....",RESET
    o=commands.getoutput("cd ../lib/ ; make clean" )
    print o

def fetchDatabase():
    """ git-pull the database """
    print RED,"git clone the database ... ",RESET
    cmd="cd %s; git clone -b v%s git@smodels.hephy.at:smodels-database ;" \
        " rm -rf smodels-database/.git smodels-database/.gitignore " % \
            (dirname, version)
    o=commands.getoutput( cmd )
    print o

def createTarball():
    """ finally create the tarball """
    print RED,"Create tarball smodels-v%s.tar.gz" % version, RESET
    o=commands.getoutput("tar czvf smodels-v%s.tar.gz %s" % (version, dirname) )
    print o

def rmExtraFiles():
    """ remove a few more files """

def convertRecipes():
    """ convert the recipes from .ipynb to .py and .html """
    print RED,"Converting the recipes", RESET
    cmd="cd %s/docs/Manual/recipes/; make convert remove_ipynbs" % dirname
    o=commands.getoutput ( cmd )
    print o

def explode ( ):
    """ explode the tarball """
    print RED,"Explode the tarball ...", RESET
    cmd="tar xzvf smodels-v%s.tar.gz" % version
    o=commands.getoutput ( cmd )

def make ( ):
    """ run make in dirname/lib """
    print RED,"Now run make in dirname/lib ...", RESET
    cmd="cd %s/lib; make" % dirname
    o=commands.getoutput ( cmd )
    print o

def runExample ( ):
    """ run Example.py """
    print RED,"Now run Example.py ...",RESET
    cmd="cd %s/; ./Example.py" % dirname
    o=commands.getoutput ( cmd )
    print o

def test ( ):
    """ test the tarball, explode it, run make, and runSModelS.py """
    print RED,"--------------------------",RESET
    print RED,"    Test the setup ...    ",RESET
    print RED,"--------------------------",RESET
    rmdir ( )
    explode ( )
    make ( )
    runExample()


def create():
    """ create a tarball for distribution """
    print RED,"Creating tarball for distribution, version",version,RESET
    makeClean()
    rmdir()
    mkdir()
    cp()
    rmpyc()
    rmExtraFiles()
    fetchDatabase()
    convertRecipes()
    createTarball()
    test ( )
##     rmdir(dirname)


if __name__ == "__main__":
    create()
