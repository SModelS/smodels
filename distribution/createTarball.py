#!/usr/bin/python

"""
.. module:: createTarball
   :synopsis: script that is meant to create the distribution tarball

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
   
"""

import sys
import commands
import os

def produceTarball():
    """ finally, create the tarball itself """
    commands.getoutput( "cd ../../; tar czvf smodels.tar.gz smodels/")

def getVersion():
    """ obtain the smodels version """
    sys.path.append("../")
    from smodels import installation
    return installation.version()

def mkdir( dirname ):
    """ create the directory which we will create the tarball with """
    print "Creating temporary directory",dirname
    o=commands.getoutput ( "mkdir -p %s" % dirname )
    print o

def rmdir( dirname ):
    """ remove the temporary directory """
    if os.path.exists(dirname):
        print "Removing temporary directory",dirname
        o=commands.getoutput ( "rm -r %s" % dirname )
        print o

def cp( dirname ):
    """ copy files to temp directory """
    print "Copying the files to",dirname
    for i in os.listdir("../"):
        if i not in [ ".git", ".gitignore", "distribution", "test" ]:
      #      print i
            o=commands.getoutput ( "cp -r ../%s %s/" % (i,dirname) )
    #print o

def rmpyc ( dirname ):
    """ remove all pyc files """
    print "Removing all pyc files ... "
    o=commands.getoutput("cd %s; rm -f *.pyc */*.pyc */*/*.pyc" % dirname )
    print o

def makeClean ():
    """ perform 'make clean' in host directory """
    print "Make clean ...."
    o=commands.getoutput("cd ../lib/ ; make clean" )
    print o

def fetchDatabase(version,dirname):
    """ git-pull the database """
    print "git pull the database """
    cmd="cd %s; git clone -b v%s git@smodels.hephy.at:smodels-database " % \
            (dirname, version)
    o=commands.getoutput( cmd )
    print o

def createTarball(version,dirname):
    """ finally create the tarball """
    print "Create tarball smodels-v%s.tar.gz" % version
    o=commands.getoutput("tar czvf smodels-v%s.tar.gz %s" % (version, dirname) )
    print o

def rmExtraFiles(dirname):
    """ remove a few more files """

def create():
    """ create a tarball for distribution """
    version=getVersion()
    dirname="smodels-v%s" % version
    print "Creating tarball for distribution, version",version
    makeClean()
    rmdir(dirname)
    mkdir(dirname)
    cp(dirname)
    rmpyc(dirname)
    rmExtraFiles(dirname)
    ## fetchDatabase(version,dirname)
    createTarball(version,dirname)
    rmdir(dirname)


if __name__ == "__main__":
    create()
