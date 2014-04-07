#!/usr/bin/env python

"""
.. module:: externalPythia6
        :synopsis: Wrapper code for pythia6

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import setPath
from externalTool import ExternalTool
from tools import logger
import tempfile

class ExternalPythia6(ExternalTool):
    def __init__(self,
                 executable_path="<install>/tools/external/pythia6/pythia_lhe",
                 test_params_path="<install>/etc/pythia6_test.cfg", 
                 src_path="<install>/pythia6/", verbose=False):
        """ 
        :param executable_path: location of executable, full path (pythia_lhe)
        :param test_params_path: location of the test config file, full path
        (external_lhe.test)
        
        """ 
        self.name="pythia6"
        self.executable_path=self.absPath ( executable_path )
        self.test_params_path=self.absPath ( test_params_path )
        self.src_path=self.absPath ( src_path )
        self.verbose=False
        self.tempdir=tempfile.mkdtemp()

    def unlink( self ): ## remove cruft
        import os
        logger.debug ( "unlinking %s " % self.tempdir )
        for File in [ "fort.61", "fort.68" ]:
            if os.path.exists ( self.tempdir + "/" + File ):
                os.unlink ( self.tempdir + "/" + File )
        if os.path.exists ( self.tempdir ):
            os.rmdir ( self.tempdir )

    def run(self, cfg_file):
        """
        Run Pythia.
        
        :params cfg_file: config file used
        :returns: stdout and stderr, or error message
        
        """
        import os, commands
        cfg=os.path.abspath ( cfg_file )
        if not os.path.exists( cfg ): 
            return "config file ``%s'' not found" % ( cfg )
        cmd="cd %s ; %s < %s" % ( self.tempdir, self.executable_path, cfg_file )
        # cmd="%s < %s" % ( self.executable_path, cfg_file )
        Out=commands.getoutput ( cmd )
        out=Out.split("\n")
        self.unlink ()
        return out

    def compile(self):
        """ compile pythia_lhe """
        logger.info("Trying to compile pythia:")
        cmd="cd %s; make" % self.src_path
        import commands
        out=commands.getoutput ( cmd )
        logger.info(out)
        

    def fetch(self, verbose=True):
        """ fetch and unpack tarball """
        import urllib, tarfile
        tempfile="/tmp/pythia.tar.gz"
        f=open( tempfile,"w")
        if verbose:
            logger.info("Fetching tarball ...")
        R=urllib.urlopen("http://smodels.hephy.at/externaltools/pythia/pythia.tar.gz")
        l=R.readlines()
        for line in l:
            f.write ( line )
        R.close()
        f.close()
        if verbose:
            logger.info("done.")
            logger.info("Untarring ...")
        tar=tarfile.open( tempfile )
        for item in tar:
            tar.extract ( item, self.src_path )
        if verbose:
            logger.info("done.")
        

    def checkInstallation(self):
        """ checks if installation of tool looks ok by
                looking for executable and running it """
        import os
        if not os.path.exists( self.executable_path ): 
            logger.error("executable ``%s'' not found" % (self.executable_path))
            return False
        if not os.path.exists( self.test_params_path ): 
            logger.error("config file ``%s'' not found" % (self.test_params_path))
            return False
        if not os.access(self.executable_path, os.X_OK ): 
            logger.error("%s is not executabloe" % self.executable)
            return False
        out=self.run ( self.test_params_path )
        out.pop()
        lines={ -1: " ********* Fraction of events that fail fragmentation cuts =  0.00000 *********" }
        for (nr, line) in lines.items():
            if out[nr].find(line)==-1:
                logger.error("Something is wrong with the setup: " + str(out))
                return False
        self.unlink()
        return True
    

if __name__ == "__main__":
    tool=ExternalPythia6()
    print("installed:" + str(tool.installDirectory()))
    print("check:" + str(tool.checkInstallation()))
    print("run:")
    tool.run ( "../etc/pythia6_test.cfg" )
