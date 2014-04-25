#!/usr/bin/env python

"""
.. module:: externalNllFast
   :synopsis: Wrapper for all nllfast versions.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import setPath
from externalTool import ExternalTool
from tools import logger

class ExternalNllFast(ExternalTool):
    def __init__ ( self, sqrts, nllfast_version, test_params, test_condition ):
        """ 
        :param sqrts: sqrt of s, in TeV, as an integer,
        :param nllfast_version: version of the nllfast tool
        :param test_params: what are the test params we need to run things with?
        :param test_condition: the line that should be the last output line when
        running executable
        :src_path: the path of the source code, for compilation
        """ 
        self.sqrts=int(sqrts)
        self.name="nllfast%d" % sqrts
        self.nllfast_version=nllfast_version
        self.cd_path=self.absPath("<install>/tools/external/nllfast/nllfast-%s/" % \
                        self.nllfast_version )
        self.executable_path=self.cd_path+"/nllfast_%dTeV" % self.sqrts
        self.test_params=test_params
        self.test_condition=test_condition
        self.src_path=self.cd_path

    def compile ( self ):
        """ try to compile tool """
        logger.info("trying to compile %s" % self.name)
        cmd="cd %s; make" % self.src_path
        import commands
        out=commands.getoutput(cmd)
        logger.info(out)
        return True
    
    def fetch(self):
        """ fetch and unpack tarball """
        import urllib, tarfile
        tempfile="/tmp/nllfast7.tar.gz"
        f=open( tempfile,"w")
        url="http://smodels.hephy.at/externaltools/nllfast%d.tar.gz" % self.sqrts
        logger.info ( "fetching tarball from "+url )    
        R=urllib.urlopen(url)
        l=R.readlines()
        for line in l:
            f.write ( line )
        R.close()
        f.close()
        tar=tarfile.open( tempfile )
        for item in tar:
            tar.extract ( item, self.src_path+ "/" )


    def unlink ( self, File ): 
        """ remove File.out """
        return
        """
        import os
        Fname="%s/%s.out" % ( self.cd_path, File )
        if os.path.exists ( Fname ):
            os.unlink ( Fname )
        """

    def run_ ( self, params ):
        """ run nllfast7
            :params params: parameters used (e.g. gg cteq5 .... ) 
             FIXME could have a fancier interface
            :returns: stdout and stderr, or error message
        """
        import commands
        cmd="cd %s; %s %s" % (self.cd_path, self.executable_path, params)
        out=commands.getoutput ( cmd )
        out=out.split("\n")
        return out

    def run ( self,    process, pdf, squarkmass, gluinomass ):
        """ run nllfast
            :params process: which process: st, sb, gg, gdcpl, sdcpl, ss, sg, tot
            :params pdf: cteq=cteq6, mstw2008 
            :params squarkmass: squarkmass, None if squark decoupled
            :params gluinomass: gluinomass, None if gluino decoupled
            :returns: stdout and stderr, or error message
        """
        if not process in [ "st", "sb", "gg", "gdcpl", "sdcpl", "ss", "sg", "tot" ]:
            return None
        if not pdf in [ "cteq", "cteq6", "mstw", "mstw2008" ]: return None
        if not squarkmass: 
            return run_("%s %s %s") % ( process, pdf, gluinomass )
        if not gluinomass: 
            return run_("%s %s %s") % ( process, pdf, squarkmass )
        return run_("%s %s %s %s") % ( process, pdf, squarkmass, gluinomass )

    def checkInstallation ( self ):
        """ checks if installation of tool looks ok by
                looking for executable and running it """
        import os
        if not os.path.exists( self.executable_path ): 
            logger.error ( "executable ``%s'' not found" % ( self.executable_path ) )
            return False
        if not os.access ( self.executable_path, os.X_OK ): 
            logger.error ( "%s is not executable" % self.executable )
            return False
        out=self.run_ ( self.test_params )
        lines={ -1: "500.         600.        0.193E+00    0.450E+00    0.497E+00" }
        if out[-1].find(self.test_condition)==-1:
            logger.error ( "Something is wrong with the setup: "+str(out) )
            return False
        self.unlink ( "gg" )
        return True


class ExternalNllFast7(ExternalNllFast):
    def __init__ ( self ):
        ExternalNllFast.__init__( self, 7, "1.2",
            test_params="gg cteq 500 600", 
            test_condition="500.     600.    0.193E+00  0.450E+00  0.497E+00" )

class ExternalNllFast8(ExternalNllFast):
    def __init__ ( self ):
        ExternalNllFast.__init__( self, 8, "2.1",
            test_params="gg cteq 500 600", 
            test_condition="500.     600.    0.406E+00  0.873E+00  0.953E+00" )

class ExternalNllFast13(ExternalNllFast):
    def __init__ ( self ):
        ExternalNllFast.__init__( self, 13, "3.0",
            test_params="gg cteq 500 600", 
            test_condition="500.     600.    0.394E+01  0.690E+01  0.731E+01" )

class ExternalNllFast14(ExternalNllFast):
    def __init__ ( self ):
        ExternalNllFast.__init__( self, 14, "4.01dcpl",
            test_params="gdcpl cteq 500 600", 
            test_condition="500.    0.235E+02  0.346E+02  0.362E+02" )

class ExternalNllFast33(ExternalNllFast):
    def __init__ ( self ):
        ExternalNllFast.__init__( self, 33, "5.01dcpl", 
            test_params="gdcpl cteq 500 600", 
            test_condition="500.    0.257E+03  0.383E+03  0.393E+03" )

nllFastTools={ 7:ExternalNllFast7(), 8:ExternalNllFast8(),
                            13:ExternalNllFast13(), 14:ExternalNllFast14(), 
                            33:ExternalNllFast33() }

if __name__ == "__main__":
    for (sqrts,tool) in nllFastTools.items():
        print("%s: installed in %s" % (tool.name,tool.installDirectory()))
