#!/usr/bin/env python3

"""
.. module:: calchepWrapper
   :synopsis: Wrapper for CalcHEP.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
from smodels.tools.wrapperBase import WrapperBase
from smodels.tools import wrapperBase
from smodels.base.smodelsLogging import logger, setLogLevel
from smodels.base.physicsUnits import fb, pb, TeV, mb
from smodels.base.crossSection import XSectionList, getXsecFromLHEFile
from smodels import installation
from smodels.decomposition.exceptions import SModelSDecompositionError as SModelSError
import os, sys, io, shutil
from typing import Union

try:
    import commands as executor  # python2
except ImportError:
    import subprocess as executor # python3

class CalcHEPWrapper(WrapperBase):
    """
    An instance of this class represents the installation of CalcHEP.
    """

    def __init__(self,
                 configFile: str="<install>/smodels/etc/calchep.cfg",
                 executablePath: str="<install>/smodels/lib/calchep/calchep.exe",
                 srcPath: str="<install>/smodels/lib/calchep/"):
        """
        :param configFile: Location of the config file, full path; copy this
                           file and provide tools to change its content and to provide a template
        :param executablePath: Location of executable, full path (calchep.exe)
        :param srcPath: Location of source code
        """

        WrapperBase.__init__(self)
        self.name = "calchep"
        self.executablePath = self.absPath(executablePath)
        self.executable = None
        self.srcPath = self.absPath(srcPath)
        self.version = self.getCalcHEPVersion()
        includeFile = f"<install>/smodels/lib/calchep/calchep{self.version}/include"
        self.includeFile = self.absPath ( includeFile )
        self.compiler = "C++"
        self.tempdir = None
        self.cfgfile = self.checkFileExists(configFile)
        self.keepTempDir = False
        self.nevents = 10000
        self.sqrts = 13
        self.secondsPerEvent = 10
        self.calchepCard = None

        self.unlink()

    def getCalcHEPVersion( self ) -> str:
        """obtain the CalcHEP version we wish to use, stored in file 'calchepversion'"""
        versionfile = f"{self.srcPath}/calchepversion"
        if not os.path.exists( versionfile ):
            print( f"[installer.py] error cannot determine calchepversion: did not find {versionfile}")
            sys.exit(-1)
        with open( versionfile, "rt") as f:
            ver = f.read()
            ver = ver.strip()
            f.close()
        return ver

    def checkFileExists(self, inputFile: str) -> str:
        """
        Check if file exists, raise an IOError if it does not.

        :returns: absolute file name if file exists.

        """
        nFile = self.absPath(inputFile)
        if not os.path.exists(nFile):
            raise IOError( f"file {nFile} does not exist" )
        return nFile

    def __str__(self):
        """
        Describe the current status

        """
        ret = f"tool: {self.name}\n"
        ret += f"executable: {self.executablePath}\n"
        ret += f"temp dir: {self.tempdir}\n"
        ret += f"nevents: {self.nevents}\n"
        return ret

    def unlink(self, unlinkdir: bool=True):
        """
        Remove temporary files.

        :param unlinkdir: remove temp directory completely

        """
        if self.tempdir == None:
            return
        if self.keepTempDir:
            logger.warning("Keeping everything in " + self.tempdir)
            return
        logger.debug("Unlinking " + self.tempdir)
        for inputFile in ["fort.61", "fort.68", "log"]:
            if os.path.exists(self.tempdir + "/" + inputFile):
                os.unlink(self.tempdir + "/" + inputFile)
        if unlinkdir:
            for inputFile in ["temp.cfg"]:
                os.unlink(self.tempdir + "/" + inputFile)
            if os.path.exists(self.tempdir):
                os.rmdir(self.tempdir)
                self.tempdir = None

    def checkInstallation ( self, compile : bool = True ):
        # super().checkInstallation(compile)
        exists = os.path.exists ( self.includeFile )
        xmldoc = self.getXmldoc()
        sleep = 0.
        while not os.path.exists ( xmldoc ) or not os.path.exists ( self.executablePath ): # if this disappears, start from scratch
            import time
            sleep += .5
            time.sleep ( sleep )
            if sleep > .5 and not os.path.exists ( xmldoc ) or not os.path.exists ( self.executablePath ):
                if compile:
                    # after a few seconds, delete, if compile is true
                    import shutil
                    p = xmldoc.find ( "share" )
                    rm = xmldoc[:p-1]
                    if False:
                        shutil.rmtree ( rm, ignore_errors = True )
                    self.compile()
                exists = False
                break

        if xmldoc == None:
            exists = False
        if exists:
            return True
        if compile:
            self.compile()
        exists = os.path.exists ( self.includeFile )
        return exists


    def run( self, slhaFile: str, lhefile: Union[str,None]=None,
             unlink: bool=True ) -> XSectionList:
        """
        Run CalcHEP.

        :param slhaFile: SLHA file
        :param unlink: clean up temporary files after run?
        :returns: List of cross sections and branching ratios

        """
        if self.maycompile:
            self.checkInstallation( compile = True )
        # Change calchep configuration file, if defined:
        calchepCard_default = self.cfgfile
        if self.calchepCard:            
            self.cfgfile = self.calchepCard

        self.xsecs = {}
        logger.debug("wrapper.run()")
        slha = self.checkFileExists(slhaFile)
        logger.debug("file check: " + slha)
        cfg = self.absPath(self.cfgfile)
        logger.debug("running with cfgfile " + str(cfg))
        ret = getXsecFromCalHEP(lheFile)

        # Reset calchep card to its default value
        if self.calchepCard:
            self.cfgfile = calchepCard_default

        if unlink:
            shutil.rmtree(self.tempDirectory())
            # print ( "rmtree", self.tempDirectory() )

        return ret

    def chmod(self):
        """
        chmod 755 on calchep executable, if it exists.
        Do nothing, if it doesnt exist.
        """
        if not os.path.exists(self.executablePath):
            logger.error( f"{self.executablePath} doesnt exist" )
            return False
        import stat

        mode = stat.S_IRWXU | stat.S_IRWXG | stat.S_IXOTH | stat.S_IROTH
        os.chmod(self.executablePath, mode)
        return True


if __name__ == "__main__":
    setLogLevel("debug")
    tool = CalcHEPWrapper()
    tool.nevents = 10
    logger.info("installed: " + str(tool.installDirectory()))
    logger.info("check: " + wrapperBase.ok(tool.checkInstallation()))
    logger.info( f"seconds per event: {tool.secondsPerEvent}" )
    slhafile = "inputFiles/slha/simplyGluino.slha"
    slhapath = os.path.join(installation.installDirectory(), slhafile)
    logger.info("slhafile: " + slhapath)
    output = tool.run(slhapath, unlink=True)
    logger.info( f"done: {output}" )
