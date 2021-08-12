#!/usr/bin/env python3

"""
.. module:: pythia8Wrapper
   :synopsis: Wrapper for pythia8.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
from smodels.tools.wrapperBase import WrapperBase
from smodels.tools import wrapperBase
from smodels.tools.smodelsLogging import logger, setLogLevel
from smodels.tools.physicsUnits import fb, pb, TeV, mb
from smodels.theory.crossSection import getXsecFromLHEFile
from smodels import installation
from smodels.tools.pythia8particles import particles
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
import os, sys, io, shutil

try:
    import commands as executor #python2
except ImportError:
    import subprocess as executor # python3

class Pythia8Wrapper(WrapperBase):
    """
    An instance of this class represents the installation of pythia8.
    """

    def __init__(self,
                 configFile="<install>/smodels/etc/pythia8.cfg",
                 executablePath="<install>/smodels/lib/pythia8/pythia8.exe",
                 srcPath="<install>/smodels/lib/pythia8/"):
        """
        :param configFile: Location of the config file, full path; copy this
                           file and provide tools to change its content and to provide a template
        :param executablePath: Location of executable, full path (pythia8.exe)
        :param srcPath: Location of source code
        """

        WrapperBase.__init__(self)
        self.name = "pythia8"
        self.executablePath = self.absPath(executablePath)
        self.executable = None
        self.srcPath = self.absPath(srcPath)
        self.compiler = "C++"
        self.tempdir = None
        self.cfgfile = self.checkFileExists(configFile)
        self.keepTempDir = False
        self.nevents = 10000
        self.sqrts = 13
        self.secondsPerEvent = 10
        self.pythiacard = None

        self.unlink()

    def checkFileExists(self, inputFile):
        """
        Check if file exists, raise an IOError if it does not.

        :returns: absolute file name if file exists.

        """
        nFile = self.absPath(inputFile)
        if not os.path.exists(nFile):
            raise IOError("file %s does not exist" % nFile)
        return nFile



    def __str__(self):
        """
        Describe the current status

        """
        ret = "tool: %s\n" % (self.name)
        ret += "executable: %s\n" % (self.executablePath)
        ret += "temp dir: %s\n" % self.tempdir
        ret += "nevents: %d\n" % self.nevents
        return ret


    def unlink(self, unlinkdir=True):
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

    def run( self, slhaFile, lhefile=None, unlink=True ):
        """
        Run pythia8.

        :param slhaFile: SLHA file
        :param lhefile: option to write LHE output to file; if None, do not write
                        output to disk. If lhe file exists, use its events for
                        xsecs calculation.
        :param unlink: clean up temporary files after run?
        :returns: List of cross sections

        """
        #Change pythia configuration file, if defined:
        if self.pythiacard:
            pythiacard_default = self.cfgfile
            self.cfgfile = self.pythiacard


        self.xsecs = {}
        logger.debug ( "wrapper.run()" )
        slha = self.checkFileExists(slhaFile)
        logger.debug ( "file check: " + slha )
        cfg = self.absPath(self.cfgfile)
        logger.debug("running with cfgfile " + str(cfg))
        lhefile = self.tempDirectory() + "/events.lhe"
        cmd = "%s -n %d -f %s -s %d -c %s -l %s" % \
             ( self.executablePath, self.nevents, slha, self.sqrts, cfg, lhefile )
        xmldoc = self.executablePath.replace ( "pythia8.exe", "xml.doc" )
        logger.debug ( "exe path=%s" % self.executablePath )
        self.checkInstallation ( compile=True )
        if os.path.exists (xmldoc ):
            logger.debug ( "xml.doc found at %s." % xmldoc )
            with open ( xmldoc ) as f:
                xmlDir = f.read()
                toadd = os.path.join ( os.path.dirname ( xmldoc ) , xmlDir.strip() )
                logger.debug ( "adding -x %s" % toadd )
                cmd += " -x %s" % toadd
        logger.debug("Now running ''%s''" % str(cmd) )
        out = executor.getoutput(cmd)
        logger.debug ( "out=%s" % out )
        if not os.path.isfile(lhefile):
            raise SModelSError( "LHE file %s not found" % lhefile )
        lheF = open(lhefile,'r')
        lhedata = lheF.read()
        lheF.close()
        os.remove(lhefile)
        if not "<LesHouchesEvents" in lhedata:
            raise SModelSError("No LHE events found in pythia output")
        if not unlink:
            tempfile = self.tempDirectory() + "/log"
            f = open( tempfile, "w")
            f.write (cmd + "\n\n\n")
            f.write (out + "\n")
            f.write (lhedata + "\n")
            f.close()
            logger.info ( "stored everything in %s" % tempfile )

        # Create memory only file object
        if sys.version[0]=="2":
            lhedata = unicode( lhedata )
        lheFile = io.StringIO(lhedata)
        ret = getXsecFromLHEFile(lheFile)

        #Reset pythia card to its default value
        if self.pythiacard:
            self.cfgfile = pythiacard_default

        if unlink:
            shutil.rmtree ( self.tempDirectory() )
            # print ( "rmtree", self.tempDirectory() )

        return ret


    def chmod(self):
        """
        chmod 755 on pythia executable, if it exists.
        Do nothing, if it doesnt exist.
        """
        if not os.path.exists ( self.executablePath ):
            logger.error("%s doesnt exist" % self.executablePath )
            return False
        import stat
        mode = stat.S_IRWXU | stat.S_IRWXG | stat.S_IXOTH | stat.S_IROTH
        os.chmod ( self.executablePath, mode )
        return True

if __name__ == "__main__":
    setLogLevel ( "debug" )
    tool = Pythia8Wrapper()
    tool.nevents=10
    logger.info("installed: " + str(tool.installDirectory()))
    logger.info("check: " + wrapperBase.ok(tool.checkInstallation()))
    logger.info("seconds per event: %d" % tool.secondsPerEvent)
    slhafile = "inputFiles/slha/simplyGluino.slha"
    slhapath = os.path.join ( installation.installDirectory(), slhafile )
    logger.info ( "slhafile: " + slhapath )
    output = tool.run(slhapath, unlink = True )
    #for i in output:
    #    print ( "%s" % i )
    logger.info ( "done: %s" % output )
