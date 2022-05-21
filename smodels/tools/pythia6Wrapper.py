#!/usr/bin/env python3

"""
.. module:: Pythia6Wrapper
   :synopsis: Wrapper for pythia6.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
from smodels.tools.wrapperBase import WrapperBase
from smodels.tools import wrapperBase
from smodels.tools.physicsUnits import pb, TeV
from smodels.tools.smodelsLogging import logger
from smodels.theory import crossSection
from smodels.theory.crossSection import LO
from smodels import installation
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
import os, sys, io

try:
    import commands as executor
except ImportError:
    import subprocess as executor
import urllib
import tarfile

class Pythia6Wrapper(WrapperBase):
    """
    An instance of this class represents the installation of pythia6. nevents keeps track of
    how many events we run. For each event we only allow a certain computation time:
    if self.secondsPerEvent * self.nevents > CPU time, we terminate Pythia.
    """

    def __init__(self,
                 configFile="<install>/smodels/etc/pythia.card",
                 executablePath="<install>/smodels/lib/pythia6/pythia_lhe",
                 srcPath="<install>/smodels/lib/pythia6/"):
        """
        :param configFile: Location of the config file, full path; copy this
                           file and provide tools to change its content and to provide a template
        :param executablePath: Location of executable, full path (pythia_lhe)
        :param srcPath: Location of source code
        """

        WrapperBase.__init__(self)
        self.name = "pythia6"
        self.executablePath = self.absPath(executablePath)
        self.executable = None
        self.srcPath = self.absPath(srcPath)
        self.compiler = "gfortran"
        self.tempdir = None
        self.cfgfile = self.checkFileExists(configFile)
        self.keepTempDir = False
        self.nevents = 10000
        self.sqrts = 8 # sqrt(s) in TeV
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
        logger.debug( "Unlinking " + self.tempdir )
        for inputFile in ["fort.61", "fort.68", "log"]:
            if os.path.exists(self.tempdir + "/" + inputFile):
                os.unlink(self.tempdir + "/" + inputFile)
        if unlinkdir:
            for inputFile in ["temp.cfg"]:
                os.unlink(self.tempdir + "/" + inputFile)
            if os.path.exists(self.tempdir):
                os.rmdir(self.tempdir)
                self.tempdir = None


    def replaceInCfgFile(self, replacements={"NEVENTS": 10000, "SQRTS":8000}):
        """
        Replace strings in the config file by other strings, similar to
        setParameter.

        This is introduced as a simple mechanism to make changes to the
        parameter file.

        :param replacements: dictionary of strings and values; the strings will
                             be replaced with the values; the dictionary keys
                             must be strings present in the config file

        """
        f = open(self.tempDirectory() + "/temp.cfg")
        lines = f.readlines()
        f.close()
        f = open(self.tempDirectory() + "/temp.cfg", "w")
        for line in lines:
            # print ( "line=%s" % line )
            for (key, value) in replacements.items():
                if key == "NEVENTS":
                    self.nevents = int(value)
                line = line.replace(key, str(value))
            f.write(line)
        f.close()


    def setParameter(self, param="MSTP(163)", value=6):
        """
        Modifies the config file, similar to .replaceInCfgFile.

        It will set param to value, overwriting possible old values.

        """
        f = open(self.tempdir + "/temp.cfg")
        lines = f.readlines()
        f.close()
        f = open(self.tempdir + "/temp.cfg", "w")
        for line in lines:  # # copy all but lines with "param"
            if not param in line:
                f.write(line)
        f.write("%s=%s\n" % (param, str(value)))
        f.close()

    def run( self, slhafile, lhefile=None, unlink=True ):
        """
        Execute pythia_lhe with n events, at sqrt(s)=sqrts.

        :param slhafile: input SLHA file
        :param lhefile: option to write LHE output to file; if None, do not write
                        output to disk. If lhe file exists, use its events for
                        xsecs calculation.
        :param unlink: Clean up temp directory after running pythia

        :returns: List of cross sections

        """
        if lhefile and os.path.isfile ( lhefile ):
            lheFile = open(lhefile, 'r')
            xsecsInfile = crossSection.getXsecFromSLHAFile(slhafile)
            loXsecs = crossSection.XSectionList()
            for xsec in xsecsInfile:
                if xsec.info.order == LO and \
                        float (xsec.info.sqrts.asNumber(TeV)) == self.sqrts:
                    loXsecs.add(xsec)
            return loXsecs

        #Change pythia card, if defined:
        if self.pythiacard:
            pythiacard_default = self.cfgfile
            self.cfgfile = self.pythiacard
        # Check if template config file exists
        if unlink:
            self.unlink()
        else:
            self.tempdir = None
        self.replaceInCfgFile({"NEVENTS": self.nevents, "SQRTS":1000 * self.sqrts})
        self.setParameter("MSTP(163)", "6")

        if unlink==False:
            logger.info ( "keeping temporary directory at %s" % self.tempDirectory() )
        r = self.checkInstallation()
        if r == False:
            logger.info ( "Installation check failed." )
            sys.exit()
        self.replaceInCfgFile({"NEVENTS": self.nevents, "SQRTS":1000 * self.sqrts})
        self.setParameter("MSTP(163)", "6")
        lhedata = self._run(slhafile, unlink=unlink )
        if not "<LesHouchesEvents" in lhedata:
            pythiadir = "%s/log" % self.tempDirectory()
            logger.error("No LHE events found in pythia output %s" % pythiadir )
            if not os.path.exists ( pythiadir ):
                logger.error ("Will dump pythia output to %s" % pythiadir )
                f=open ( pythiadir, "w" )
                for line in lhedata:
                    f.write ( line )
                f.close()
            raise SModelSError( "No LHE events found in %s" % pythiadir )

        #Reset pythia card to its default value
        if self.pythiacard:
            self.cfgfile = pythiacard_default

        #if not unlink:
        #    lhefile = self.tempdir + "/events.lhe"
        # Generate file object with lhe events
        if lhefile:
            lheFile = open(lhefile, 'w')
            lheFile.write(lhedata)
            lheFile.close()
            lheFile = open(lhefile, 'r')
        else:
            # Create memory only file object
            if sys.version[0]=="2":
                lhedata = unicode ( lhedata )
            lheFile = io.StringIO(lhedata)
        return crossSection.getXsecFromLHEFile ( lheFile )


    def _run(self, slhaFile, cfgfile=None, unlink=True, do_compile=False,
            do_check=True ):
        """
        Really Run Pythia.

        :param slhaFile: SLHA file
        :param cfgfile: optionally supply a new config file; if not supplied,
                        use the one supplied at construction time;
                        this config file will not be touched or copied;
                        it will be taken as is
        :param unlink: clean up temporary files after run?
        :param do_compile: if true, we try to compile binary if it isnt installed.
        :param do_check: check installation, before running
        :returns: stdout and stderr, or error message

        """
        if do_check:
            self.checkInstallation( do_compile )
        slha = self.checkFileExists(slhaFile)
        cfg = self.absPath(cfgfile)
        ck_cfg = self.checkFileExists(cfg)
        logger.debug("running with " + str(cfg))
        import shutil
        shutil.copy(slha, self.tempDirectory() + "/fort.61")
        cmd = "cd %s ; %s < %s" % \
             (self.tempDirectory(), self.executablePath, cfg)
        logger.debug("Now running " + str(cmd))
        out = executor.getoutput(cmd)
        if unlink:
            self.unlink( unlinkdir=True )
        else:
            f = open(self.tempDirectory() + "/log", "w")
            f.write (cmd + "\n\n\n")
            f.write (out + "\n")
            f.close()
        return out

if __name__ == "__main__":
    #from smodelsLogging import setLogLevel
    #setLogLevel ( "debug" )
    tool = Pythia6Wrapper()
    print("[Pythia6Wrapper] installed: " + str(tool.installDirectory()))
    td_exists = wrapperBase.ok(os.path.exists(tool.tempDirectory()))
    print("[Pythia6Wrapper] temporary directory: %s: %s" % (str(tool.tempDirectory()),
                                           td_exists))
    print("[Pythia6Wrapper] check: " + wrapperBase.ok(tool.checkInstallation()))
    print("[Pythia6Wrapper] seconds per event: %d" % tool.secondsPerEvent)
    tool.replaceInCfgFile({"NEVENTS": 1, "SQRTS":8000})
    tool.setParameter("MSTP(163)", "6")
    slhafile = "inputFiles/slha/simplyGluino.slha"
    slhapath = os.path.join ( installation.installDirectory(), slhafile )
    # print ( "[Pythia6Wrapper] slhapath:",slhapath )
    xsec = tool.run(slhapath, unlink=True )
    # n = len (output.split("\n"))
    # print ( "n: ",n )
    isok = abs ( xsec[0].value.asNumber ( pb ) - 2.80E-01  ) < 1e-5
    print("[Pythia6Wrapper] run: " + wrapperBase.ok (isok))
    # tool.unlink()
