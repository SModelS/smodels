#!/usr/bin/env python

"""
.. module:: externalPythia8
   :synopsis: Wrapper for pythia8.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
from smodels.tools.externalTool import ExternalTool
from smodels.tools import externalTool
from smodels.tools.smodelsLogging import logger
from smodels import installation
import os
try:
    import commands as executor
except ImportError:
    import subprocess as executor
import urllib
import tarfile

class ExternalPythia8(ExternalTool):
    """
    An instance of this class represents the installation of pythia8.
    
    """
    def __init__(self,
                 configFile="<install>/lib/pythia8/pythia8.cfg",
                 executablePath="<install>/lib/pythia8/pythia8.exe",
                 srcPath="<install>/lib/pythia8/"):
        """ 
        :param configFile: Location of the config file, full path; copy this
        file and provide tools to change its content and to provide a template
        :param executablePath: Location of executable, full path (pythia_lhe)
        
        nevents - Keep track of how many events we run over for each event we
        only allow a certain computation time if
        self.secondsPerEvent * self.nevents > CPU time, we terminate Pythia.
        
        """
        ExternalTool.__init__(self)
        self.name = "pythia8"
        self.executablePath = self.absPath(executablePath)
        self.executable = None
        self.srcPath = self.absPath(srcPath)
        self.tempdir = None
        self.cfgfile = self.checkFileExists(configFile)
        self.keepTempDir = False
        self.nevents = None
        self.secondsPerEvent = 10

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
            logger.warn("Keeping everything in " + self.tempdir)
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


    def complain ( self ):
        import sys
        logger.error("please fix manually, e.g. try 'make' in smodels/lib, " \
               " or file a complaint at smodels-users@lists.oeaw.ac.at" )
        sys.exit(0)

    def run(self, slhaFile, cfgfile=None, do_unlink=True, do_compile=False,
            do_check=True ):
        """
        Run Pythia.
        
        :param slhaFile: SLHA file
        :param cfgfile: optionally supply a new config file; if not supplied,
                        use the one supplied at construction time; 
                        this config file will not be touched or copied;  
                        it will be taken as is
        :param do_unlink: clean up temporary files after run?
        :param do_compile: if true, we try to compile binary if it isnt installed.
        :param do_check: check installation, before running 
        :returns: stdout and stderr, or error message
        
        """
        if do_check:
            ci=self.checkInstallation()
            if not ci:
                if not do_compile:
                    logger.error("couldnt find pythia8 binary." )
                    self.complain()
                logger.warning("couldnt find pythia8 binary. I have been asked to try to compile it, though. Lets see.")
                # self.compile()
                self.complain()
        slha = self.checkFileExists(slhaFile)
        cfg = self.absPath(cfgfile)
        logger.debug("running with " + str(cfg))
        import shutil
        shutil.copy(slha, self.tempDirectory() + "/fort.61")
        cmd = "cd %s ; %s < %s" % \
             (self.tempDirectory(), self.executablePath, cfg)
        logger.debug("Now running " + str(cmd))
        out = executor.getoutput(cmd)
        # out = subprocess.check_output ( cmd, shell=True, universal_newlines=True )
        if do_unlink:
            self.unlink( unlinkdir=True )
        else:
            f = open(self.tempDirectory() + "/log", "w")
            f.write (cmd + "\n\n\n")
            f.write (out + "\n")
            f.close()
        return out

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


    def compile(self):
        """
        Compile pythia_lhe.
        
        """
        logger.info("Trying to compile pythia in %s" % self.srcPath )
        cmd = "cd %s; make" % self.srcPath
        outputMessage = executor.getoutput(cmd)
        #outputMessage = subprocess.check_output ( cmd, shell=True, 
        #                                          universal_newlines=True )
        logger.info(outputMessage)


    def fetch(self):
        """
        Fetch and unpack tarball.
        
        """
        tempFile = "/tmp/pythia.tar.gz"
        fileHandle = open(tempFile, "w")
        logger.debug("Fetching tarball...")
        url = "http://smodels.hephy.at/externaltools/pythia/pythia.tar.gz"
        link = urllib.urlopen(url)
        lines = link.readlines()
        for line in lines:
            fileHandle.write(line)
        link.close()
        fileHandle.close()
        logger.debug("... done.")
        logger.debug("Untarring...")
        tar = tarfile.open(tempFile)
        for item in tar:
            tar.extract(item, self.srcPath)
        logger.debug("... done.")


    def checkInstallation(self, fix=False ):
        """
        Check if installation of tool is correct by looking for executable and
        running it.

        :param fix: should it try to fix the situation, if something is wrong?

        :returns: True, if everything is ok
        
        """
        if not os.path.exists(self.executablePath):
            logger.error("Executable '%s' not found. Maybe you didn't compile " \
                         "the external tools in smodels/lib?", self.executablePath)
            if fix:
                self.compile()
            return False
        if not os.access(self.executablePath, os.X_OK):
            logger.error("%s is not executable", self.executable)
            self.chmod()
            return False
        slhaFile = "/inputFiles/slha/gluino_squarks.slha"
        slhaPath = installation.installDirectory() + slhaFile
        try:
            output = self.run(slhaPath, "<install>/lib/pythia8/pythia8.cfg",
                    do_compile=False, do_check=False ) 
            output = output.split("\n")
            print ( "output0=",output[0] )
        except Exception as e:
            logger.error("Something is wrong with the setup: exception %s", e)
            return False
        return True


if __name__ == "__main__":
    tool = ExternalPythia8()
    print("installed: " + str(tool.installDirectory()))
    td_exists = externalTool.ok(os.path.exists(tool.tempDirectory()))
    print("temporary directory: %s: %s" % (str(tool.tempDirectory()),
                                           td_exists))
    print("check: " + externalTool.ok(tool.checkInstallation()))
    print("seconds per event: %d" % tool.secondsPerEvent)
    tool.replaceInCfgFile({"NEVENTS": 1, "SQRTS":8000})
    tool.setParameter("MSTP(163)", "6")
    slhafile = "inputFiles/slha/simplyGluino.slha"
    slhapath = os.path.join ( installation.installDirectory(), slhafile )
    print ( "slhafile=", slhapath )
    output = tool.run(slhapath)
    isok = (len (output.split("\n")) > 570)
    print("run: " + externalTool.ok (isok))
    # tool.unlink()
