#!/usr/bin/env python

"""
.. module:: externalPythia8
   :synopsis: Wrapper for pythia8.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
from smodels.tools.externalTool import ExternalTool
from smodels.tools import externalTool
from smodels.tools.smodelsLogging import logger, setLogLevel
from smodels.tools.physicsUnits import fb, pb, TeV
from smodels import installation
import os
try:
    import commands as executor #python2 
except ImportError:
    import subprocess as executor # python3
    
setLogLevel ( "debug" )

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
        self.nevents = 10000
        self.sqrts = 13
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

    def complain ( self ):
        import sys
        logger.error("please fix manually, e.g. try 'make' in smodels/lib, " \
               " or file a complaint at smodels-users@lists.oeaw.ac.at" )
        sys.exit(0)

    def run(self, slhaFile, cfgfile=None, do_unlink=True, do_compile=False,
            do_check=False ):
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
        self.xsecs = {}
        logger.debug ( "started .run" )
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
        logger.debug ( "file check: " + slha )
        if cfgfile != None:
            self.cfgfile = cfgfile
        cfg = self.absPath(self.cfgfile)
        logger.debug("running with cfgfile " + str(cfg))
        cmd = "%s -n %d -f %s -s %d -c %s" % \
             ( self.executablePath, self.nevents, slha, self.sqrts, cfg )
        logger.debug("Now running " + str(cmd))
        out = executor.getoutput(cmd)
        # out = subprocess.check_output ( cmd, shell=True, universal_newlines=True )
        if not do_unlink:
            tempfile = self.tempDirectory() + "/log"
            f = open( tempfile, "w")
            f.write (cmd + "\n\n\n")
            f.write (out + "\n")
            f.close()
            logger.debug ( "stored everything in %s" % tempfile )
        self.parse ( out )
        return self.xsecs

    def getProcess ( self, token ):
        return (100021,100021)

    def getXSec ( self, token ):
        return 3.0 * fb

    def parseLine ( self, line ):
        tmp = line.split ( "|" )
        tokens = [ x.strip() for x in tmp[1:-1]  ]
        process = self.getProcess ( tokens[0] )
        xsec = self.getXSec ( tokens[-1] )
        self.xsecs[ process ] = xsec

    def parse( self, out ):
        lines = out.split ("\n")
        # print ( "Parsing outfile with %d lines." % len(lines) )
        begin,end=1,0
        for ctr,line in enumerate(lines):
            if "--  PYTHIA Event and Cross Section Statistics" in line:
                begin=ctr+7
            if "End PYTHIA Event and Cross Section Statistics" in line:
                end=ctr-3
        if begin > end:
            logger.error ( "could not parse pythia8 file correctly." )
            return None
        # print ( "need to look at %d to %d" % ( begin,end ) )
        for line in lines[begin:end]:
            self.parseLine ( line )

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


    def checkInstallation(self, fix=False ):
        """
        Check if installation of tool is correct by looking for executable

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
        return True

if __name__ == "__main__":
    tool = ExternalPythia8()
    tool.nevents=10
    logger.info("installed: " + str(tool.installDirectory()))
    logger.info("check: " + externalTool.ok(tool.checkInstallation()))
    logger.info("seconds per event: %d" % tool.secondsPerEvent)
    slhafile = "inputFiles/slha/simplyGluino.slha"
    slhapath = os.path.join ( installation.installDirectory(), slhafile )
    logger.info ( "slhafile: " + slhapath )
    output = tool.run(slhapath, do_unlink = False)
    logger.info ( "done: %s" % output )
