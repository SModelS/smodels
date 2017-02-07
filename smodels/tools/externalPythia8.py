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
from smodels.tools.physicsUnits import fb, pb, TeV, mb
from smodels.theory.crossSection import XSection, XSectionInfo, LO, XSectionList
from smodels import installation
from smodels.tools.pythia8particles import particles
import os, sys

try:
    import commands as executor #python2 
except ImportError:
    import subprocess as executor # python3
    
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
        self.xsecs = {}
        logger.debug ( "started .run" )
        slha = self.checkFileExists(slhaFile)
        logger.debug ( "file check: " + slha )
        cfg = self.absPath(self.cfgfile)
        logger.debug("running with cfgfile " + str(cfg))
        cmd = "%s -n %d -f %s -s %d -c %s" % \
             ( self.executablePath, self.nevents, slha, self.sqrts, cfg )
        xmldoc = self.executablePath.replace ( "pythia8.exe", "xml.doc" )
        if os.path.exists (xmldoc ):
            logger.info ( "xml.doc found at %s." % xmldoc )
            with open ( xmldoc ) as f:
                xmlDir = f.read()
                logger.debug ( "adding -x %s" % xmlDir )
                cmd += " -x %s" % xmlDir.strip()
        logger.debug("Now running ''%s''" % str(cmd) )
        out = executor.getoutput(cmd)
        if not unlink:
            tempfile = self.tempDirectory() + "/log"
            f = open( tempfile, "w")
            f.write (cmd + "\n\n\n")
            f.write (out + "\n")
            f.close()
            logger.debug ( "stored everything in %s" % tempfile )
        self.parse ( out )
        ret = XSectionList()
        for key,value in self.xsecs.items():
            xsec = XSection()
            xsec.info = XSectionInfo( self.sqrts * TeV, LO, 
                                      "%d TeV (LO)" % self.sqrts )
            xsec.value = value
            xsec.pid = key
            ret.add ( xsec )
        self.xsecs.clear()
        return ret

    def getProcess ( self, token ):
        """ obtain the process identifier from the pythia8 output """
        parrow = token.find ( "->" ) 
        if parrow < 1:
            logger.error ( "cannot parse pythia8 output" )
            return None
        process = token[parrow+2:-4].strip()
        s1 = process.find ( " " )
        p=[ None, None ]
        p[0] = process [ : s1 ]
        s2 = process [ s1+1 : ].find ( " " )
        if s2 == -1:
            s2 = len( process )
        else:
            s2 += s1+1
        p[1] = process [ s1+1 : s2 ]

        """
        print ( "process >>%s<<" % ( process ) )
        print ( "cut off >>%s<<" % process[s1+1 : ] )
        print ( "s2t off >>%s<<" % process[s1+1 : ] )
        print ( "p2 >>%s<<" % p2 )
        """
        for key, value in particles.items():
            particles[key+"bar"]=-value
        ret = [ None, None ]
        for i in [0,1]:
            if p[i] in particles.keys():
                ret[i]=particles[p[i] ]
            else:
                logger.error ( "particle >>%s<< unknown. " \
                        "Check pythia8particles.py." % p[i] )
                sys.exit()
        return tuple ( ret )

    def getXSec ( self, token ):
        """ obtain the cross section from the pythia8 output """
        value, error = map ( float, token.split () )
        # print ( "xsec: >>%s<<" % (value * mb) )
        return (value * mb).asUnit ( fb )

    def parseLine ( self, line ):
        tmp = line.split ( "|" )
        tokens = [ x.strip() for x in tmp[1:-1]  ]
        process = self.getProcess ( tokens[0] )
        xsec = self.getXSec ( tokens[-1] )
        if xsec < 1e-8 * fb:
            return
        if not process in self.xsecs.keys():
            self.xsecs [ process ] = 0. * fb
        self.xsecs[ process ] += xsec

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
            logger.error ( "could not parse pythia8 output file correctly." )
            with open("/tmp/pythia8.log", "w" ) as f:
                f.write ( out )
                logger.error ( "see /tmp/pythia8.log for more details" )
            return None
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
    setLogLevel ( "debug" )
    tool = ExternalPythia8()
    tool.nevents=10
    logger.info("installed: " + str(tool.installDirectory()))
    logger.info("check: " + externalTool.ok(tool.checkInstallation()))
    logger.info("seconds per event: %d" % tool.secondsPerEvent)
    slhafile = "inputFiles/slha/simplyGluino.slha"
    slhapath = os.path.join ( installation.installDirectory(), slhafile )
    logger.info ( "slhafile: " + slhapath )
    output = tool.run(slhapath, unlink = True )
    #for i in output:
    #    print ( "%s" % i )
    logger.info ( "done: %s" % output )
