#!/usr/bin/env python

"""
.. module:: externalNllFast
   :synopsis: Wrapper for all nllfast versions.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
try:
    import commands as executor
except ImportError:
    import subprocess as executor
import os
from smodels.tools.externalTool import ExternalTool
from smodels.tools.smodelsLogging import logger


class ExternalNllFast(ExternalTool):
    """
    An instance of this class represents the installation of nllfast.
    
    """
    def __init__(self, sqrts, nllfastVersion, testParams, testCondition):
        """
        :param sqrts: sqrt of s, in TeV, as an integer,
        :param nllfastVersion: version of the nllfast tool
        :param testParams: what are the test params we need to run things with?
        :param testCondition: the line that should be the last output line when
        running executable
        :srcPath: the path of the source code, for compilation
        
        """
        ExternalTool.__init__(self)
        self.sqrts = int(sqrts)
        self.name = "nllfast%d" % sqrts
        self.nllfastVersion = nllfastVersion
        path = "<install>/lib/nllfast/nllfast-"
        location = path + self.nllfastVersion + "/"
        self.cdPath = self.absPath(location)
        self.executablePath = self.cdPath + "/nllfast_%dTeV" % self.sqrts
        self.testParams = testParams
        self.testCondition = testCondition
        self.srcPath = self.cdPath
        self.executable = ""


    def compile(self):
        """
        Try to compile nllfast.
        
        """
        logger.info("Trying to compile %s", self.name)
        cmd = "cd %s; make" % self.srcPath
        out = executor.getoutput(cmd)
        # out = subprocess.check_output ( cmd, shell=True, universal_newlines=True )
        logger.info(out)
        return True


    def fetch(self):
        """
        Fetch and unpack tarball.
        
        """
        import urllib, tarfile
        tempfile = "/tmp/nllfast7.tar.gz"
        f = open(tempfile, "w")
        url = "http://smodels.hephy.at/externaltools/nllfast%d.tar.gz" \
                % self.sqrts
        logger.info("fetching tarball from " + url)
        R = urllib.urlopen(url)
        l = R.readlines()
        for line in l:
            f.write(line)
        R.close()
        f.close()
        tar = tarfile.open(tempfile)
        for item in tar:
            tar.extract(item, self.srcPath + "/")


    def unlink(self, inputFile):
        """
        Remove inputFile.out 
        
        """
        return
        # fname = "%s/%s.out" % (self.cdPath, inputFile)
        # if os.path.exists(fname):
        #     os.unlink(fname)


    def run_(self, params):
        """
        Execute nllfast7.
        
        :params params: parameters used (e.g. gg cteq5 .... )
        :returns: stdout and stderr, or error message
        
        """
        cmd = "cd %s; %s %s" % (self.cdPath, self.executablePath, params)
        out = executor.getoutput(cmd)
        # out = subprocess.check_output ( cmd, shell=True, universal_newlines=True )
        out = out.split("\n")
        return out


    def run(self, process, pdf, squarkmass, gluinomass):
        """
        Execute nllfast.
        
        :params process: which process: st, sb, gg, gdcpl, sdcpl, ss, sg, tot
        :params pdf: cteq=cteq6, mstw2008 
        :params squarkmass: squarkmass, None if squark decoupled
        :params gluinomass: gluinomass, None if gluino decoupled
        :returns: stdout and stderr, or error message
        
        """
        processes = ["st", "sb", "gg", "gdcpl", "sdcpl", "ss", "sg", "tot"]
        if not process in processes:
            return None
        if not pdf in ["cteq", "cteq6", "mstw", "mstw2008"]:
            return None
        if not squarkmass:
            return self.run_("%s %s %s") % (process, pdf, gluinomass)
        if not gluinomass:
            return self.run_("%s %s %s") % (process, pdf, squarkmass)
        return self.run_("%s %s %s %s") % \
                (process, pdf, squarkmass, gluinomass)


    def checkInstallation(self):
        """
        Checks if installation of tool is valid by looking for executable and
        executing it.
        
        """
        if not os.path.exists(self.executablePath):
            logger.error("Executable '%s' not found. Maybe you didn't compile " \
                         "the external tools in smodels/lib?", self.executablePath)
            return False
        if not os.access(self.executablePath, os.X_OK):
            logger.error("%s is not executable", self.executable)
            return False
        out = self.run_(self.testParams)
        if out[-1].find(self.testCondition) == -1:
            logger.error("Setup invalid: " + str(out))
            return False
        self.unlink("gg")
        return True


class ExternalNllFast7(ExternalNllFast):
    """
    An instance of this class represents the installation of nllfast 7.
    
    """
    def __init__(self):
        ExternalNllFast.__init__(self, 7, "1.2",
                                 testParams="gg cteq 500 600",
                                 testCondition="500.     600.    0.193E+00  "
                                 "0.450E+00  0.497E+00")


class ExternalNllFast8(ExternalNllFast):
    """
    An instance of this class represents the installation of nllfast 8.
    
    """
    def __init__(self):
        ExternalNllFast.__init__(self, 8, "2.1",
                                 testParams="gg cteq 500 600",
                                 testCondition="500.     600.    0.406E+00  "
                                 "0.873E+00  0.953E+00")

class ExternalNllFast13(ExternalNllFast):
    """
    An instance of this class represents the installation of nllfast 8.
    
    """
    def __init__(self):
        ExternalNllFast.__init__(self, 13, "3.1",
                                 testParams="gg cteq 500 600",
                                 testCondition="600.    0.394E+01  0.690E+01  "
                                 "0.731E+01    0.394E+00" )

nllFastTools = { 7 : ExternalNllFast7(),
                 8 : ExternalNllFast8(),
                13 : ExternalNllFast13() }


if __name__ == "__main__":
    for (sqrts, tool) in nllFastTools.items():
        print("%s: installed in %s" % (tool.name, tool.installDirectory()))
