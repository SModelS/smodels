#!/usr/bin/env python

"""
.. module:: externalNllFast
   :synopsis: Wrapper for all nllfast versions.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from . import setPath  # pylint: disable=W0611
from .externalTool import ExternalTool
import commands
import os
import logging

logger = logging.getLogger(__name__)


class ExternalNllFast(ExternalTool):
    """
    An instance of this class represents the installation of nllfast.
    
    """
    def __init__(self, sqrts, nllfast_version, test_params, test_condition):
        """
        :param sqrts: sqrt of s, in TeV, as an integer,
        :param nllfast_version: version of the nllfast tool
        :param test_params: what are the test params we need to run things with?
        :param test_condition: the line that should be the last output line when
        running executable
        :src_path: the path of the source code, for compilation
        
        """
        self.sqrts = int(sqrts)
        self.name = "nllfast%d" % sqrts
        self.nllfast_version = nllfast_version
        path = "<install>/tools/external/nllfast/nllfast-"
        location = path + self.nllfast_version + "/"
        self.cd_path = self.absPath(location)
        self.executable_path = self.cd_path + "/nllfast_%dTeV" % self.sqrts
        self.test_params = test_params
        self.test_condition = test_condition
        self.src_path = self.cd_path
        self.executable = ""


    def compile(self):
        """
        Try to compile nllfast.
        
        """
        logger.info("Trying to compile %s", self.name)
        cmd = "cd %s; make" % self.src_path
        out = commands.getoutput(cmd)
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
            tar.extract(item, self.src_path + "/")


    def unlink(self, inputFile):
        """
        Remove inputFile.out 
        
        """
        return
        # fname = "%s/%s.out" % (self.cd_path, inputFile)
        # if os.path.exists(fname):
        #     os.unlink(fname)


    def run_(self, params):
        """
        Execute nllfast7.
        
        :params params: parameters used (e.g. gg cteq5 .... )
        :returns: stdout and stderr, or error message
        
        """
        cmd = "cd %s; %s %s" % (self.cd_path, self.executable_path, params)
        out = commands.getoutput(cmd)
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
        if not os.path.exists(self.executable_path):
            logger.error("Executable '%s' not found", self.executable_path)
            return False
        if not os.access(self.executable_path, os.X_OK):
            logger.error("%s is not executable", self.executable)
            return False
        out = self.run_(self.test_params)
        if out[-1].find(self.test_condition) == -1:
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
                                 test_params="gg cteq 500 600",
                                 test_condition="500.     600.    0.193E+00  "
                                 "0.450E+00  0.497E+00")


class ExternalNllFast8(ExternalNllFast):
    """
    An instance of this class represents the installation of nllfast 8.
    
    """
    def __init__(self):
        ExternalNllFast.__init__(self, 8, "2.1",
                                 test_params="gg cteq 500 600",
                                 test_condition="500.     600.    0.406E+00  "
                                 "0.873E+00  0.953E+00")


class ExternalNllFast13(ExternalNllFast):
    """
    An instance of this class represents the installation of nllfast 13.
    
    """
    def __init__(self):
        ExternalNllFast.__init__(self, 13, "3.0",
                                 test_params="gg cteq 500 600",
                                 test_condition="500.     600.    0.394E+01  "
                                 "0.690E+01  0.731E+01")


class ExternalNllFast14(ExternalNllFast):
    """
    An instance of this class represents the installation of nllfast 14.
    
    """

    def __init__(self):
        ExternalNllFast.__init__(self, 14, "4.01dcpl",
                                 test_params="gdcpl cteq 500 600",
                                 test_condition="500.    0.235E+02  0.346E+02 "
                                 " 0.362E+02")


class ExternalNllFast33(ExternalNllFast):
    """
    An instance of this class represents the installation of nllfast 33.
    
    """
    def __init__(self):
        ExternalNllFast.__init__(self, 33, "5.01dcpl",
                                 test_params="gdcpl cteq 500 600",
                                 test_condition="500.    0.257E+03  0.383E+03"
                                 "  0.393E+03")


nllFastTools = {7:ExternalNllFast7(),
                8:ExternalNllFast8(),
                13:ExternalNllFast13(),
                14:ExternalNllFast14(),
                33:ExternalNllFast33()}


if __name__ == "__main__":
    for (sqrts, tool) in nllFastTools.items():
        print("%s: installed in %s" % (tool.name, tool.installDirectory()))
