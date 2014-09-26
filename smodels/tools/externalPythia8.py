#!/usr/bin/env python

"""
.. module:: externalPythia8
   :synopsis: Wrapper for pythia8.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
import ctypes.util
from smodels.tools.externalTool import ExternalTool
from smodels.tools import externalTool
import logging
logger = logging.getLogger(__name__)


class ExternalPythia8(ExternalTool):
    """
    An instance of this class represents the installation of pythia8.
    
    """
    def __init__(self):
        self.name = "pythia8"
        lib = ctypes.util.find_library(self.name)
        self.executablePath = lib


    def compile(self):
        """
        Compile pythia_lhe.
        
        """
        logger.info("Trying to compile pythia8.")


    def checkInstallation(self):
        """
        Check if installation of tool is correct by looking for executable and
        running it.
        
        """
        for i in [ "HepMC", "LHAPDF", "lhapdfdummy", "pythia8" ]:
            try:
                a = ctypes.util.find_library (i)
                if a == None:
                    logger.error("could not find library lib%s.so", i)
                    return False
                d = ctypes.CDLL("lib%s.so" % i, ctypes.RTLD_GLOBAL)
            except OSError, e:
                logger.error("could not load library lib%s.so: %s", i, e)
                return False
        return True


if __name__ == "__main__":
    tool = ExternalPythia8()
    print("installed: " + str(tool.installDirectory()))
    # print("temporary directory: %s: %s" % (str(tool.tempDirectory()),
    #                                      td_exists))
    print("check: " + externalTool.ok(tool.checkInstallation()))
    # print("run:")
    # slhafile = SModelS.installDirectory() + "/inputFiles/slha/andrePT4.slha"
    # out = tool.run(slhafile)
