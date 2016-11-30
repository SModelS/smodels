#!/usr/bin/env python

"""
.. module:: externalPythonTools
   :synopsis: This module is to check the installation of python tools, 
              i.e. unum, scipy, numpy, pyslha.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from smodels.tools.smodelsLogging import logger

class ExternalPythonTool(object):
    """
    An instance of this class represents the installation of unum.
    As it is python-only, we need this only for installation,
    not for running (contrary to nllfast or pythia).    
    
    """
    def __init__(self, importname):
        """
        Initializes the ExternalPythonTool object. Useful for installation. 
        """
        self.name = importname
        self.pythonPath = ""
        try:
            i = __import__(importname)
            self.pythonPath = i.__file__.replace("/__init__.pyc", "")
        except ImportError as e:
            logger.error("could not find %s: %s" % (importname, e))


    def pathOfExecutable (self):
        """
        Just returns the pythonPath variable
        """
        return self.pythonPath


    def installDirectory(self):
        """
        Just returns the pythonPath variable
        """
        return self.pythonPath


    def checkInstallation(self):
        """
        The check is basically done in the constructor
        """
        if self.pythonPath == "":
            return False
        return True


pythonTools = { "unum" : ExternalPythonTool("unum"),
                "numpy": ExternalPythonTool("numpy"),
                "pyslha": ExternalPythonTool("pyslha"),
                "scipy": ExternalPythonTool("scipy") }


if __name__ == "__main__":
    for (name, tool) in pythonTools.items():
        print("%s: installed in %s" % (name, str(tool.installDirectory())))
