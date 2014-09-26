#!/usr/bin/env python

"""
.. module:: externalPythonTools
   :synopsis: this module is to check the installation of python tools, 
              i.e. unum, scipy, numpy.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os
import inspect
import logging

logger = logging.getLogger(__name__)

class ExternalPythonTool(object):
    """
    An instance of this class represents the installation of unum.
    As it is python-only, we need this only for installation,
    not for running (contrary to nllfast or pythia).    
    """
    def __init__(self, importname ):
        """
        Initializes the ExternalPythonTool object. Useful for installation. 
        """
        self.name=importname
        self.python_path = ""
        try:
            i=__import__(importname)
            self.python_path = i.__file__.replace("/__init__.pyc","")
        except ImportError,e:
            logger.error("could not find %s: %s" % (importname,e) )

    def pathOfExecutable ( self ):
        """
        Just returns the python_path variable
        """
        return self.python_path

    def installDirectory(self):
        """
        Just returns the python_path variable
        """
        return self.python_path

    def checkInstallation(self):
        """
        The check is basically done in the constructor
        """
        if self.python_path == "":
            return False
        return True

pythonTools = { "unum": ExternalPythonTool("unum"), 
                "numpy": ExternalPythonTool("numpy"),
                "scipy": ExternalPythonTool("scipy") }

if __name__ == "__main__":
    import setPath
    from smodels import tools # so we have a logging handler
    for (name,tool) in pythonTools.items():
        print("%s: installed in %s" % ( name, str(tool.installDirectory())) )
