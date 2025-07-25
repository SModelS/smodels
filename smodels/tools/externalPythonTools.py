#!/usr/bin/env python3

"""
.. module:: externalPythonTools
   :synopsis: This module is to check the installation of python tools, 
              i.e. unum, scipy, numpy, pyslha.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from smodels.base.smodelsLogging import logger

class ExternalPythonTool(object):
    """
    An instance of this class represents the installation of a python package.
    As it is python-only, we need this only for installation, not for running
    (contrary to nllfast or pythia).    
    
    """
    def __init__(self, importname, optional=False ):
        """
        Initializes the ExternalPythonTool object. Useful for installation. 
        :params optional: optional package, not needed for core SModelS.
        """
        self.name = importname
        self.optional = optional
        self.pythonPath = ""
        try:
            i = __import__(importname)
            self.pythonPath = i.__file__.replace("/__init__.pyc", "")
        except ImportError as e:
            if optional:
                logger.debug(f"could not find {importname}: {e} (but its not necessary for smodels, so dont worry)")
            else:
                logger.error(f"could not find {importname}: {e}")

    def compile ( self ):
        import sys
        cmd = ["install",self.name]
        if sys.prefix == sys.base_prefix:
            # we are *not* in a virtual env, so add '--user'
            cmd = ["install","--user",self.name]
        try:
            import pip
            pip.main( cmd )
            return
        except (ImportError,AttributeError):
            pass
        try:
            import pip._internal
            pip._internal.main( cmd )
            return
        except (ImportError,AttributeError):
            pass
        try:
            from setuptools.command import easy_install
            easy_install.main(["-U","--user",self.name])
            return
        except (ImportError,AttributeError):
            pass



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
                "scipy": ExternalPythonTool("scipy"),
                "pyhf": ExternalPythonTool("pyhf",optional=True),
                "plotly": ExternalPythonTool("plotly",optional=True),
                "pandas": ExternalPythonTool("pandas",optional=True),
                "ipython": ExternalPythonTool("IPython",optional=True), }


if __name__ == "__main__":
    for (name, tool) in pythonTools.items():
        print(f"{name}: installed in {str(tool.installDirectory())}")
