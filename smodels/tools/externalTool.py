#!/usr/bin/env python

"""
.. module:: externalTool
   :synopsis: Wrapper code for external tools: base class

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os
import inspect
from smodels import installation


class ExternalTool(object):
    """
    An instance of this class represents the installation of an external tool.
    
    An external tool encapsulates a tool that is executed via
    commands.getoutput. It defines how the tool is tested for proper
    installation and how the tool is executed.
    
    """
    def __init__(self):
        self.executablePath = ""
        self.tempdir = ""


    def installDirectory(self):
        """
        :returns: the installation directory of the tool
        
        """
        t = self.executablePath
        p = t.rfind("/")
        if p == -1:
            return ""
        return self.executablePath[:p]


    def pathOfExecutable(self):
        """
        :returns: path of executable
        
        """
        return self.executablePath


    def basePath(self):
        """
        Get the base installation path.
        
        """
        return os.path.dirname(inspect.getabsfile(self.basePath))


    def absPath(self, path):
        """
        Get the absolute path of <path>, replacing <install> with the
        installation directory.
        
        """
        if path == None:
            return self.tempdir + "/temp.cfg"
        installdir = installation.installDirectory()
        path = path.replace("<install>", installdir)
        path = path.replace(".egg/smodels", ".egg/")
        path = os.path.abspath(path)
        return path


def ok(b):
    """
    :returns: 'ok' if b is True, else, return 'error'.
    
    """
    if b:
        return "ok"
    return "error"
