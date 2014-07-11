#!/usr/bin/env python

"""
.. module:: externalTool
   :synopsis: Wrapper code for external tools: base class

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os
import inspect

class ExternalTool(object):
    """
    An instance of this class represents the installation of an external tool.
    
    An external tool encapsulates a tool that is executed via
    commands.getoutput. It defines how the tool is tested for proper
    installation and how the tool is executed.
    
    """
    def __init__(self):
        self.executable_path = ""


    def installDirectory(self):
        """
        TODO: write docstring
        
        """
        T = self.executable_path
        P = T.rfind("/")
        if P == -1:
            return ""
        return self.executable_path[:P]


    def pathOfExecutable(self):
        """
        :returns: path of executable
        
        """
        return self.executable_path


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
        if None:
            return None
        installdir = self.basePath()
        installdir = installdir.replace("tools", "")
        path = path.replace("<install>", installdir)
        path = os.path.abspath(path)
        return path

    def ok(self,B):
        """
        return 'ok', B is True. Else, return 'error'
        
        """
        if B:
            return "ok"
        return "error"
