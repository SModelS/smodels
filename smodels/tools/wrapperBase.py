#!/usr/bin/env python3

"""
.. module:: wrapperBase
   :synopsis: Wrapper code for external tools: base class

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os
import inspect
from smodels.base.smodelsLogging import logger
from smodels import installation
import sys
try:
    import commands as executor
except ImportError:
    import subprocess as executor

class WrapperBase(object):
    """
    An instance of this class represents the installation of an external tool.

    An external tool encapsulates a tool that is executed via
    commands.getoutput. The wrapper defines how the tool is tested for proper
    installation and how the tool is executed.

    """
    #defaulttempdir = "./" ## the default directory for temp dirs
    defaulttempdir = "/tmp/" ## the default directory for temp dirs

    def __init__(self):
        self.executablePath = ""
        self.tempdir = ""
        self.maycompile = True

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

    def chmod(self):
        """
        chmod 755 on executable, if it exists.  Do nothing, if it doesnt exist.
        """
        if not os.path.exists ( self.executablePath ):
            logger.error(f"{self.executablePath} doesnt exist" )
            return False
        import stat
        mode = stat.S_IRWXU | stat.S_IRWXG | stat.S_IXOTH | stat.S_IROTH
        os.chmod ( self.executablePath, mode )
        return True

    def compile ( self ):
        """
        Try to compile the tool.
        """
        if not self.maycompile:
            logger.error("Asking to compile, but auto-compilation turned off for %s", self.name )
            return
        logger.info( f"Trying to compile {self.name}" )
        cmd = f"cd {self.srcPath}; make"
        out = executor.getoutput(cmd)
        # out = subprocess.check_output ( cmd, shell=True, universal_newlines=True )
        logger.debug(out)
        if not os.path.exists ( self.executablePath ):
            if self.maycompile: ## should have worked
                logger.error ( f"Compilation of {self.name} failed. Is the {self.compiler} compiler installed?" )
            sys.exit()
        logger.info ( f"Compilation of {self.name} succeeded!" )
        return True


    def checkInstallation(self, compile=True ):
        """
        Checks if installation of tool is correct by looking for executable and
        executing it. If check is False and compile is True, then try and compile it.

        :returns: True, if everything is ok

        """
        if not os.path.exists(self.executablePath):
            if compile:
                logger.warning( f"{self.name} executable not found. Trying to compile it now. This may take a while." )
                self.compile()
            else:
                logger.warning( f"{self.name} executable not found." )
                self.complain()
                return False
        if not os.path.exists(self.executablePath):
            if self.maycompile: ## should have worked
                logger.error(f"Compilation of {self.name} failed. Is a according compiler installed?" )
            self.complain()
        if not os.access(self.executablePath, os.X_OK):
            logger.warning(f"{self.executable} is not executable Trying to chmod")
            self.chmod()
        return True

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
        if path is None:
            return os.path.abspath(self.tempDirectory() + "/temp.cfg")
        installdir = installation.installDirectory()
        path = path.replace("<install>", installdir)
        # path = path.replace(".egg/smodels", ".egg/")
        path = os.path.abspath(path)
        return path

    def tempDirectory(self):
        """
        Return the temporary directory name.

        """
        import tempfile
        import shutil
        if self.tempdir in [ None, "" ]:
            self.tempdir = tempfile.mkdtemp( prefix="xsec",
                                             dir=self.defaulttempdir )
            # self.tempdir = tempfile.mkdtemp()
            shutil.copy(self.cfgfile, self.tempdir + "/temp.cfg")
        return self.tempdir

    def complain ( self ):
        logger.error("please fix manually, e.g. try 'make' in smodels/lib, " \
               " or file a complaint at smodels-users@lists.oeaw.ac.at" )
        sys.exit(0)



def ok(b):
    """
    :returns: 'ok' if b is True, else, return 'error'.

    """
    if b:
        return "ok"
    return "error"
