#!/usr/bin/env python3

"""
.. module:: wrapperBase
   :synopsis: Wrapper code for external tools: base class

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os
import inspect
from smodels.tools.smodelsLogging import logger
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
            logger.error("%s doesnt exist" % self.executablePath )
            return False
        import stat
        mode = stat.S_IRWXU | stat.S_IRWXG | stat.S_IXOTH | stat.S_IROTH
        os.chmod ( self.executablePath, mode )
        return True

    def compile(self):
        """
        Try to compile the tool.
        """
        if not self.maycompile:
            logger.error("Asking to compile, but auto-compilation turned off for %s", self.name )
            return
        logger.debug("Trying to compile %s", self.name)
        cmd = "cd %s; make" % self.srcPath
        out = executor.getoutput(cmd)
        # out = subprocess.check_output ( cmd, shell=True, universal_newlines=True )
        logger.debug(out)
        if not os.path.exists ( self.executablePath ):
            if self.maycompile: ## should have worked
                logger.error ( "Compilation of %s failed. Is the %s compiler installed?" % ( self.name, self.compiler ) )
            sys.exit()
        logger.info ( "Compilation of %s succeeded!" % ( self.name ) )
        return True


    def checkInstallation(self, compile=True ):
        """
        Checks if installation of tool is correct by looking for executable and
        executing it. If check is False and compile is True, then try and compile it.

        :returns: True, if everything is ok

        """
        if not os.path.exists(self.executablePath):
            if compile:
                logger.warning("%s executable not found. Trying to compile it now. This may take a while." % self.name )
                self.compile()
            else:
                logger.warning("%s executable not found." % self.name )
                self.complain()
                return False
        if not os.path.exists(self.executablePath):
            if self.maycompile: ## should have worked
                logger.error("Compilation of %s failed. Is a according compiler installed?" % self.name )
            self.complain()
        if not os.access(self.executablePath, os.X_OK):
            logger.warning("%s is not executable Trying to chmod" % self.executable)
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
        if path == None:
            return self.tempDirectory() + "/temp.cfg"
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
