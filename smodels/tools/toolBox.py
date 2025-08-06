#!/usr/bin/env python3

"""
.. module:: toolBox
   :synopsis: Contains a singleton-like class that keeps track of all external
      "HEP" tools, such as pythia, nllfast, etc. 
      Used primarily for installation and deployment.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from smodels.tools import pythia6Wrapper
from smodels.tools import pythia8Wrapper
from smodels.tools import nllFastWrapper
from smodels.tools import externalPythonTools
from smodels.base.smodelsLogging import logger
from smodels.base.smodelsLogging import colors

class ToolBox(object):
    """
    A singleton-like class that keeps track of all external tools.
    Intended to make installation and deployment easier.
    
    """
    __shared_state = {"tools" : {}}


    def __init__(self):
        """
        Constructor creates the singleton.
        """
        # instead of making this a singleton, we introduce
        self.__dict__ = self.__shared_state
        if len(self.__shared_state["tools"]) == 0:
            self.initSingleton()


    def initSingleton(self):
        """
        Initializes singleton instance (done only once for the entire class).
        
        """
        self.add(pythia6Wrapper.Pythia6Wrapper())
        self.add(pythia8Wrapper.Pythia8Wrapper())
        for tool in nllFastWrapper.nllFastTools.values():
                self.add(tool)
        for tool in externalPythonTools.pythonTools.values():
                self.add(tool)
        

    def add(self, instance):
        """
        Adds a tool by passing an instance to this method.
        
        """
        self.tools[instance.name] = instance


    def listOfTools(self):
        """
        Returns a simple list with the tool names.
        
        """
        return self.tools.keys()


    def installationOk(self, ok ):
        """
        Returns color coded string to signal installation issues.
        """

        if ok == True:
            ret = f"{colors.green}installation ok!{colors.reset}"
            return ret
        ret = f"{colors.red}problem with installation"
        if type(ok) == str:
            ret += f" ({ok})"
        ret += colors.reset
        return ret


    def checkInstallation(self, make=False, printit=True, longL=False ):
        """
        Checks if all tools listed are installed properly, 
        returns True if everything is ok, False otherwise.
        """

        ret = "%sThe following tools have been found in the Toolbox:%s\n" % \
               ( colors.yellow, colors.reset )
        hasMade = False
        allOk=True
        maxl = 45
        if longL: maxl=75
        for(name, instance) in self.tools.items():
            ok = instance.checkInstallation()
            if not ok:
                allOk=False
            exe = instance.pathOfExecutable()
            if len(exe) > maxl + 4:
                exe = "... " + instance.pathOfExecutable()[-maxl:]
            ret += ( "%-12s [%-"+str(maxl+5)+"s]:  %s\n" ) % (name, exe,
                                            self.installationOk(ok ))
            if not ok and make:
                instance.compile()
                hasMade = True
        if make and hasMade:
            ret += "Check again:\n"
            r = self.checkInstallation(make=False, printit=False)
            ret += str(r)
            return r
        if printit:
            print (ret)
        return allOk


    def compile(self):
        """
        Tries to compile and install tools that are not yet marked
        as 'installed'.
        """
        for(name, instance) in self.tools.items():
            installOk = instance.checkInstallation()
            if installOk == True:
                continue
            logger.info("Installation of " + str(name) + " not correct. \
                        Trying to compile.")
            instance.compile()


    def get(self, tool, verbose=True):
        """
        Gets instance of tool from the toolbox.
        
        """
        if not tool in self.tools:
            if verbose:
                logger.error(f"Asking for non-existent tool ``{tool}''")
            return None
        return self.tools[tool]

def main ( args ):
    tmp = ToolBox()
    if args.make:
        tmp.compile()
    if args.colors:
        colors.on = True
    tmp.checkInstallation( printit=True, longL = args.long )
