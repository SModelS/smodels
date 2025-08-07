#!/usr/bin/env python3

"""
.. module:: runtime
    :synopsis: Tools to gather info about the runtime environment,
               ( nCPUs() ), or obtain file type ( filetype() ). Pointer
               to model file is also kept here.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from typing import Union, Text, Dict
import os

## place to keep the pointer to the model file (default = mssm)
modelFile="smodels.share.models.mssm"

_experimental = { "truncatedgaussians": False,
                  "spey": False } ## experimental features

_deltas_rel_default = .2 ## the default relative error on the signal strength

def checkForIncompatibleModuleVersions() -> bool:
    """ for 3.1.0, at release time, we have the problem that scipy 1.16
    does not work with pyhf 0.7.6. we warn about such situations.

    :returns: false if incompatible version found, else true
    """
    from importlib.metadata import version
    if version("scipy")[:4] in [ "1.16", "1.17" ]:
        print ( f"SModelS warning: scipy v{version('scipy')} is installed, but it seems " \
                 "this version does not work with pyhf 0.7.6. We recommend scipy 1.15.x " \
                 "and pyhf 0.7.6. You have been warned." )
        return False
    return True
        
checkForIncompatibleModuleVersions()

def returnGitCommitId():
    """ if we are in a git-managed repository, return the git commit id """
    try:
        import subprocess
        # Check if we're inside a Git repository
        subprocess.run(
            ['git', 'rev-parse', '--is-inside-work-tree'],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True
        )

        # Get the current commit hash
        commit = subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'],
            stderr=subprocess.DEVNULL
        ).decode('utf-8').strip()

        return commit

    except subprocess.CalledProcessError:
        # Not a Git repo, or git command failed—do nothing
        return None

def printEnvironmentInfo( args : Dict ):
    """ very simple method that prints out info relevant to debugging
        machine-dependent problems
    :param args: a dictionary with options, currently only a "colors" option
    is available
    :returns: true if all depenendencies are met
    """
    from smodels.base.smodelsLogging import colors
    colors.on = True if "colors" in args and args["colors"] == True else False
    import importlib, platform

    modules = [ "scipy", "sympy", "numpy",
        "pyslha", "unum", "pyhf" ]

    print("Environment Information:")
    print(f"Operating System: {colors.green}{platform.system()} {platform.release()}{colors.reset}")
    print(f"Python Version: {colors.green}{platform.python_version()}{colors.reset}")
    print(f"Machine Architecture: {colors.green}{platform.machine()}{colors.reset}")
    print(f"Processor: {colors.green}{platform.processor()}{colors.reset}")
    print("\nModule Versions:")

    depsMet = True
    for module_name in modules:
        try:
            module = importlib.import_module(module_name)
            version = getattr(module, '__version__', 'Unknown version attribute')
            print(f"{module_name:<12}: {colors.green}{version}{colors.reset}")
        except ImportError:
            print(f"{module_name:<12}: Not installed")
            depsMet = False
    gitId = returnGitCommitId()
    if gitId != None:
        print ( f"\ngit commit: {gitId}" )
    return depsMet

def filetype ( filename : os.PathLike ) -> Union[Text,None]:
    """ obtain information about the filetype of an input file,
        currently only used to discriminate between slha and lhe
        files.

        :returns: filetype as string("slha" or "lhe"),
                  None if file does not exist, or filetype is unknown.
    """
    import os
    if not os.path.exists ( filename ):
        return None
    if filename.endswith(".slha"):
        return "slha"
    if filename.endswith(".SLHA"):
        return "slha"
    if filename.endswith(".lhe" ):
        return "lhe"
    if filename.endswith(".LHE" ):
        return "lhe"
    try:
        with open ( filename, "rt" ) as f:
            for line in f:
                if "<LesHouchesEvents" in line:
                    return "lhe"
                if "<event>" in line:
                    return "lhe"
                if "block " in line.lower():
                    return "slha"
    except UnicodeDecodeError:
        ## a binary file??
        return None
    return None

def experimentalFeature( feature : str ) -> Union[None,bool]:
    """ method to check if a certain experimental feature is enabled.
    can be turned on and off via options:experimental in parameters.ini.
    :param feature: ask for feature

    :returns: None if feature does not exist, else boolean
    """
    if not feature in _experimental:
        return None
    return _experimental[feature]

def nCPUs():
    """ obtain the number of *available* CPU cores on the machine, for several
        platforms and python versions. """
    try:
        # next few lines taken from
        # https://stackoverflow.comhttps//stackoverflow.com/questions/1006289/how-to-find-out-the-number-of-cpus-using-python/questions/1006289/how-to-find-out-the-number-of-cpus-using-python
        import re
        with open('/proc/self/status') as f:
            m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$', f.read())
        if m:
            res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
            if res > 0:
                return res
    except IOError:
        pass
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except ImportError:
        pass
    try:
        import psutil
        return psutil.NUM_CPUS
    except ImportError:
        pass
    try:
        import os
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))
        if res>0: return res
    except ImportError:
        pass
    return None

if __name__ == "__main__":
    printEnvironmentInfo()
    # print ( f"This machine has {nCPUs()} CPUs" )
