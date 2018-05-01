#!/usr/bin/env python

"""
.. module:: runtime
    :synopsis: Tools to gather info about runtime enviroment,
               ( nCPUs() ), or obtain file type ( filetype() ). Pointer
               to the Flong calculator is also kept here. 

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from smodels.installation import installDirectory
from smodels.tools.smodelsLogging import logger
import imp,os


def filetype ( filename ):
    """ obtain information about the filetype of an input file,
        currently only used to discriminate between slha and lhe
        files.

        :returns: filetype as string("slha" or "lhe"),
                  None if file does not exist, or filetype is unknown.
    """
    import os
    if not os.path.exists ( filename ):
        return None
    if filename[-5:].lower() == ".slha":
        return "slha"
    if filename[-4:].lower() == ".lhe":
        return "lhe"
    with open ( filename ) as f:
        for line in f:
            if "<LesHouchesEvents" in line:
                return "lhe"
            if "<event>" in line:
                return "lhe"
            if "block " in line.lower():
                return "slha"
    return None

def nCPUs():
    """ obtain the number of CPU cores on the machine, for several
        platforms and python versions. """
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except ImportError as e:
        pass
    try:
        import psutil
        return psutil.NUM_CPUS
    except ImportError as e:
        pass
    try:
        import os
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))
        if res>0: return res
    except ImportError as e:
        pass
    return None

def setFlongCalc(path=os.path.abspath("%ssmodels/tools/flongCalc.py" %installDirectory()),
             method="FlongCalculator"):
    """
    Define/replace the function used by the SLHA decomposer to compute
    the fraction of long-lived and prompt decays.
    It can be used to apply (external) user-defined function for
    the calculation.  
 
    :param method: Name of the function.
    :param path: Path to the .py file defining the function.
    """
     
    from smodels.theory import slhaDecomposer
     
    try:
        mod = imp.load_source(os.path.basename(path).replace('.py',''),path)        
        fCalc = getattr(mod,method)
        logger.debug("Using Flong calculator: %s from %s" %(method,path))
    except:
        from smodels.tools.flongCalc import FlongCalculator as fCalc 
        logger.warning("Could not load Flong calculator %s from %s. Using default." %(method,path))
     
    slhaDecomposer.fCalc = fCalc
    
    
    

if __name__ == "__main__":
    print ( "This machine has %d CPUs" % nCPUs() )