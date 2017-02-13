#!/usr/bin/env python

"""
.. module:: nllFast
   :synopsis: This module provides methods to access the nllfast grid and
              compute k-factors (when available) to SUSY pair 
              production cross sections.

.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
#
try:
    import commands as executor
except ImportError:
    import subprocess as executor
import os
from smodels.tools import toolBox
import pyslha
from smodels.tools.physicsUnits import TeV
import numpy
import operator

from smodels import installation
from smodels.tools.smodelsLogging import logger

squarks = [1000001,
           2000001,
           1000002,
           2000002,
           1000003,
           2000003,
           1000004,
           2000004]
antisquarks = map(operator.neg, squarks)
third = [1000005,
         2000005,
         1000006,
         2000006]
gluinos = [1000021]

def getKfactorsFor(pIDs, sqrts, slhafile, pdf='cteq'):
    """
    Read the NLLfast grid and returns a pair of k-factors (NLO and NLL) for the
    pair.

    :returns: k-factors = None, if NLLfast does not contain the process; uses
              the slhafile to obtain the SUSY spectrum.
    
    """
    if not os.path.isfile(slhafile):
        logger.error("SLHA file %s not found", slhafile)
        return False

    box = toolBox.ToolBox()
    # Set up NLLfast run, the old way
    sqrtS = float(sqrts/TeV)
    energy = str(int(sqrtS)) + 'TeV'
    toolname = "nllfast%d" % int(sqrtS)
    tool = box.get(toolname)
    # Get process name (in NLLfast notation)
    process = tool.getProcessName(pIDs)
    if not process:
        # Return k-factors = None, if NLLfast does not have the process
        return (None, None)

    # Obtain relevant masses
    readfile = pyslha.readSLHAFile(slhafile)
    masses=readfile.blocks['MASS']
    check_pids=squarks+gluinos+third
    for check in check_pids:
        if not check in masses.entries:
            logger.error ( "cannot compute k factor for pdgid %d: " \
              " no particle mass given. will set mass to inf." % check )
            masses.entries[check]=1.e10

    gluinomass = abs(masses.entries[1000021])
    squarkmass = sum([abs(masses.entries[pid])
                      for pid in squarks]) / 8.
    pid1, pid2 = sorted(pIDs)
    if pid1 in antisquarks and pid2 in squarks:
        squarkmass = (abs(masses.entries[abs(pid1)]) +
                      abs(masses.entries[pid2])) / 2.
    elif pid1 in squarks and pid2 in squarks:
        squarkmass = (abs(masses.entries[pid1]) + abs(masses.entries[pid2])) / 2.
    elif abs(pid1) == pid2 and pid2 in third:
        squarkmass = abs(masses.entries[abs(pid1)])

    if tool == None:
        logger.warning("No NLLfast data for sqrts = " + str(sqrts))
        return (None, None)
    nllpath = tool.installDirectory()
    tool.pathOfExecutable()
    tool.checkInstallation()
    nll_output = tool.compute ( energy, pIDs, pdf, squarkmass, gluinomass )

    # If run was successful, return k-factors:
    if "K_NLO" in nll_output:
        # NLLfast ran ok, try to get the k-factors
        kFacs = tool.getKfactorsFrom(nll_output)
        if not kFacs or min(kFacs) <= 0.:
            logger.warning("Error obtaining k-factors")
            return (None, None)
        else:
            return kFacs
    # If run was not successful, check for decoupling error messages:
    elif not "too low/high" in nll_output.lower():
        logger.warning("Error running NLLfast")
        return (None, None)

    # Check for decoupling cases with a decoupling grid (only for sb and gg)
    doDecoupling = False
    if "too low/high gluino" in nll_output.lower():        
        if gluinomass > 500. and process == 'sb': 
            doDecoupling = True
            dcpl_mass = gluinomass
    elif "too low/high squark" in nll_output.lower():
        if squarkmass > 500. and process == 'gg':
            doDecoupling = True
            dcpl_mass = squarkmass

    # If process do not have decoupled grids, return None:
    if not doDecoupling:
        logger.warning("Masses out of NLLfast grid for " + process)
        return (None, None)

    # Obtain k-factors from the NLLfast decoupled grid
    kfacs = tool.getDecoupledKfactors(process,energy,pdf,min(gluinomass,squarkmass))
    # Decoupling limit is satisfied, do not interpolate
    if not kfacs:
        logger.warning("Error obtaining k-factors from the NLLfast decoupled grid for " + process)
        return (None, None)
    elif dcpl_mass/min(gluinomass,squarkmass) > 10.:    
        return kfacs
    # Interpolate between the non-decoupled and decoupled grids
    else:
        kFacsVector = [[10.*min(gluinomass,squarkmass),kfacs]]  #First point for interpolation (decoupled grid)
        kfacs = None        
        while not kfacs and dcpl_mass > 500.:
            dcpl_mass -= 100.  # Reduce decoupled mass, until NLLfast produces results
            if process == 'sb': nllinput = (process, pdf, squarkmass, dcpl_mass)
            else:  nllinput = (process, pdf, dcpl_mass, gluinomass)
            nll_output = tool.runForDecoupled ( energy, nllinput )
            #nll_run = "./nllfast_" + energy + " %s %s %s %s" % nllinput
            #nll_output = tool.run(nll_run )
            kfacs = tool.getKfactorsFrom(nll_output)        
        kFacsVector.append([dcpl_mass, kfacs]) #Second point for interpolation (non-decoupled grid)

    if len(kFacsVector) < 2:
        logger.warning("Not enough points for interpolation in the decoupling "
                       "limit")
        return (None, None)
    else:
        # Interpolate k-factors
        kFacs = tool.interpolateKfactors(kFacsVector,
                        max(squarkmass, gluinomass))
    return kFacs

if __name__ == "__main__":
    """
    Calculate k factors for a pid pair.
    """
    slhaF = installation.installDirectory() + "inputFiles/slha/T1.slha"
    # tool = toolBox.ToolBox().get ( "nllfast13" )
    kNLO, kNLL = getKfactorsFor((1000021, 1000021), 13.*TeV, slhaF)
    print("nlo, nll = " + str(kNLO) + ", " + str(kNLL))
