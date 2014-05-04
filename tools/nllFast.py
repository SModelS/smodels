#!/usr/bin/env python

"""
.. module:: nllFast
   :synopsis: This module provides methods to access the nllfast grid and
   compute k-factors (when available) to SUSY pair production cross-sections.

.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import commands
import os
from . import setPath
from . import toolBox
from . import modpyslha as pyslha
from .physicsUnits import rmvunit
from .physicsUnits import addunit
import numpy
import operator

import SModelS
import logging

logger = logging.getLogger(__name__)

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

    # Get process name (in NLLfast notation)
    process = getProcessName(pIDs)
    if not process:
        # Return k-factors = None, if NLLfast does not have the process
        return (None, None)

    # Obtain relevant masses
    readfile = pyslha.readSLHAFile(slhafile)
    gluinomass = abs(readfile.blocks['MASS'].entries[1000021])
    squarkmass = sum([abs(readfile.blocks['MASS'].entries[pid])
                      for pid in squarks]) / 8.
    pid1, pid2 = sorted(pIDs)
    if pid1 in antisquarks and pid2 in squarks:
        squarkmass = (abs(readfile.blocks['MASS'].entries[abs(pid1)]) +
                      abs(readfile.blocks['MASS'].entries[pid2])) / 2
    elif pid1 in squarks and pid2 in squarks:
        squarkmass = (abs(readfile.blocks['MASS'].entries[pid1]) +
                      abs(readfile.blocks['MASS'].entries[pid2])) / 2
    elif abs(pid1) == pid2 and pid2 in third:
        squarkmass = abs(readfile.blocks['MASS'].entries[abs(pid1)])

    # Set up NLLfast run, the old way
    sqrtS = float(rmvunit(sqrts, 'TeV'))
    energy = str(int(sqrtS)) + 'TeV'
    toolname = "nllfast%d" % int(sqrtS)
    box = toolBox.ToolBox()
    tool = box.get(toolname)
    if tool == None:
        logger.warning("No NLLfast data for sqrts = " + str(sqrts))
        return (None, None)
    nllpath = tool.installDirectory()
    tool.pathOfExecutable()
    tool.checkInstallation()
    if process == "st":
        nll_run = "./nllfast_" + energy + " %s %s %s" % \
                  (process, pdf, squarkmass)
    else:
        nll_run = "./nllfast_" + energy + " %s %s %s %s" % \
                  (process, pdf, squarkmass, gluinomass)

    # Run NLLfast
    nll_output = runNLLfast(nll_run, nllpath)

    # If run was successful, return k-factors:
    if "K_NLO" in nll_output:
        # NLLfast ran ok, try to get the k-factors
        kFacs = getKfactorsFrom(nll_output)
        if not kFacs or min(kFacs) <= 0.:
            logger.warning("Error obtaining k-factors")
            return (None, None)
        else:
            return kFacs
    # If run was not successful, check for decoupling error messages:
    elif not "too low/high" in nll_output.lower():
        logger.warning("Error running NLLfast")
        return (None, None)

    # Deal with decoupling regimes
    gluino_dcp, squark_dcp = False, False
    if "too low/high gluino" in nll_output.lower():
        gluino_dcp = True
    elif "too low/high squark" in nll_output.lower():
        squark_dcp = True

    # If produced particles are too heavy, return k-factors = 1
    if gluino_dcp and ('g' in process or gluinomass < 500.):
        logger.warning("Gluino mass out of NLLfast grid for " + process)
        return (None, None)
    elif squark_dcp and (process != 'gg' or squarkmass < 500.):
        logger.warning("Squark mass out of NLLfast grid for " + process)
        return (None, None)

    # If virtual particles are too heavy, interpolate to the decoupling limit
    kFacsVector = []
    xmass = max(squarkmass, gluinomass)  # To be interpolated
    while len(kFacsVector) < 2 and xmass > 500.:
        xmass -= 100.  # Reduce decoupled mass, until NLLfast produces results
        if squark_dcp:
            nll_run = "./nllfast_" + energy + " %s %s %s %s" % \
                      (process, pdf, xmass, gluinomass)
        elif gluino_dcp:
            nll_run = "./nllfast_" + energy + " %s %s %s %s" % \
                      (process, pdf, squarkmass, xmass)
        nll_output = runNLLfast(nll_run, nllpath)
        if "K_NLO" in nll_output:
            kfacs = getKfactorsFrom(nll_output)
            if kfacs and min(kfacs) > 0.:
                kFacsVector.append([xmass, kfacs])

    if len(kFacsVector) < 2:
        logger.warning("Not enough points for interpolation in the decoupling "
                       "limit")
        return (None, None)
    else:
        if kFacsVector[0][0] / min(squarkmass, gluinomass) > 10.:
            # Decoupling limit is satisfied, do not interpolate
            kFacs = kFacsVector[0][1]
        else:
            # Interpolate k-factors
            kFacs = interpolateKfactors(kFacsVector,
                                        max(squarkmass, gluinomass))
    return kFacs



def getProcessName(pIDs):
    """
    Return the process name (in NLLfast notation) for the pair production of
    pIDs.
    
    :returns: None, if the particle ID pair is not contained in NLLfast
    
    """
    pid1, pid2 = sorted(pIDs)
    process = None

    # Obtain the type of process:
    #  - gluino-gluino production = gg
    #  - squark-antisquark = sb
    #  - squark-squark = ss
    #  - squark-gluino = sg
    #  - antistop-stop
    #  - antisbottom-sbottom = st
    if pid1 in antisquarks and pid2 in squarks:
        process = 'sb'
    elif pid1 in squarks and pid2 in squarks:
        process = 'ss'
    elif pid1 == pid2 and pid1 in gluinos:
        process = 'gg'
    elif pid1 in squarks and pid2 == 1000021 or \
            pid2 in squarks and pid1 == 1000021:
        process = 'sg'
    elif abs(pid1) == pid2 and pid2 in third:
        process = 'st'

    return process


def runNLLfast(nll_run, nllpath):
    """
    Execute NLLfast with command nll_run at nllpath.
    
    :returns: NLLfast output as a string
    
    """
    current_dir = os.getcwd()
    os.chdir(nllpath)
    nll_output = commands.getoutput(nll_run)
    os.chdir(current_dir)
    return nll_output


def getKfactorsFrom(output):
    """
    Read NLLfast output and return the k-factors.
    
    """
    if not output:
        return False
    else:
        lines = output.split('\n')
        il = 0
        line = lines[il]
        while not "K_NLO" in line and il < len(lines) - 2:
            if "process" in line:
                process = line[line.find("process:") + 8:].replace(" ", "")
            il += 1
            line = lines[il]
        if not process:
            return False
        # Line with header
        line = lines[il]
        # Count number of mass entries
        nmass = line.count('GeV')
        # Line with values
        line = lines[il + 2]
        data = [eval(x) for x in line.split()]
        if len(data) != nmass + 11:
            return False
        else:
            kFacs = data[9 + nmass:]

    return kFacs


def interpolateKfactors(kFacsVector, xval):
    """
    Interpolate a list of k-factor  values from
    kFacsVector = [[x0,[k1,k2,..]], [x1,[k1,k2,..],...].
    
    :returns: list of interpolated k-factor values at x-value xval
    
    """
    kFacs = []

    xpts = [x[0] for x in kFacsVector]
    ypts = [x[1] for x in kFacsVector]
    coeffs = numpy.matrix.transpose(numpy.polyfit(xpts, ypts, len(xpts) - 1))
    for ik in range(len(ypts[0])):
        kfac = 0.
        for ip, coeff in enumerate(coeffs[ik]):
            kfac += coeff * xval ** (len(xpts) - 1 - ip)
        if kfac <= 0.:
            kfac = 1.
        kFacs.append(kfac)

    return kFacs


if __name__ == "__main__":
    """
    Calculate k factors for a pid pair.
    
    """
    slhaF = SModelS.installDirectory() + "inputFiles/slha/T1.slha"
    kNLO, kNLL = getKfactorsFor((1000021, 1000021), addunit(13., "TeV"), slhaF)
    print("nlo, nll = " + str(kNLO) + ", " + str(kNLL))
