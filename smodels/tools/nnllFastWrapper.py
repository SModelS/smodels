#!/usr/bin/env python3

"""
.. module:: nnllFastWrapper
   :synopsis: This module provides methods to access the nnllfast grid and
              compute k-factors (when available) to SUSY pair
              production cross sections.

.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
from __future__ import print_function
import operator
import pyslha

try:
    import commands as executor
except ImportError:
    import subprocess as executor

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

import os
from smodels.tools.wrapperBase import WrapperBase
from smodels.tools.smodelsLogging import logger


class NnllfastWrapperWrapper(WrapperBase):
    """
    An instance of this class represents the installation of nnllfast.

    """
    def __init__(self, sqrts, nnllfastVersion, testParams, testCondition):
        """
        :param sqrts: sqrt of s, in TeV, as an integer,
        :param nnllfastVersion: version of the nnllfast tool
        :param testParams: what are the test params we need to run things with?
        :param testCondition: the line that should be the last output line when
        running executable
        :srcPath: the path of the source code, for compilation

        """
        WrapperBase.__init__(self)
        self.sqrts = int(sqrts)
        self.name = "nnllfast%d" % sqrts
        self.nnllfastVersion = nnllfastVersion
        path = "<install>/smodels/lib/nnllfast/nnllfast-"
        location = path + self.nnllfastVersion + "/"
        self.cdPath = self.absPath(location)
        self.executablePath = self.cdPath + "/nnllfast_%dTeV" % self.sqrts
        self.testParams = testParams
        self.testCondition = testCondition
        self.srcPath = self.cdPath
        self.compiler = "gfortran"
        self.executable = ""

    def _interpolateKfactors( self, kFacsVector, xval):
        """
        Interpolate a list of k-factor  values from
        kFacsVector = [[x0,[k1,k2,..]], [x1,[k1,k2,..],...].
        FIXME what is xval?

        :returns: list of interpolated k-factor values at x-value xval

        """
        import numpy
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

    def _getKfactorsFrom( self, output ):
        """
        Read NNLLfast output and return the k-factors.

        """
        if not output:
            return False
        else:
            lines = output.split('\n')
            il = 0
            line = lines[il]
            process = False
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

    def _run ( self, this ):
        """
        Run. Code taken from nnllFast.runNNLLfast
        Return the process name (in NNLLfast notation) for the pair production of
        pIDs.

        :returns: None, if the particle ID pair is not contained in NNLLfast

        """
        current_dir = os.getcwd()
        os.chdir ( self.cdPath )
        nnll_output = executor.getoutput( this )
        os.chdir(current_dir)
        return nnll_output

    def _compute ( self, energy, pIDs, squarkmass, gluinomass ):
        process = self._getProcessName(pIDs)
        if process == "st":
            nnll_run = "./nnllfast_" + energy + " %s %s" % \
                  (process, squarkmass)
        else:
            nnll_run = "./nnllfast_" + energy + " %s %s %s" % \
                  (process, squarkmass, gluinomass)
        return self._run( nnll_run )

    def _getProcessName(self, pIDs):
        """
        Return the process name (in NNLLfast notation) for the pair production of
        pIDs.

        :returns: None, if the particle ID pair is not contained in NNLLfast
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
        elif abs(pid1) in squarks and abs(pid2) in squarks:
            process = 'ss'
        elif pid1 == pid2 and pid1 in gluinos:
            process = 'gg'
        elif abs(pid1) in squarks and pid2 == 1000021 or \
                abs(pid2) in squarks and pid1 == 1000021:
            process = 'sg'
        elif abs(pid1) == pid2 and pid2 in third:
            process = 'st'
        return process

    def _getDecoupledKfactors( self, process, energy, mass ):
        """
        Compute k-factors in the decoupled (squark or gluino) regime for the process.
        If a decoupled grid does not exist for the process, return None
        """

        if process != 'sb' and process != 'gg': return None
        elif process == 'sb': process_dcpl = 'sdcpl'
        elif process == 'gg': process_dcpl = 'gdcpl'
        nnll_run = "./nnllfast_" + energy + " %s %s" % \
                          (process_dcpl, mass)
        e = energy.replace ( "TeV", "" ).replace ( "*", "" )
        # tool = toolBox.ToolBox().get ( "nnllfast%d" % int ( e ) )
        print ( "run nnllfast", nnll_run )
        nnll_output = self._run( nnll_run )
        if "K_NLO" in nnll_output:
            return self._getKfactorsFrom(nnll_output)
        else: return None

    def _runForDecoupled ( self, energy, nnllinput ):
        nnll_run = "./nnllfast_" + energy + " %s %s %s" % nnllinput
        return self._run ( nnll_run )

    def getKfactorsFor( self, pIDs, slhafile ):
        """
        Read the NNLLfast grid and returns a pair of k-factors (NNLL) for
        the PIDs pair.

        :returns: k-factors = None, if NNLLfast does not contain the process; uses
                  the slhafile to obtain the SUSY spectrum.

        """
        if not os.path.isfile(slhafile):
            logger.error("SLHA file %s not found", slhafile)
            return False

        energy = str(int(self.sqrts)) + 'TeV'
        # Get process name (in NNLLfast notation)
        process = self._getProcessName(pIDs)
        if not process:
            # Return k-factors = None, if NNLLfast does not have the process
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

        #if tool == None:
        #    logger.warning("No NNLLfast data for sqrts = " + str(sqrts))
        #    return (None, None)
        nnllpath = self.installDirectory()
        # self.pathOfExecutable()
        self.checkInstallation()
        nnll_output = self._compute ( energy, pIDs, squarkmass, gluinomass )

        # If run was successful, return k-factors:
        if "K_NLO" in nnll_output:
            # NNLLfast ran ok, try to get the k-factors
            kFacs = self._getKfactorsFrom(nnll_output)
            if not kFacs or min(kFacs) <= 0.:
                logger.warning("Error obtaining k-factors")
                return (None, None)
            else:
                return kFacs
        # If run was not successful, check for decoupling error messages:
        elif not "too low/high" in nnll_output.lower():
            logger.warning("Error running NNLLfast")
            return (None, None)

        # Check for decoupling cases with a decoupling grid (only for sb and gg)
        doDecoupling = False
        if "too low/high gluino" in nnll_output.lower():
            if gluinomass > 500. and process == 'sb':
                doDecoupling = True
                dcpl_mass = gluinomass
        elif "too low/high squark" in nnll_output.lower():
            if squarkmass > 500. and process == 'gg':
                doDecoupling = True
                dcpl_mass = squarkmass

        # If process do not have decoupled grids, return None:
        if not doDecoupling:
            logger.warning("Masses of (q,g)=(%s,%s) out of NNLLfast grid for %s, %s" % ( squarkmass, gluinomass, process, energy ))
            return (None, None)

        # Obtain k-factors from the NNLLfast decoupled grid
        kfacs = self._getDecoupledKfactors(process,energy,min(gluinomass,squarkmass))
        # Decoupling limit is satisfied, do not interpolate
        if not kfacs:
            logger.warning("Error obtaining k-factors from the NNLLfast decoupled grid for " + process)
            return (None, None)
        elif dcpl_mass/min(gluinomass,squarkmass) > 10.:
            return kfacs
        # Interpolate between the non-decoupled and decoupled grids
        else:
            kFacsVector = [[10.*min(gluinomass,squarkmass),kfacs]]  #First point for interpolation (decoupled grid)
            kfacs = None
            while not kfacs and dcpl_mass > 500.:
                dcpl_mass -= 100.  # Reduce decoupled mass, until NNLLfast produces results
                if process == 'sb': nnllinput = (process, squarkmass, dcpl_mass)
                else:  nnllinput = (process, dcpl_mass, gluinomass)
                nnll_output = self._runForDecoupled ( energy, nnllinput )
                kfacs = self._getKfactorsFrom(nnll_output)
            kFacsVector.append([dcpl_mass, kfacs]) #Second point for interpolation (non-decoupled grid)

        if len(kFacsVector) < 2:
            logger.warning("Not enough points for interpolation in the decoupling "
                           "limit")
            return (None, None)
        else:
            # Interpolate k-factors
            kFacs = self._interpolateKfactors(kFacsVector,
                            max(squarkmass, gluinomass))
        return kFacs

class NnllfastWrapperWrapper13(NnllfastWrapperWrapper):
    """
    An instance of this class represents the installation of nnllfast 13.

    """
    def __init__(self):
        NnllfastWrapperWrapper.__init__(self, 13, "1.1",
                                 testParams="gg 500 600",
                                 testCondition="600.    0.394E+01  0.690E+01  "
                                 "0.731E+01    0.394E+00" )

nnllFastTools = { 13 : NnllfastWrapperWrapper13() }


if __name__ == "__main__":
    for (sqrts, tool) in nnllFastTools.items():
        print("%s: installed in %s" % (tool.name, tool.installDirectory()))
