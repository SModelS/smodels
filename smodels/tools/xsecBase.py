#!/usr/bin/env python3

"""
.. module:: xsecBasis
   :synopsis: Computation of reference ("theory") production cross sections.

.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Théo Reymermier <theo.reymermier@gmail.com>

"""

from __future__ import print_function
import sys
import os, copy
current = os.getcwd()
sys.path.append(current)


from smodels import installation
from smodels.tools import toolBox
from smodels.base import runtime
from smodels.base.physicsUnits import pb, TeV, GeV
from smodels.base import crossSection
from smodels.base.crossSection import LO, NLO, NLL
from smodels.base.smodelsLogging import logger, setLogLevel
from smodels.decomposition.exceptions import SModelSDecompositionError as SModelSError
import subprocess
from concurrent.futures import ProcessPoolExecutor

import pyslha
import math
try:
    import cStringIO as io
except ImportError as e:
    import io



class XSecBase:
    """ cross section computer class, what else? """
    def __init__ ( self, maxOrder,slha_folder_name, maycompile=True):
        """
        :param maxOrder: maximum order to compute the cross section, given as an integer
                    if maxOrder == LO, compute only LO pythia xsecs
                    if maxOrder == NLO, apply NLO K-factors from NLLfast (if available)
                    if maxOrder == NLL, apply NLO+NLL K-factors from NLLfast (if available)
        :param nevents: number of events for pythia run
        :param pythiaVersion: pythia6 or pythia8 (integer)
        :param maycompile: if True, then tools can get compiled on-the-fly
        """
        self.resummino_bin = "./smodels/lib/resummino/resummino-3.1.2/bin/resummino"
        self.input_file_original = "smodels/etc/ff1a240db6c1719fe9f299b3390d49d32050c4f1003286d2428411eca45bd50c.in"
        self.slha_folder_name = slha_folder_name
        self.maxOrder = maxOrder
        self.countNoXSecs = 0
        self.countNoNLOXSecs = 0
        self.maycompile = maycompile

    def addXSecToFile( self, xsecs, slhafile, comment=None, complain=True):
        """
        Write cross sections to an SLHA file.

        :param xsecs: a XSectionList object containing the cross sections
        :param slhafile: target file for writing the cross sections in SLHA format
        :param comment: optional comment to be added to each cross section block
        :param complain: complain if there are already cross sections in file

        """

        if not os.path.isfile(slhafile):
            line = f"SLHA file {slhafile} not found."
            logger.error( line )
            raise SModelSError( line )
        if len(xsecs) == 0:
            self.countNoXSecs+=1
            if self.countNoXSecs < 3:
                logger.warning("No cross sections available.")
            if self.countNoXSecs == 3:
                logger.warning("No cross sections available (will quench such warnings in future).")
            return False
        # Check if file already contain cross section blocks
        xSectionList = crossSection.getXsecFromSLHAFile(slhafile)
        if xSectionList and complain:
            logger.info("SLHA file already contains XSECTION blocks. Adding "
                           "only missing cross sections.")

        # Write cross sections to file, if they do not overlap with any cross
        # section in the file
        outfile = open(slhafile, 'r')
        lastline = outfile.readlines()[-1]
        lastline = lastline.strip()
        outfile.close()

        outfile = open(slhafile, 'a')
        if lastline != "":
            outfile.write("\n" )
        nxsecs = 0
        for xsec in xsecs:
            writeXsec = True
            for oldxsec in xSectionList:
                if oldxsec.info == xsec.info and set(oldxsec.pid) == set(xsec.pid):
                    writeXsec = False
                    break
            if writeXsec:
                nxsecs += 1
                outfile.write( self.xsecToBlock(xsec, (2212, 2212), comment) + "\n")
        outfile.close()

        return nxsecs 

    def xsecToBlock( self, xsec, inPDGs=(2212, 2212), comment=None, xsecUnit = pb):
        """
        Generate a string for a XSECTION block in the SLHA format from a XSection
        object.

        :param inPDGs: defines the PDGs of the incoming states
                    (default = 2212,2212)

        :param comment: is added at the end of the header as a comment
        :param xsecUnit: unit of cross sections to be written (default is pb).
                        Must be a Unum unit.

        """
        if type(xsec) != type(crossSection.XSection()):
            logger.error("Wrong input")
            raise SModelSError()
        # Sqrt(s) in GeV
        header = "XSECTION  " + str(xsec.info.sqrts / GeV)
        for pdg in inPDGs:
            # PDGs of incoming states
            header += " " + str(pdg)
        # Number of outgoing states
        header += " " + str(len(xsec.pid))
        for pid in xsec.pid:
            # PDGs of outgoing states
            header += " " + str(pid)
        if comment:
            header += " # " + str(comment)  # Comment
        entry = "  0  " + str(xsec.info.order) + "  0  0  0  0  " + \
                str( f"{float(xsec.value / xsecUnit):16.8E}" ) + " SModelSv" + \
                    installation.version()

        return "\n" + header + "\n" + entry

class ArgsStandardizer:
    """ simple class to collect all argument manipulators """

    def getSSMultipliers ( self, multipliers ):
        if type ( multipliers ) == str:
            if multipliers in [ "", "None", "none", "no", "{}" ]:
                return None
            if multipliers.count("{") != 1 or multipliers.count("}") != 1:
                logger.error ( "need to pass signal strengh multipliers as dictionary with tuple of pids as keys" )
            if multipliers.count("(") != multipliers.count(")"):
                logger.error ( "need to pass signal strengh multipliers as dictionary with tuple of pids as keys" )
            return eval(multipliers)
        return multipliers

    def getInputFiles ( self, args ):
        """ geth the names of the slha files to run over """
        inputPath  = args.filename.strip()
        if not os.path.exists( inputPath ):
            logger.error( f"Path {inputPath} does not exist." )
            sys.exit(1)
        inputFiles = []
        if os.path.isfile ( inputPath ):
            inputFiles = [ inputPath ]
        else:
            files = os.listdir ( inputPath )
            for f in files:
                inputFiles.append ( os.path.join ( inputPath, f ) )
        import random
        random.shuffle ( inputFiles )
        return inputFiles

    def checkAllowedSqrtses ( self, order, sqrtses ):
        """ check if the sqrtses are 'allowed' """
        if order == 0: return
        allowedsqrtses=[7, 8, 13, 13.6]
        for sqrts in sqrtses:
            if not sqrts in allowedsqrtses:
                logger.error("Cannot compute NLO or NLL xsecs for sqrts = %d "
                        "TeV! Available are: %s TeV." % (sqrts, allowedsqrtses ))
                sys.exit(-2)

    def getOrder ( self, args ):
        """ retrieve the order in perturbation theory from argument list """
        if args.NLL:
            return NLL
        if args.NLO:
            return NLO
        return LO

    def getjson ( self, args ):
        """ retrieve the path to the json file from argument list """
        json = args.conf

        if json == 'default':
            return None
        else:
            if os.path.exists(json):
                return json
            else:
                return logger.error("Path does not exist.")

    def queryCrossSections ( self, filename ):
        if os.path.isdir ( filename ):
            logger.error ( "Cannot query cross sections for a directory." )
            sys.exit(-1)
        xsecsInfile = crossSection.getXsecFromSLHAFile(filename)
        if xsecsInfile:
            print ( "1" )
        else:
            print ( "0" )

    def getSqrtses ( self, args ):
        """ extract the sqrtses from argument list """
        if args.sqrts is None or len(args.sqrts) == 0:
            return {13}
        sqrtses = [item for sublist in args.sqrts for item in sublist]
        sqrtses.sort()
        sqrtses = set(sqrtses)
        return sqrtses

    def getParticles ( self, args ):
        """ extract the particles from argument list, default to None, then channels are chosen by the json file """
        if not hasattr ( args, "particles" ) or len(args.particles) == 0:
            return None
        particles = [int(item) for sublist in args.particles for item in sublist]
        if len(particles) == 0:
            particles= None
            return particles
        particles.sort()
        particles = set(particles)
        return particles
    
    def checkNCPUs ( self, ncpus, inputFiles ):
        if ncpus < -1 or ncpus == 0:
            logger.error ( "Weird number of CPUs given: %d" % ncpus )
            sys.exit()
        if ncpus == -1:
            ncpus = runtime.nCPUs()
        ncpus = min ( len(inputFiles), ncpus )
        if ncpus == 1:
            logger.info ( "We run on a single cpu" )
        else:
            logger.info ( "We run on %d cpus" % ncpus )
        return ncpus

    def tempDir ( self, args ):
        ret = "/tmp/"
        if hasattr ( args, "tempdir" ):
            ret = args.tempdir
        return ret


    def checkXsec_limit (self,args ):
        if args.xseclimit == None:
            return None
        if args.xseclimit < 0:
            logger.error("Xsec_limit cannot be negative. Set to 0 if you want no limitation.")
            sys.exit()
        if args.xseclimit > 1000:
            logger.warn("Are you sure to use a limit that high ? you might get errors.")
        return args.xseclimit

    def getPythiaVersion ( self, args ):
        pythiaVersion = 8

        if hasattr(args, 'pythia6' ) and args.pythia6 == True:
            pythiaVersion = 6
            if hasattr(args, 'pythia8') and args.pythia8 == True:
                logger.error ( "cannot both use pythia6 and pythia8 for LO xsecs." )
                sys.exit()
        return pythiaVersion

    def writeToFile ( self, args ):
        toFile = None
        if args.tofile:
            toFile="highest"
        if args.alltofile:
            if toFile=="highest":
                logger.warning ( "Specified both --tofile and --alltofile. Will use "\
                              "--alltofile" )
            toFile="all"
        return toFile
