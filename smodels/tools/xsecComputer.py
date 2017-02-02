#!/usr/bin/env python

"""
.. module:: xsecComputer
   :synopsis: Computation of reference ("theory") production cross sections.

.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
from __future__ import print_function
from smodels import installation
from smodels.tools import toolBox, runtime
from smodels.tools.physicsUnits import pb, TeV, GeV
from smodels.theory import crossSection
from smodels.theory.crossSection import LO, NLO, NLL
from smodels.tools import nllFast
from smodels.tools.smodelsLogging import logger, setLogLevel
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
import os, copy
import pyslha
try:
    import cStringIO as io
except ImportError as e:
    import io
import sys

class XSecComputer:
    """ cross section computer class, what else? """
    def __init__ ( self, maxOrder, nevents, pythiaVersion ):
        """
        :param maxOrder: maximum order to compute the cross section, given as an integer
                    if maxOrder == LO, compute only LO pythia xsecs
                    if maxOrder == NLO, apply NLO K-factors from NLLfast (if available)
                    if maxOrder == NLL, apply NLO+NLL K-factors from NLLfast (if available)
        :param nevents: number of events for pythia run
        :param pythiaVersion: pythia6 or pythia8 (integer)
        """
        self.maxOrder = self._checkMaxOrder ( maxOrder )
        if nevents < 1:
            logger.error ( "Supplied nevents < 1" )
            sys.exit()
        self.nevents = nevents
        if pythiaVersion not in [ 6, 8 ]:
            logger.error ( "Unknown pythia version %s. Allowed values: 6, 8" % \
                            ( pythiaVersion ) )
            sys.exit()
        self.pythiaVersion = pythiaVersion 

    def _checkSLHA ( self, slhafile ):
        if not os.path.isfile(slhafile):
            logger.error("SLHA file %s not found.", slhafile)
            raise SModelSError()
        try:
            f=pyslha.readSLHAFile(slhafile)
        except pyslha.ParseError as e:
            logger.error("File cannot be parsed as SLHA file: %s" % e )
            raise SModelSError()

    def _checkSqrts ( self, sqrts ):
        if type(sqrts)==type(float) or type(sqrts)==type(int):
            logger.warning("sqrt(s) given as scalar, will add TeV as unit." )
            sqrts=float(sqrts)*TeV
        return sqrts

    def _checkMaxOrder ( self, maxOrder ):
        smaxorder={ "LO": 0, "NLO": 1, "NLL": 2 }
        if maxOrder in smaxorder.keys():
            logger.warning("maxorder given as string, please supply integer.")
            maxOrder=smaxorder[maxOrder]
        return maxOrder

    def addHigherOrders ( self, sqrts, slhafile ):
        """ add higher order xsecs """
        xsecs = copy.deepcopy ( self.loXsecs )
        wlabel = str(int(sqrts / TeV)) + ' TeV'
        if self.maxOrder == LO:
            wlabel += ' (LO)'
        elif self.maxOrder == NLO:
            wlabel += ' (NLO)'
        elif self.maxOrder == NLL:
            wlabel += ' (NLO+NLL)'
        for ixsec, xsec in enumerate(xsecs):
            xsecs[ixsec].info.label = wlabel
            xsecs[ixsec].info.order = self.maxOrder

        if self.maxOrder > 0:
            pIDs = self.loXsecs.getPIDpairs()
            for pID in pIDs:
                k = 0.
                kNLO, kNLL = nllFast.getKfactorsFor(pID, sqrts, slhafile)
                if self.maxOrder == NLO and kNLO:
                    k = kNLO
                elif self.maxOrder == NLL and kNLL and kNLO:
                    k = kNLO * kNLL
                elif self.maxOrder > 2 and kNLL and kNLO:
                    logger.warning("Unkown xsec order, using NLL+NLO k-factor, "
                                   "if available")
                    k = kNLO * kNLL
                k = float(k)
                for i, xsec in enumerate(xsecs):
                    if set(xsec.pid) == set(pID):
                        # Apply k-factor
                        xsecs[i] = xsec * k

        # Remove zero cross sections
        while len(xsecs) > 0 and xsecs.getMinXsec() == 0. * pb:
            for xsec in xsecs:
                if xsec.value == 0. * pb:
                    xsecs.delete(xsec)
                    break
        if self.maxOrder > 0 and len(xsecs) == 0:
            logger.warning("No NLO or NLL cross sections available.")

        #for i in xsecs:
        #    logger.error ( "xsec=%s (%s)" % (i,type(i)) )
        return xsecs

    def compute ( self, sqrts, slhafile,  lhefile=None, unlink=True, loFromSlha=None, 
                  pythiacard=None ):
        """
        Run pythia and compute SUSY cross sections for the input SLHA file.

        :param sqrts: sqrt{s} to run Pythia, given as a unum (e.g. 7.*TeV)
        :param slhafile: SLHA file
        :param lhefile: LHE file. If None, do not write pythia output to file. If
                        file does not exist, write pythia output to this file name. If
                        file exists, read LO xsecs from this file (does not run pythia).
        :param unlink: Clean up temp directory after running pythia

        :param loFromSlha: If True, uses the LO xsecs from the SLHA file to compute the
                           higher order xsecs
        :param pythiaCard: Optional path to pythia.card. If None, uses /etc/pythia.card

        :returns: XSectionList object

        """
        sqrts = self._checkSqrts( sqrts )
        self._checkSLHA ( slhafile )

        if lhefile:
            if os.path.isfile(lhefile):
                logger.warning("Using LO cross sections from " + lhefile)
                logger.error ( "Cross section retrieval from lhefile currently not implemented" )
                sys.exit()
            else:
                logger.info("Writing pythia LHE output to " + lhefile)
        if loFromSlha:
            logger.info("Using LO cross sections from " + slhafile)
            xsecsInfile = crossSection.getXsecFromSLHAFile(slhafile)
            loXsecs = crossSection.XSectionList()
            for xsec in xsecsInfile:
                if xsec.info.order == 0 and xsec.info.sqrts == sqrts:
                    loXsecs.add(xsec)
                    
        else:
            logger.info("get LO cross sections from pythia%d" % self.pythiaVersion )
            tool = toolBox.ToolBox().get("pythia%d" % self.pythiaVersion )
            tool.nevents = self.nevents
            tool.sqrts = sqrts / TeV
            tool.pythiacard = pythiacard
            loXsecs = tool.run( slhafile, lhefile, unlink=unlink )
        self.loXsecs = loXsecs
        self.loXsecs.sort()
        self.xsecs = self.addHigherOrders ( sqrts, slhafile )
        self.xsecs.sort()
        #for xsec in self.loXsecs:
        #    logger.debug ( "now writing out xsecs: %s" % xsec.value )
        logger.debug ( "how many NLL xsecs? %d" % len(self.xsecs) )
        return self.xsecs

    def computeForOneFile ( self, sqrtses, inputFile, unlink,
                            lOfromSLHA, tofile, pythiacard=None ):
        """
        compute the cross sections for one file.
        :param sqrtses: list of sqrt{s} tu run pythia, as a unum (e.g. 7*TeV)
            
        """
        if tofile:
            logger.info("Computing SLHA cross section from %s, adding to "
                        "SLHA file." % inputFile )
            for s in sqrtses:
                ss = s*TeV 
                self.compute( ss, inputFile, unlink= unlink, 
                              loFromSlha= lOfromSLHA, pythiacard=pythiacard )
                if tofile == "all":
                    comment = str(self.nevents)+" evts, pythia%d [pb]"%\
                                              self.pythiaVersion
                    self.addXSecToFile(self.loXsecs, inputFile, comment )
                comment = str(self.nevents)+" events, [pb], pythia%d for LO"%\
                                              self.pythiaVersion
                self.addXSecToFile( self.xsecs, inputFile, comment)
        else:
            logger.info("Computing SLHA cross section from %s." % inputFile )
            print()
            print( "     Cross sections:" )
            print( "=======================" )
            for s in sqrtses:
                ss = s*TeV 
                self.compute( ss, inputFile, unlink=unlink, loFromSlha=lOfromSLHA )
                for xsec in self.xsecs: 
                    print( "%s %20s:  %.3e pb" % \
                            ( xsec.info.label,xsec.pid,xsec.value/pb ) )
            print()

    def computeForBunch ( self, sqrtses, inputFiles, unlink,
                            lOfromSLHA, tofile, pythiacard=None ):
        """ compute xsecs for a bunch of slha files """
        # computer = XSecComputer( order, nevents, pythiaVersion )
        for inputFile in inputFiles:
            logger.debug ( "computing xsec for %s" % inputFile )
            self.computeForOneFile ( sqrtses, inputFile, unlink, lOfromSLHA, 
                                     tofile, pythiacard=pythiacard )

    def addXSecToFile( self, xsecs, slhafile, comment=None, complain=True):
        """
        Write cross sections to an SLHA file.
        
        :param xsecs: a XSectionList object containing the cross sections
        :param slhafile: target file for writing the cross sections in SLHA format
        :param comment: optional comment to be added to each cross section block
        :param complain: complain if there are already cross sections in file
        
        """
        
        if not os.path.isfile(slhafile):
            logger.error("SLHA file not found.")
            raise SModelSError()
        if len(xsecs) == 0:
            logger.warning("No cross sections available.")
            return False
        # Check if file already contain cross section blocks
        xSectionList = crossSection.getXsecFromSLHAFile(slhafile)
        if xSectionList and complain:
            logger.info("SLHA file already contains XSECTION blocks. Adding "
                           "only missing cross sections.")

        # Write cross sections to file, if they do not overlap any cross section in
        # the file
        outfile = open(slhafile, 'a')
        for xsec in xsecs:
            writeXsec = True
            for oldxsec in xSectionList:
                if oldxsec.info == xsec.info and set(oldxsec.pid) == set(xsec.pid):
                    writeXsec = False
                    break
            if writeXsec:
                outfile.write( self.xsecToBlock(xsec, (2212, 2212), comment) + "\n")
        outfile.close()

        return True

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
                str( "%16.8E" % (xsec.value / xsecUnit) ) + " SModelSv" + \
                     installation.version()

        return "\n" + header + "\n" + entry


class ArgsStandardizer:
    """ simple class to collect all argument manipulators """

    def getInputFiles ( self, args ):
        """ geth the names of the slha files to run over """
        inputPath  = args.filename.strip()
        if not os.path.exists( inputPath ):
            logger.error( "Path %s does not exist." % inputPath )
            sys.exit(1)
        inputFiles = []
        if os.path.isfile ( inputPath ):
            inputFiles = [ inputPath ]
        else:
            files = os.listdir ( inputPath )
            for f in files:
                inputFiles.append ( os.path.join ( inputPath, f ) )
        return inputFiles

    def checkAllowedSqrtses ( self, order, sqrtses ):
        """ check if the sqrtses are 'allowed' """
        if order == 0: return
        allowedsqrtses=[7, 8, 13]
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
        sqrtses = [item for sublist in args.sqrts for item in sublist]
        if len(sqrtses) == 0:
            sqrtses = [8,13]
        sqrtses.sort()
        sqrtses = set(sqrtses)
        return sqrtses

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
                logger.warn ( "Specified both --tofile and --alltofile. Will use "\
                              "--alltofile" )
            toFile="all"
        return toFile

def main(args):
    canonizer = ArgsStandardizer()
    setLogLevel ( args.verbosity )
    if args.query:
        return canonizer.queryCrossSections ( args.filename )
    sqrtses = canonizer.getSqrtses ( args )
    order = canonizer.getOrder ( args )
    canonizer.checkAllowedSqrtses ( order, sqrtses )
    inputFiles = canonizer.getInputFiles ( args )
    ncpus = canonizer.checkNCPUs ( args.ncpus, inputFiles )
    pythiaVersion = canonizer.getPythiaVersion ( args )

    pythiacard = None
    if hasattr(args, 'pythiacard'):
        pythiacard = args.pythiacard

    children = []
    for i in range(ncpus):
        pid = os.fork()
        chunk = inputFiles [ i::ncpus ]
        if pid < 0:
            logger.error ( "fork did not succeed! Pid=%d" % pid ) 
            sys.exit()
        if pid == 0:
            logger.debug ( "chunk #%d: pid %d (parent %d)." % 
                       ( i, os.getpid(), os.getppid() ) )
            logger.debug ( " `-> %s" % " ".join ( chunk ) )
            computer = XSecComputer( order, args.nevents, pythiaVersion )
            toFile = canonizer.writeToFile ( args )
            computer.computeForBunch (  sqrtses, chunk, not args.keep,
                                args.LOfromSLHA, toFile, pythiacard=pythiacard )
            os._exit ( 0 )
        if pid > 0:
            children.append ( pid )
    for child in children:
        r = os.waitpid ( child, 0 )
        logger.debug ( "child %d terminated: %s" % (child,r) )
    logger.debug ( "all children terminated." )
