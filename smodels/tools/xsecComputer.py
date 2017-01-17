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
from smodels.tools import nllFast
from smodels.tools.smodelsLogging import logger, setLogLevel
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
import os
import pyslha
try:
    import cStringIO as io
except ImportError as e:
    import io
import sys

LO= 0  ## simple variables used to increase readability the perturbation order 
NLO=1
NLL=2

def computeXSec(sqrts, maxOrder, nevts, slhafile, lhefile=None, unlink=True, loFromSlha=None, pythiacard=None ):
    """
    Run pythia and compute SUSY cross sections for the input SLHA file.

    :param sqrts: sqrt{s} to run Pythia, given as a unum (e.g. 7.*TeV)
    :param maxOrder: maximum order to compute the cross section, given as an integer
                if maxOrder == 0, compute only LO pythia xsecs
                if maxOrder == 1, apply NLO K-factors from NLLfast (if available)
                if maxOrder == 2, apply NLO+NLL K-factors from NLLfast (if available)
    :param nevts: number of events for pythia run
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
    if not os.path.isfile(slhafile):
        logger.error("SLHA file %s not found.", slhafile)
        raise SModelSError()
    try:
        f=pyslha.readSLHAFile(slhafile)
    except pyslha.ParseError as e:
        logger.error("File cannot be parsed as SLHA file: %s" % e )
        raise SModelSError()

    if type(sqrts)==type(float) or type(sqrts)==type(int):
        logger.warning("sqrt(s) given as scalar, will add TeV as unit." )
        sqrts=float(sqrts)*TeV

    smaxorder={ "LO": 0, "NLO": 1, "NLL": 2 }
    if maxOrder in smaxorder.keys():
        logger.warning("maxorder given as string, please supply integer.")
        maxOrder=smaxorder[maxOrder]

    if lhefile:
        if os.path.isfile(lhefile):
            logger.warning("Using LO cross sections from " + lhefile)
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
        if not lhefile or not os.path.isfile(lhefile):
            lheFile = runPythia(slhafile, nevts, sqrts / TeV, lhefile, unlink=unlink,pythiacard=pythiacard)
        else:
            lheFile = open(lhefile, 'r')
        loXsecs = crossSection.getXsecFromLHEFile(lheFile)
    xsecs = loXsecs
    wlabel = str(int(sqrts / TeV)) + ' TeV'
    if maxOrder == 0:
        wlabel += ' (LO)'
    elif maxOrder == 1:
        wlabel += ' (NLO)'
    elif maxOrder >= 2:
        wlabel += ' (NLO+NLL)'
    for ixsec, xsec in enumerate(xsecs):
        xsecs[ixsec].info.label = wlabel
        xsecs[ixsec].info.order = maxOrder

    if maxOrder > 0:
        pIDs = loXsecs.getPIDpairs()
        for pID in pIDs:
            k = 0.
            kNLO, kNLL = nllFast.getKfactorsFor(pID, sqrts, slhafile)
            if maxOrder == 1 and kNLO:
                k = kNLO
            elif maxOrder == 2 and kNLL and kNLO:
                k = kNLO * kNLL
            elif maxOrder > 2 and kNLL and kNLO:
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
    if maxOrder > 0 and len(xsecs) == 0:
        logger.warning("No NLO or NLL cross sections available.")
        
    return xsecs


def addXSecToFile(xsecs, slhafile, comment=None, complain=True):
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
            outfile.write(xsecToBlock(xsec, (2212, 2212), comment) + "\n")
    outfile.close()

    return True


def xsecToBlock(xsec, inPDGs=(2212, 2212), comment=None, xsecUnit = pb):
    """
    Generate a string for a XSECTION block in the SLHA format from a XSection
    object.

    :param inPDGs: defines the PDGs of the incoming states
                   (default = 2212,2212)

    :param comment: is added at the end of the header as a comment
    :param xsecUnit: unit of cross sections to be written (default is pb). Must be a Unum unit.

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
        header += "   # " + str(comment)  # Comment
    entry = "  0  " + str(xsec.info.order) + "  0  0  0  0  " + \
            str("%16.8E" % (xsec.value / xsecUnit) ) + " SModelS " + installation.version()

    return "\n" + header + "\n" + entry


def runPythia(slhafile, nevts, sqrts, lhefile=None, unlink=True, pythiacard=None ):
    """
    Execute pythia_lhe with n events, at sqrt(s)=sqrts.

    :param slhafile: input SLHA file
    :param nevts: number of events to be generated
    :param sqrts: center of mass sqrt{s} (in TeV)
    :param lhefile: option to write LHE output to file; ff None, do not write
                    output to disk.
    :param unlink: Clean up temp directory after running pythia
    :param pythiaCard: Optional path to pythia.card. If None, uses /etc/pythia.card
    
    :returns: file object with the LHE events

    """
    box = toolBox.ToolBox()
    tool = box.get("pythia6")
    #Change pythia card, if defined:
    if pythiacard:
        pythiacard_default = tool.cfgfile
        tool.cfgfile = pythiacard
    # Check if template config file exists
    tool.unlink()
    tool.replaceInCfgFile({"NEVENTS": nevts, "SQRTS":1000 * sqrts})
    tool.setParameter("MSTP(163)", "6")

    if unlink==False:
        logger.info ( "keeping temporary directory at %s" % tool.tempDirectory() )
    r = tool.checkInstallation()
    if r == False:
        logger.info ( "Installation check failed." )
        sys.exit()
    #logger.info ( "cfgfile=%s" % tool.cfgfile )
    #logger.info ( "executable=%s" % tool.executable )
    #logger.info ( "tempdir=%s" % tool.tempdir )
    #logger.info ( "nevents=%s" % tool.nevents )
    tool.replaceInCfgFile({"NEVENTS": nevts, "SQRTS":1000 * sqrts})
    tool.setParameter("MSTP(163)", "6")
    lhedata = tool.run(slhafile, do_check=False, do_unlink=unlink )
    if not "<LesHouchesEvents" in lhedata:
        pythiadir = "%s/log" % tool.tempDirectory()
        logger.error("No LHE events found in pythia output %s" % pythiadir )
        if not os.path.exists ( pythiadir ):
            logger.error ("Will dump pythia output to %s" % pythiadir )
            f=open ( pythiadir, "w" )
            for line in lhedata:
                f.write ( line )
            f.close()
        raise SModelSError( "No LHE events found in %s" % pythiadir )

    #Reset pythia card to its default value
    if pythiacard:
        tool.cfgfile = pythiacard_default

    # Generate file object with lhe events
    if lhefile:
        lheFile = open(lhefile, 'w')
        lheFile.write(lhedata)
        lheFile.close()
        lheFile = open(lhefile, 'r')
    else:
        # Create memory only file object
        lheFile = io.StringIO(lhedata)

    return lheFile

def computeForOneFile ( sqrtses, order, nevents, inputFile, unlink,
                        lOfromSLHA, tofile, pythiacard=None ):
    """ compute the cross sections for one file """
    if tofile:
        logger.info("Computing SLHA cross section from %s, adding to "
                    "SLHA file." % inputFile )
        for s in sqrtses:
            ss = s*TeV 
            xsecs = computeXSec( ss, order, nevents, inputFile, 
                                 unlink= unlink, loFromSlha= lOfromSLHA, pythiacard=pythiacard)
            comment = "Nevts: " + str(nevents) + " xsec unit: pb"
            addXSecToFile(xsecs, inputFile, comment)
    else:
        logger.info("Computing SLHA cross section from %s." % inputFile )
        print()
        print( "     Cross sections:" )
        print( "=======================" )
        for s in sqrtses:
            ss = s*TeV 
            xsecs = computeXSec(ss, order, nevents, inputFile, \
                        unlink=unlink, loFromSlha=lOfromSLHA )
            for xsec in xsecs: 
                print( "%s %20s:  %.3e pb" % ( xsec.info.label,xsec.pid,xsec.value/pb ) )
        print()

def queryCrossSections ( filename ):
    if os.path.isdir ( filename ):
        logger.error ( "Cannot query cross sections for a directory." )
        sys.exit(-1)
    xsecsInfile = crossSection.getXsecFromSLHAFile(filename)
    if xsecsInfile:
        print ( "1" )
    else:
        print ( "0" )

def getOrder ( args ):
    """ retrieve the order in perturbation theory from argument list """
    if args.NLL:
        return 2
    if args.NLO:
        return 1
    return 0

def getSqrtses ( args ):
    """ extract the sqrtses from argument list """
    sqrtses = [item for sublist in args.sqrts for item in sublist]
    if len(sqrtses) == 0:
        sqrtses = [8,13]
    sqrtses.sort()
    sqrtses = set(sqrtses)
    return sqrtses

def checkAllowedSqrtses ( order, sqrtses ):
    """ check if the sqrtses are 'allowed' """
    if order == 0: return
    allowedsqrtses=[7, 8, 13]
    for sqrts in sqrtses:
        if not sqrts in allowedsqrtses:
            logger.error("Cannot compute NLO or NLL xsecs for sqrts = %d "
                    "TeV! Available are: %s TeV." % (sqrts, allowedsqrtses ))
            sys.exit(-2)

def computeForBunch ( sqrtses, order, nevents, inputFiles, unlink,
                        lOfromSLHA, tofile, pythiacard=None ):
    """ compute xsecs for a bunch of slha files """
    for inputFile in inputFiles:
        logger.debug ( "computing xsec for %s" % inputFile )
        computeForOneFile ( sqrtses, order, nevents, inputFile, unlink, 
                            lOfromSLHA, tofile, pythiacard=pythiacard )

def getInputFiles ( args ):
    """ geth the names of the slha files to run over """
    inputPath  = args.filename.strip()
    if not os.path.exists( inputPath ):
        logger.error("Path '%s' does not exist.", inputFile)
        sys.exit(1)
    inputFiles = []
    if os.path.isfile ( inputPath ):
        inputFiles = [ inputPath ]
    else:
        files = os.listdir ( inputPath )
        for f in files:
            inputFiles.append ( os.path.join ( inputPath, f ) )
    return inputFiles

def main(args):
    setLogLevel ( args.verbosity )
    if args.query:
        return queryCrossSections ( args.filename )
    sqrtses = getSqrtses ( args )
    order = getOrder ( args )
    checkAllowedSqrtses ( order, sqrtses )
    inputFiles = getInputFiles ( args )
    ncpus = args.ncpus

    if hasattr(args, 'pythiacard'):
        pythiacard = args.pythiacard
    else:
        pythiacard = None
    if ncpus < -1 or ncpus == 0:
        logger.error ( "Weird number of CPUs given: %d" % ncpus )
        sys.exit()
    if ncpus == -1:
        ncpus = runtime.nCPUs()
    logger.info ( "We run on %d cpus" % ncpus )
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
            computeForBunch (  sqrtses, order, args.nevents, chunk, not args.keep,
                               args.LOfromSLHA, args.tofile, pythiacard=pythiacard)
            os._exit ( 0 )
        if pid > 0:
            children.append ( pid )
    for child in children:
        r = os.waitpid ( child, 0 )
        logger.debug ( "child %d terminated: %s" % (child,r) )
    logger.debug ( "all children terminated." )
