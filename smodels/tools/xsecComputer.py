#!/usr/bin/env python

"""
.. module:: xsecComputer
   :synopsis: Computation of reference ("theory") production cross sections.

.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
from smodels import installation
from smodels.tools import toolBox
from smodels.tools.physicsUnits import fb, TeV
from smodels.theory import crossSection
from smodels.tools import nllFast
import os
import cStringIO
import logging
import argparse
import types
import sys

logger = logging.getLogger(__name__)

LO=0
NLO=1
NLL=2

def computeXSec(sqrts, maxOrder, nevts, slhafile, lhefile=None, unlink=True, loFromSlha=None ):
    """
    Run pythia and compute SUSY cross-sections for the input SLHA file.

    :param sqrts: sqrt{s} to run Pythia
    :param maxOrder: maximum order to compute the cross-section
                if maxOrder = 0, compute only LO pythia xsecs
                if maxOrder = 1, apply NLO K-factors from NLLfast (if available)
                if maxOrder = 2, apply NLO+NLL K-factors from NLLfast (if available)
    :param nevts: number of events for pythia run
    :param slhafile: SLHA file
    :param lhefile: LHE file. If None, do not write pythia output to file. If
                    file does not exist, write pythia output to this file name. If
                    file exists, read LO xsecs from this file (does not run pythia).
    :param unlink: Clean up temp directory after running pythia

    :param loFromSlha: If True, uses the LO xsecs from the SLHA file to compute the
                       higher order xsecs

    :returns: XSectionList object

    """

    if not os.path.isfile(slhafile):
        logger.error("SLHA file %s not found.", slhafile)
        sys.exit()
    if lhefile:
        if os.path.isfile(lhefile):
            logger.warning("Using LO cross-sections from " + lhefile)
        else:
            logger.info("Writing pythia LHE output to " + lhefile)
    if loFromSlha:
        logger.info("Using LO cross-sections from " + slhafile)
        xsecsInfile = crossSection.getXsecFromSLHAFile(slhafile)
        loXsecs = crossSection.XSectionList()
        for xsec in xsecsInfile:
            if xsec.info.order == 0: loXsecs.add(xsec)
    else:
        if not lhefile or not os.path.isfile(lhefile):
            lheFile = runPythia(slhafile, nevts, sqrts / TeV, lhefile, unlink=unlink)
        else:
            lheFile = open(lhefile, 'r')
        loXsecs = crossSection.getXsecFromLHEFile(lheFile)

    xsecs = loXsecs
    wlabel = str(int(sqrts / TeV)) + ' TeV'
    if maxOrder == 0:
        wlabel += ' (LO)'
    elif maxOrder == 1:
        wlabel += ' (NLO)'
    elif maxOrder == 2:
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
            elif maxOrder > 2:
                logger.warning("Unkown xsec order, using NLL+NLO k-factor, "
                               "if available")
                k = kNLO * kNLL
            k = float(k)
            for i, xsec in enumerate(xsecs):
                if set(xsec.pid) == set(pID):
                    # Apply k-factor
                    xsecs[i] = xsec * k

    # Remove zero cross-sections
    while len(xsecs) > 0 and xsecs.getMinXsec() == 0. * fb:
        for xsec in xsecs:
            if xsec.value == 0. * fb:
                xsecs.delete(xsec)
                break
    if maxOrder > 0 and len(xsecs) == 0:
        logger.warning("No NLO or NLL cross-sections available.")
    return xsecs


def addXSecToFile(xsecs, slhafile, comment=None, complain=True):
    """
    Write cross-sections to an SLHA file.
    
    :param xsecs: a XSectionList object containing the cross-sections
    :param slhafile: target file for writing the cross-sections in SLHA format
    :param comment: optional comment to be added to each cross-section block
    :param complain: complain if there are already cross sections in file

    """
    if not os.path.isfile(slhafile):
        logger.error("SLHA file not found.")
        import sys
        sys.exit()
    if len(xsecs) == 0:
        logger.warning("No cross-sections available.")
        import sys
        sys.exit()
    # Check if file already contain cross-section blocks
    xSectionList = crossSection.getXsecFromSLHAFile(slhafile)
    if xSectionList and complain:
        logger.warning("SLHA file already contains XSECTION blocks. Adding "
                       "only missing cross-sections.")

    # Write cross-sections to file, if they do not overlap any cross-section in
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


def xsecToBlock(xsec, inPDGs=(2212, 2212), comment=None):
    """
    Generate a string for a XSECTION block in the SLHA format from a XSection
    object.

    :param inPDGs: defines the PDGs of the incoming states
                   (default = 2212,2212)

    :param comment: is added at the end of the header as a comment

    """
    if type(xsec) != type(crossSection.XSection()):
        logger.error("Wrong input")
        import sys
        sys.exit()
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
    entry = "0  " + str(xsec.info.order) + "  0  0  0  0  " + \
            str("%16.8E" % (xsec.value / fb) ) + " SModelS " + installation.version()

    return "\n" + header + "\n" + entry


def runPythia(slhafile, nevts, sqrts, lhefile=None, unlink=True ):
    """
    Execute pythia_lhe with n events, at sqrt(s)=sqrts.

    :param slhafile: input SLHA file
    :param nevts: number of events to be generated
    :param sqrts: center of mass sqrt{s} (in TeV)
    :param lhefile: option to write LHE output to file; ff None, do not write
                    output to disk.
    :param unlink: Clean up temp directory after running pythia
    :returns: file object with the LHE events

    """
    box = toolBox.ToolBox()
    tool = box.get("pythia6")
    # Check if template config file exists
    tool.reset()
    tool.replaceInCfgFile({"NEVENTS": nevts, "SQRTS":1000 * sqrts})
    tool.setParameter("MSTP(163)", "6")

    if unlink==False:
        logger.info ( "keeping temporary directory at %s" % tool.tempDirectory() )
    lhedata = tool.run(slhafile, do_unlink=unlink )
    if not "<LesHouchesEvents" in lhedata:
        logger.error("LHE events not found in pythia output")
        import sys
        sys.exit()

    # Generate file object with lhe events
    if lhefile:
        lheFile = open(lhefile, 'w')
        lheFile.write(lhedata)
        lheFile.close()
        lheFile = open(lhefile, 'r')
    else:
        # Create memory only file object
        lheFile = cStringIO.StringIO(lhedata)

    return lheFile


if __name__ == "__main__":
    """
    Compute the cross section of a given SLHA file.

    """
    desc = "compute the cross section of an SLHA file"
    argparser = argparse.ArgumentParser(description=desc)
    argparser.add_argument('file', type=types.StringType, nargs=1,
                           help="SLHA file to compute cross section "
                           "for")
    argparser.add_argument('-s', '--sqrts', nargs='+', action='append',
                           help="sqrt(s) TeV. Can supply more than one value.",
                           type=int, default=[])
    argparser.add_argument('-e', '--nevents', type=int, default=100,
                           help="number of events to be simulated.")
    argparser.add_argument('-f', '--tofile', action='store_true',
                           help="write cross sections to file")
    #argparser.add_argument('-S', '--slha', action='store_true',
    #                       help="input file is an SLHA file")
    #argparser.add_argument('-L', '--lhe', action='store_true',
    #                       help="input file is an LHE file")
    argparser.add_argument('-k', '--keep', action='store_true',
                           help="do not unlink temporary directory")
    argparser.add_argument('-n', '--NLO', action='store_true',
                           help="compute at the NLO level (default is LO)")
    argparser.add_argument('-N', '--NLL',
                           help="compute at the NLL level (takes precedence "
                           "over NLL, default is LO)",
                           action='store_true')
    args = argparser.parse_args()



    sqrtses = [item for sublist in args.sqrts for item in sublist]
    if len(sqrtses) == 0:
        sqrtses = [8]
    sqrtses.sort()
    sqrtses = set(sqrtses)
    order = 0
    if args.NLO:
        order = 1
    if args.NLL:
        order = 2
    if order > 0:
        for sqrts in sqrtses:
            allowedsqrtses=[7, 8, 13, 14, 33, 100]
            if not sqrts in allowedsqrtses:
                logger.error("Cannot compute NLO or NLL xsecs for sqrts = %d "
                             "TeV! Available are: %s TeV." % 
                             (sqrts, allowedsqrtses ))
                sys.exit(0)
                # sqrtses.remove(sqrts)
    inputFile = args.file[0]
    if not os.path.exists(inputFile):
        logger.error("File '%s' does not exist.", inputFile)
        sys.exit(1)
    #if inputFile[-5:].lower() == ".slha" or args.slha:
    if args.tofile:
        logger.info("Computing SLHA cross section from %s, adding to "
                    "SLHA file." % inputFile )
        for s in sqrtses:
            ss = s*TeV 
            xsecs = computeXSec( ss, order, args.nevents, inputFile, 
                                 unlink=(not args.keep) )
            comment = "Nevts: " + str(args.nevents) + " xsec unit: fb"
            addXSecToFile(xsecs, inputFile, comment)
        sys.exit(0)
    else:
        logger.info("Computing SLHA cross section from %s." % inputFile )
        print
        print "     Cross sections:"
        print "======================="
        for s in sqrtses:
            ss = s*TeV 
            xsecs = computeXSec(ss, order, args.nevents, inputFile, \
                        unlink=(not args.keep) )
            for xsec in xsecs: 
                print "%s %20s:  %5.2f fb" % ( xsec.info.label,xsec.pid,xsec.value/fb )
        print
        sys.exit(0)
    #if inputFile[-4:].lower() == ".lhe" or args.lhe:
    #    if args.tofile:
    #        logger.error("Compute LHE section, and add to file. TODO: I guess "
    #                     "we do not need this case?")
    #        sys.exit(0)
    #    else:
    #        logger.error("Compute LHE section, print out, but do not add to "
    #                     "file. TODO: not yet implemented.")
    #        sys.exit(0)

