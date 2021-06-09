#!/usr/bin/env python3

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
    def __init__ ( self, maxOrder, nevents, pythiaVersion, maycompile=True ):
        """
        :param maxOrder: maximum order to compute the cross section, given as an integer
                    if maxOrder == LO, compute only LO pythia xsecs
                    if maxOrder == NLO, apply NLO K-factors from NLLfast (if available)
                    if maxOrder == NLL, apply NLO+NLL K-factors from NLLfast (if available)
        :param nevents: number of events for pythia run
        :param pythiaVersion: pythia6 or pythia8 (integer)
        :param maycompile: if True, then tools can get compiled on-the-fly
        """
        self.maxOrder = self._checkMaxOrder ( maxOrder )
        self.countNoXSecs = 0
        self.countNoNLOXSecs = 0
        self.maycompile = maycompile
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
        except (pyslha.ParseError,ValueError) as e:
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
            nllfast = toolBox.ToolBox().get("nllfast%d" % sqrts.asNumber(TeV) )
            nllfast.maycompile = self.maycompile
            for pID in pIDs:
                k = 0.
                kNLO, kNLL = nllfast.getKfactorsFor(pID, slhafile)
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
            self.countNoNLOXSecs += 1
            if self.countNoNLOXSecs < 3:
                logger.warning("No NLO or NLL cross sections available.")
            if self.countNoNLOXSecs == 3:
                logger.warning("No NLO or NLL cross sections available (will quench such warnings in future).")

        #for i in xsecs:
        #    logger.error ( "xsec=%s (%s)" % (i,type(i)) )
        return xsecs

    def match ( self, pids, theorypid ):
        """ do the pids given by the user match the
            pids of the theorypred? """
        spids = list(pids)
        stpids = list(theorypid)
        if len(spids)!=len(stpids):
            return False
        jokerpids = []
        for spid in spids:
            if type(spid)==int:
                if not spid in stpids:
                    return False
                else:
                    stpids.remove(spid)
            else:
                if type(spid)!=str:
                    logger.error ( "I have a pid of type %s. Dont know what to do." % \
                                    type(spid) )
                    sys.exit()
                jokerpids.append ( spid )
        jokerpids.sort( key=len, reverse=True ) # the longer, the more constraining
        if stpids == []: # no wildcards
            #print ( "tpid", theorypid, "matched", pids, "no wildcards" )
            return True
        import fnmatch
        for jpid in jokerpids:
            nomatch=[]
            hasMatched = False ## one jpid should match only one
            for i,stpid in enumerate(stpids):
                if not fnmatch.fnmatch ( str(stpid), jpid ):
                    nomatch.append ( stpid )
                else:
                    for j in stpids[i+1:]:
                        nomatch.append ( j )
                    break
            #print ( "jpid", jpid," no matches", nomatch )
            stpids = nomatch
        if len(stpids)>0:
            return False
        #print ( "tpid", theorypid, "matched", pids )
        return True


    def applyMultipliers ( self, xsecs, ssmultipliers ):
        """
        apply the given multipliers to the cross sections """
        for pids in ssmultipliers.keys():
            if type(pids) not in [ list, tuple ]:
                logger.error ( "signal strength multipliers need to be supplied as a dictionary, with the keys being tuples of the mothers' pids, e.g. { (1000021, -1000001 ): 0.9, .... }" )
                sys.exit()
            if len(pids) != 2:
                logger.warning ( "currently we always only have two mothers, so why are the signal strength multipliers given for %d mothers?" % len(pids) )
            known_pids = [ 3000006 ] ## here we can define some exceptional pids
            for pid in pids:
                if type(pid)==int and (abs(pid) < 1000000 or ( abs(pid) > 3000000 and abs(pid) not in known_pids) ):
                    logger.warning ( "signal strength multiplier for pid %d supplied. what does that mean?" % pid )
        for x in xsecs:
            for pids, multiplier in ssmultipliers.items():
                if self.match ( pids, x.pid ):
                    x.value = x.value * multiplier
                    break
        # return xsecs
        newx = crossSection.XSectionList()
        for x in xsecs:
            if x.value.asNumber(pb)>0.:
                newx.add(x)
        return newx

    def getPythia ( self ):
        """ returns the pythia tool that is configured to be used """
        ret= toolBox.ToolBox().get("pythia%d" % self.pythiaVersion )
        ret.maycompile = self.maycompile
        return ret

    def compute ( self, sqrts, slhafile,  lhefile=None, unlink=True, loFromSlha=None,
                  pythiacard=None, ssmultipliers=None ):
        """
        Run pythia and compute SUSY cross sections for the input SLHA file.

        :param sqrts: sqrt{s} to run Pythia, given as a unum (e.g. 7.*TeV)
        :param slhafile: SLHA file
        :param lhefile: LHE file. If None, do not write pythia output to file. If
                     file does not exist, write pythia output to this file name. If
                     file exists, read LO xsecs from this file (does not run pythia).
        :param unlink: Clean up temp directory after running pythia
        :param loFromSlha: If True, uses the LO xsecs from the SLHA file to
                           compute the higher order xsecs
        :param pythiaCard: Optional path to pythia.card. If None, uses
                           smodels/etc/pythia.card
        :param ssmultipliers: optionally supply signal strengh multipliers,
                given as dictionary of the tuple of the mothers' pids as keys and
                multipliers as values, e.g { (1000001,1000021):1.1 }.
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
                if xsec.info.order == 0 and abs ( xsec.info.sqrts - sqrts ).asNumber(TeV) < 1e-5:
                    loXsecs.add(xsec)
            if ssmultipliers != None:
                logger.warning ( "supplied signal strength multipliers, but LO xsecs are taken from file. Will not apply them" )
        else:
            logger.info("get LO cross sections from pythia%d" % self.pythiaVersion )
            tool = self.getPythia()
            tool.tempdir = None # reset, to make sure it works in parallel mode
            tool.nevents = self.nevents
            tool.sqrts = sqrts / TeV
            tool.pythiacard = pythiacard
            loXsecs = tool.run( slhafile, lhefile, unlink=unlink )
            if ssmultipliers != None:
                loXsecs = self.applyMultipliers ( loXsecs, ssmultipliers )
        self.loXsecs = loXsecs
        self.loXsecs.sort()
        self.xsecs = self.addHigherOrders ( sqrts, slhafile )
        self.xsecs.sort()
        #for xsec in self.loXsecs:
        #    logger.debug ( "now writing out xsecs: %s" % xsec.value )
        logger.debug ( "how many NLL xsecs? %d" % len(self.xsecs) )
        return self.xsecs

    def computeForOneFile ( self, sqrtses, inputFile, unlink,
                            lOfromSLHA, tofile, pythiacard = None,
                            ssmultipliers = None, comment = None ):
        """
        Compute the cross sections for one file.

        :param sqrtses: list of sqrt{s} tu run pythia, as a unum (e.g. [7*TeV])
        :param inputFile: input SLHA file to compute xsecs for
        :param unlink: if False, keep temporary files
        :param lofromSLHA: try to obtain LO xsecs from SLHA file itself
        :param tofile: False, True, "all": write results to file, if "all" also write lower xsecs to file.
        :param pythiacard: optionally supply your own runcard
        :param ssmultipliers: optionally supply signal strengh multipliers,
                              given as dictionary of the tuple of the mothers' pids as keys and
                              multipliers as values, e.g { (1000001,1000021):1.1 }.
        :param comment: an optional comment that gets added to the slha file.

        :returns: number of xsections that have been computed
        """
        
        nXSecs = 0 ## count the xsecs we are adding
        if tofile:
            logger.info("Computing SLHA cross section from %s, adding to "
                        "SLHA file." % inputFile )
            complain = True ## dont complain about already existing xsecs,
            # if we were the ones writing them
            for s in sqrtses:
                ss = s*TeV
                self.compute( ss, inputFile, unlink= unlink, loFromSlha= lOfromSLHA,
                              pythiacard=pythiacard, ssmultipliers = ssmultipliers )
                if tofile == "all":
                    xcomment = str(self.nevents)+" evts, pythia%d [pb]"%\
                                              self.pythiaVersion
                    nXSecs += self.addXSecToFile(self.loXsecs, inputFile, xcomment, complain )
                    complain = False
                xcomment = str(self.nevents)+" events, [pb], pythia%d for LO"%\
                                              self.pythiaVersion
                if tofile != False:
                    nXSecs += self.addXSecToFile( self.xsecs, inputFile, xcomment, complain)
                    complain = False
            if nXSecs > 0: ## only add if we actually added xsecs
                self.addMultipliersToFile ( ssmultipliers, inputFile )
            self.addCommentToFile ( comment, inputFile )
        else:
            logger.info("Computing SLHA cross section from %s." % inputFile )
            print()
            print( "     Cross sections:" )
            print( "=======================" )
            for s in sqrtses:
                ss = s*TeV
                self.compute( ss, inputFile, unlink=unlink, loFromSlha=lOfromSLHA,
                              ssmultipliers = ssmultipliers )
                for xsec in self.xsecs:
                    nXSecs += 1
                    print( "%s %20s:  %.3e pb" % \
                            ( xsec.info.label,xsec.pid,xsec.value/pb ) )
            print()
        return nXSecs

    def computeForBunch ( self, sqrtses, inputFiles, unlink,
                          lOfromSLHA, tofile, pythiacard=None,
                          ssmultipliers=None   ):
        """ compute xsecs for a bunch of slha files """
        for inputFile in inputFiles:
            logger.debug ( "computing xsec for %s" % inputFile )
            self.computeForOneFile ( sqrtses, inputFile, unlink, lOfromSLHA,
                      tofile, pythiacard=pythiacard, ssmultipliers = ssmultipliers )

    def addCommentToFile ( self, comment, slhaFile ):
        """ add the optional comment to file """
        if comment in [ None, "" ]:
            return
        if not os.path.isfile(slhaFile ):
            logger.error("SLHA file %s not found." % slhaFile )
            raise SModelSError()
        outfile = open(slhaFile, 'a')
        outfile.write ( "# %s\n" % comment )
        outfile.close()

    def addMultipliersToFile ( self, ssmultipliers, slhaFile ):
        """ add the signal strength multipliers to the SLHA file """
        if ssmultipliers in [ None, {} ]:
            return
        if not os.path.isfile(slhaFile ):
            logger.error("SLHA file %s not found." % slhaFile )
            raise SModelSError()
        tokens = []
        for k,v in ssmultipliers.items():
            tokens.append ( "%s:%.4g" % ( k, v ) )
        newline = "# Signal strength multipliers: " + ", ".join ( tokens )
        with open(slhaFile, 'r' ) as r:
            lines = r.readlines()
            r.close()

        rewrite = []
        for line in lines:
            if "Signal strength multipliers" in line:
                if ( line.strip() == newline ):
                    logger.debug ( "Signal strength multipliers have alread been applied." )
                else:
                    logger.error ( "Different signal strength multipliers have alread been applied!!!" )
                    rewrite.append ( line+" ERROR inconsistent!" )
            else:
                if not "produced at step" in line:
                    rewrite.append ( line )
        outfile = open(slhaFile, 'w')
        for line in rewrite:
            outfile.write ( line )
        if line != "\n": ## last line not an empty newline?
            outfile.write ( "\n" )
        outfile.write ( newline )
        outfile.write ( "\n" )
        outfile.close()

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

        # Write cross sections to file, if they do not overlap any cross section in
        # the file
        outfile = open(slhafile, 'a')
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
                str( "%16.8E" % (xsec.value / xsecUnit) ) + " SModelSv" + \
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
                logger.warning ( "Specified both --tofile and --alltofile. Will use "\
                              "--alltofile" )
            toFile="all"
        return toFile

def main(args):
    canonizer = ArgsStandardizer()
    setLogLevel ( args.verbosity )
    if not hasattr ( args, "noautocompile" ):
        args.noautocompile = False
    if args.query:
        return canonizer.queryCrossSections ( args.filename )
    if args.colors:
        from smodels.tools.colors import colors
        colors.on = True
    sqrtses = canonizer.getSqrtses ( args )
    order = canonizer.getOrder ( args )
    canonizer.checkAllowedSqrtses ( order, sqrtses )
    inputFiles = canonizer.getInputFiles ( args )
    ncpus = canonizer.checkNCPUs ( args.ncpus, inputFiles )
    pythiaVersion = canonizer.getPythiaVersion ( args )
    ssmultipliers = None
    if hasattr ( args, "ssmultipliers" ):
        ssmultipliers = canonizer.getSSMultipliers ( args.ssmultipliers )
        if ssmultipliers != None:
            for pids,multiplier in ssmultipliers.items():
                if type(pids) != tuple:
                    logger.error ( "keys of ssmultipliers need to be supplied as tuples" )
                    sys.exit()
                if type(multiplier) not in [ int, float ]:
                    logger.error ( "values of ssmultipliers need to be supplied as ints or floats" )
                    sys.exit()

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
            computer = XSecComputer( order, args.nevents, pythiaVersion, \
                                     not args.noautocompile )
            toFile = canonizer.writeToFile ( args )
            computer.computeForBunch (  sqrtses, chunk, not args.keep,
                          args.LOfromSLHA, toFile, pythiacard=pythiacard, \
                        ssmultipliers = ssmultipliers )
            os._exit ( 0 )
        if pid > 0:
            children.append ( pid )
    for child in children:
        r = os.waitpid ( child, 0 )
        logger.debug ( "child %d terminated: %s" % (child,r) )
    logger.debug ( "all children terminated." )
