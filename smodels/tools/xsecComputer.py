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
from smodels.tools import toolBox
from smodels.base.physicsUnits import pb, TeV, GeV
from smodels.base import crossSection, runtime
from smodels.base.crossSection import LO, NLO, NLL
from smodels.base.smodelsLogging import logger, setLogLevel
from smodels.decomposition.exceptions import SModelSDecompositionError as SModelSError
from smodels.tools.xsecBase import XSecBase, ArgsStandardizer
import os, copy
import pyslha
try:
    import cStringIO as io
except ImportError as e:
    import io
import sys

class XSecComputer(XSecBase):
    """ cross section computer class, what else? """
    def __init__ ( self, maxOrder, nevents, pythiaVersion, maycompile=True,
                   defaulttempdir : str = "/tmp/" ):
        """
        :param maxOrder: maximum order to compute the cross section, given as an integer
                    if maxOrder == LO, compute only LO pythia xsecs
                    if maxOrder == NLO, apply NLO K-factors from NLLfast (if available)
                    if maxOrder == NLL, apply NLO+NLL K-factors from NLLfast (if available)
        :param nevents: number of events for pythia run
        :param pythiaVersion: pythia6 or pythia8 (integer)
        :param maycompile: if True, then tools can get compiled on-the-fly
        :param defaulttempdir: the default temp directory
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
            logger.error ( f"Unknown pythia version {pythiaVersion}. Allowed values: 6, 8" )
            sys.exit()
        self.pythiaVersion = pythiaVersion
        self.defaulttempdir = defaulttempdir

    def _checkSLHA ( self, slhafile ):
        if not os.path.isfile(slhafile):
            logger.error( f"SLHA file {slhafile} not found." )
            raise SModelSError()
        try:
            f=pyslha.readSLHAFile(slhafile)
        except (pyslha.ParseError,ValueError) as e:
            logger.error( f"File {slhafile} cannot be parsed as SLHA file: {e}" )
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
        s = float(sqrts/TeV)
        if abs (s % 1) < 1e-5:
            s = int(s)
        wlabel = str(s) + ' TeV'
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
            nllfastnr = sqrts.asNumber(TeV)
            if nllfastnr % 1 == 0:
                nllfastnr = int(nllfastnr)
                nllfast = toolBox.ToolBox().get( f"nllfast{nllfastnr}" )
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
            else:
                logger.error ( f"nllfast{nllfastnr} not yet available! will compute LO only!" )

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
                    logger.error ( f"I have a pid of type {type(spid)}. Dont know what to do." % \
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
                logger.warning ( f"currently we always only have two mothers, so why are the signal strength multipliers given for {len(pids)} mothers?" )
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
        ret.defaulttempdir = self.defaulttempdir
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
        logger.debug ( f"how many NLL xsecs? {len(self.xsecs)}" )
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
            logger.info(f"Computing SLHA cross section from {inputFile}." )
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
            logger.debug ( f"computing xsec for {inputFile}" )
            self.computeForOneFile ( sqrtses, inputFile, unlink, lOfromSLHA,
                      tofile, pythiacard=pythiacard, ssmultipliers = ssmultipliers )

    def addCommentToFile ( self, comment, slhaFile ):
        """ add the optional comment to file """
        if comment in [ None, "" ]:
            return
        if not os.path.isfile(slhaFile ):
            logger.error(f"SLHA file {slhaFile} not found." )
            raise SModelSError()
        outfile = open(slhaFile, 'a')
        outfile.write ( f"# {comment}\n" )
        outfile.close()

    def addMultipliersToFile ( self, ssmultipliers, slhaFile ):
        """ add the signal strength multipliers to the SLHA file """
        if ssmultipliers in [ None, {} ]:
            return
        if not os.path.isfile(slhaFile ):
            logger.error(f"SLHA file {slhaFile} not found." )
            raise SModelSError()
        tokens = []
        for k,v in ssmultipliers.items():
            tokens.append ( f"{k}:{v:.4g}" )
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

def main(args):
    canonizer = ArgsStandardizer()
    setLogLevel ( args.verbosity )
    if not hasattr ( args, "noautocompile" ):
        args.noautocompile = False
    if args.query:
        return canonizer.queryCrossSections ( args.filename )
    if args.colors:
        from smodels.base.smodelsLogging import colors
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
            logger.debug ( f"chunk #{i}: pid {os.getpid()} (parent {os.getppid()})." )
            logger.debug ( " `-> {' '.join(chunk)}" )
            computer = XSecComputer( order, args.nevents, pythiaVersion, \
                                   not args.noautocompile, canonizer.tempDir(args) )
            toFile = canonizer.writeToFile ( args )
            computer.computeForBunch (  sqrtses, chunk, not args.keep,
                          args.LOfromSLHA, toFile, pythiacard=pythiacard, \
                        ssmultipliers = ssmultipliers )
            os._exit ( 0 )
        if pid > 0:
            children.append ( pid )
    for child in children:
        r = os.waitpid ( child, 0 )
        logger.debug ( f"child {child} terminated: {r}" )
    logger.debug ( "all children terminated." )
