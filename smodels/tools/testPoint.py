#!/usr/bin/env python

from smodels.tools import ioObjects
from smodels.tools import coverage, runtime
from smodels.theory import slhaDecomposer
from smodels.theory import lheDecomposer
from smodels.theory.theoryPrediction import theoryPredictionsFor
import smodels.tools.printer as prt
import logging
import os
import sys
from smodels.tools.physicsUnits import GeV, fb

log = logging.getLogger(__name__)

def testPoint(inputFile, outputDir, parser, databaseVersion, listOfExpRes):
    """ Minimum value of cross-section for an element to be considered eligible
        for decomposition.  Too small sigmacut leads to too large decomposition
        time.  """
    sigmacut = parser.getfloat("parameters", "sigmacut") * fb
    """ Minimum value for considering two states non-degenerate (only used for
        mass compression) """
    minmassgap = parser.getfloat("parameters", "minmassgap") * GeV
    inputType = parser.get("options", "inputType").lower()


    outputFile = os.path.join(outputDir, os.path.basename(inputFile))+'.smodels'

    if os.path.exists(outputFile):
        log.warning("Removing old output file " + outputFile)
    outfile = open(outputFile, 'w')
    outfile.close()

    """ Check input file for errors """
    inputStatus = ioObjects.FileStatus()
    if parser.getboolean("options", "checkInput"):
        inputStatus.checkFile(inputType, inputFile, sigmacut)
    """ Initialize output status and exit if there were errors in the input """
    outputStatus = ioObjects.OutputStatus(inputStatus.status, inputFile,
            dict(parser.items("parameters")), databaseVersion, outputFile)
    if outputStatus.status < 0: return

    """ Setup output printers """
    stdoutPrinter = prt.TxTPrinter(output = 'stdout')
    summaryPrinter = prt.SummaryPrinter(output = 'file', filename = outputFile)

    """
    Decompose input file
    ====================
    """
    try:
        """ Decompose input SLHA file, store the output elements in smstoplist """
        if inputType == 'slha':
            smstoplist = slhaDecomposer.decompose(inputFile, sigmacut,
                    doCompress=parser.getboolean("options", "doCompress"),
                    doInvisible=parser.getboolean("options", "doInvisible"),
                    minmassgap=minmassgap)
        else:
            smstoplist = lheDecomposer.decompose(inputFile,
                    doCompress=parser.getboolean("options", "doCompress"),
                    doInvisible=parser.getboolean("options", "doInvisible"),
                    minmassgap=minmassgap)
    except:
        """ Update status to fail, print error message and exit """
        outputStatus.updateStatus(-1)
        return

    """ Print Decomposition output.
        If no topologies with sigma > sigmacut are found, update status, write
        output file, stop running """
    if not smstoplist:
        outputStatus.updateStatus(-3)
        return

    outLevel = 0
    if parser.getboolean("stdout", "printDecomp"):
        outLevel = 1
        outLevel += parser.getboolean("stdout", "addElmentInfo")
    stdoutPrinter.addObj(smstoplist,outLevel)

    """
    Compute theory predictions
    ====================================================
    """

    """ Get theory prediction for each analysis and print basic output """
    allPredictions = []
    for expResult in listOfExpRes:
        theorypredictions = theoryPredictionsFor(expResult, smstoplist)
        if not theorypredictions: continue
        if parser.getboolean("stdout", "printResults"):
            stdoutPrinter.addObj(theorypredictions)

        allPredictions += theorypredictions._theoryPredictions

    """ Define result list that collects all theoryPrediction objects."""
    maxcond = parser.getfloat("parameters", "maxcond")
    results = ioObjects.ResultList(allPredictions,maxcond)
    if not parser.getboolean("file", "expandedSummary"):
        results.useBestResult()

    outLevel = 0
    if not results.isEmpty():
        outputStatus.updateStatus(1)
        outLevel = 1
        outLevel += parser.getboolean("file", "addConstraintInfo")
    else:
        outputStatus.updateStatus(0) # no results after enforcing maxcond
    stdoutPrinter.addObj(outputStatus)
    summaryPrinter.addObj(outputStatus)
    summaryPrinter.addObj(results,outLevel)
    if parser.getboolean("stdout", "printResults"):
        stdoutPrinter.addObj(results,outLevel)

    if parser.getboolean("options", "testCoverage"):
        """ Look for missing topologies, add them to the output file """
        uncovered = coverage.Uncovered(smstoplist)
        summaryPrinter.addObj(uncovered.missingTopos)
        stdoutPrinter.addObj(uncovered.missingTopos,2)
        summaryPrinter.addObj(uncovered,2)
        stdoutPrinter.addObj(uncovered,2) 
    stdoutPrinter.close()
    summaryPrinter.close()

def runSingleFile ( inputFile, outputDir, parser, databaseVersion, listOfExpRes,
                    crashReport=True):
    try:
        testPoint( inputFile, outputDir, parser, databaseVersion,
                             listOfExpRes )
    except Exception,e:
        if crashReport:
            e.inputFile = inputFile
            raise e

def runSetOfFiles ( inputFiles, outputDir, parser, databaseVersion, listOfExpRes,
                    crashReport=True ):
    for inputFile in inputFiles:
        runSingleFile( inputFile, outputDir, parser, databaseVersion,
                       listOfExpRes )

def runAllFiles ( fileList, inDir, outputDir, parser, databaseVersion, listOfExpRes ):
    """ run over all files """
    if len( fileList ) == 0:
        log.error ( "no files given." )
        return
    if len(fileList ) == 1:
        runSingleFile ( fileList[0], outputDir, parser, databaseVersion, listOfExpRes )
        return

    """ loop over input files and run SModelS """
    ncpusAll = runtime.nCPUs()
    ncpus = parser.getint("parameters", "ncpus")
    if ncpus < 1 or ncpus > ncpusAll: ncpus = ncpusAll
    log.info ("we run on %d cores" % ncpus )

    cleanedList = []
    for f in fileList:
        tmp = os.path.join(inDir, f )
        if not os.path.isfile ( tmp ):
            log.info ( "%s does not exist or is not a file. Skipping it." % tmp )
            continue
        cleanedList.append ( tmp )

    runOnSingleCore = False
    if runOnSingleCore:
        runSetOfFiles ( cleanedList, outputDir, parser, databaseVersion,
                        listOfExpRes )
        return

    ### now split up for every fork
    chunkedFiles = [ cleanedList[x::ncpus] for x in range(ncpus ) ]
    children = []
    for (i,chunk) in enumerate ( chunkedFiles ):
        pid=os.fork()
        print "Forking: ",i,pid,os.getpid()
        if pid == 0:
            log.info ("chunk #%d: pid %d (parent %d)." %
                    ( i, os.getpid(), os.getppid() ) )
            log.info ( " `-> %s" % " ".join ( chunk ) )
            runSetOfFiles ( chunk, outputDir, parser, databaseVersion, listOfExpRes )
            sys.exit()
            # return
            # continue
        if pid < 0:
            log.error ( "fork did not succeed! Pid=%d" % pid )
            sys.exit()
        if pid > 0:
            children.append ( pid )
    for child in children:
        r = os.waitpid ( child, 0 )
        log.info ( "child %d terminated: %s" % (child,r) )
    log.info ( "all children terminated" )

