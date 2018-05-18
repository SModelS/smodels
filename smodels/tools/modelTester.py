#!/usr/bin/env python

"""
.. module:: tools.modelTester
   :synopsis: Functions to test (a set of) points, handling decomposition,
              result and coverage checks, parallelisation.

.. moduleauthor:: Ursula Laa <ursula.laa@lpsc.in2p3.fr>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from smodels.tools import ioObjects
from smodels.tools import coverage, runtime
from smodels.theory import slhaDecomposer
from smodels.theory import lheDecomposer
from smodels.theory.theoryPrediction import theoryPredictionsFor, TheoryPredictionList
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools import crashReport, timeOut
from smodels.tools.printer import MPrinter
import os
import sys
import time
try:
    from ConfigParser import SafeConfigParser
except ImportError as e:
    from configparser import ConfigParser
from smodels.tools.physicsUnits import GeV, fb, TeV
from smodels.experiment.exceptions import DatabaseNotFoundException
from smodels.experiment.databaseObj import Database, ExpResultList
from smodels.tools.smodelsLogging import logger

def testPoint(inputFile, outputDir, parser, databaseVersion, listOfExpRes):
    """
    Test model point defined in input file (running decomposition, check
    results, test coverage)

    :parameter inputFile: path to input file
    :parameter outputDir: path to directory where output is be stored
    :parameter parser: ConfigParser storing information from parameters file
    :parameter databaseVersion: Database version (printed to output file)
    :parameter listOfExpRes: list of ExpResult objects to be considered
    :returns: output of printers
    """

    """Get run parameters and options from the parser"""
    sigmacut = parser.getfloat("parameters", "sigmacut") * fb
    minmassgap = parser.getfloat("parameters", "minmassgap") * GeV
    inputType = runtime.filetype ( inputFile )


    """Setup output printers"""
    masterPrinter = MPrinter()
    masterPrinter.setPrinterOptions(parser)
    masterPrinter.setOutPutFiles(os.path.join(outputDir, os.path.basename(inputFile)))

    """ Add list of analyses loaded to printer"""
    masterPrinter.addObj(ExpResultList(listOfExpRes))

    """Check input file for errors"""
    inputStatus = ioObjects.FileStatus()
    if parser.getboolean("options", "checkInput"):
        inputStatus.checkFile( inputFile, sigmacut)
    """Initialize output status and exit if there were errors in the input"""
    outputStatus = ioObjects.OutputStatus(inputStatus.status, inputFile,
            dict(parser.items("parameters")), databaseVersion)
    masterPrinter.addObj(outputStatus)
    if outputStatus.status < 0:
        return masterPrinter.flush()


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
    except SModelSError as e:
        print ( "Exception %s %s" % ( e, type(e) ) )
        """ Update status to fail, print error message and exit """
        outputStatus.updateStatus(-1)
        return masterPrinter.flush()

    """ Print Decomposition output.
        If no topologies with sigma > sigmacut are found, update status, write
        output file, stop running """
    if not smstoplist:
        outputStatus.updateStatus(-3)
        return masterPrinter.flush()

    masterPrinter.addObj(smstoplist)


    """
    Compute theory predictions
    ====================================================
    """

    """ Get theory prediction for each analysis and print basic output """
    allPredictions = []
    for expResult in listOfExpRes:
        # logger.error ("FIXME maybe need to set useBestDataset to FALSE")
        theorypredictions = theoryPredictionsFor(expResult, smstoplist,useBestDataset=True )
        if not theorypredictions: continue
        allPredictions += theorypredictions._theoryPredictions

    """Compute chi-square and likelihood"""
    if parser.getboolean("options","computeStatistics"):
        for theoPred in allPredictions:
            theoPred.computeStatistics()

    """ Define result list that collects all theoryPrediction objects."""
    maxcond = parser.getfloat("parameters", "maxcond")
    results = ioObjects.ResultList(allPredictions,maxcond)

    if not results.isEmpty():
        outputStatus.updateStatus(1)
        masterPrinter.addObj(results)
    else:
        outputStatus.updateStatus(0) # no results after enforcing maxcond

    if parser.getboolean("options", "testCoverage"):
        """ Testing coverage of model point, add results to the output file """
        if  parser.has_option("options","coverageSqrts"):
            sqrts = parser.getfloat("options", "coverageSqrts")*TeV
        else:
            sqrts = None
        uncovered = coverage.Uncovered(smstoplist,sqrts=sqrts)
        masterPrinter.addObj(uncovered)

    return masterPrinter.flush()

def runSingleFile(inputFile, outputDir, parser, databaseVersion, listOfExpRes,
                    timeout, development, parameterFile):
    """
    Call testPoint on inputFile, write crash report in case of problems

    :parameter inputFile: path to input file
    :parameter outputDir: path to directory where output is be stored
    :parameter parser: ConfigParser storing information from parameter.ini file
    :parameter databaseVersion: Database version (printed to output file)
    :parameter listOfExpRes: list of ExpResult objects to be considered
    :parameter crashReport: if True, write crash report in case of problems
    :param timeout: set a timeout for one model point (0 means no timeout)
    :returns: output of printers
    """
    try:
        with timeOut.Timeout(timeout):
            return testPoint(inputFile, outputDir, parser, databaseVersion,
                             listOfExpRes)
    except Exception as e:
        crashReportFacility = crashReport.CrashReport()

        if development:
            print(crashReport.createStackTrace())
            raise e
        else:
            print(crashReport.createStackTrace())
            crashReportFacility.createCrashReportFile( inputFile, parameterFile )
            print(crashReportFacility.createUnknownErrorMessage())
    return None

def runSetOfFiles(inputFiles, outputDir, parser, databaseVersion, listOfExpRes,
                    timeout, development, parameterFile):
    """
    Loop over all input files in inputFiles with testPoint

    :parameter inputFiles: list of input files to be tested
    :parameter outputDir: path to directory where output is be stored
    :parameter parser: ConfigParser storing information from parameter.ini file
    :parameter databaseVersion: Database version (printed to output file)
    :parameter listOfExpRes: list of ExpResult objects to be considered
    :parameter development: turn on development mode (e.g. no crash report)
    :parameter parameterFile: parameter file, for crash reports
    :returns: printers output
    """
    a={}
    n=len(inputFiles)
    for i,inputFile in enumerate(inputFiles):
        txt=""
        if n>100:
            txt="[%d/%d] " % ( i, n )
        logger.info ( "Start testing %s%s" % (txt, os.path.relpath ( inputFile ) ) )
        a[inputFile] = runSingleFile(inputFile, outputDir, parser, databaseVersion,
                                  listOfExpRes, timeout, development, parameterFile)
    return a

def _cleanList ( fileList, inDir ):
    """ clean up list of files """
    cleanedList = []
    for f in fileList:
        tmp = os.path.join(inDir, f )
        if not os.path.isfile ( tmp ):
            logger.info ( "%s does not exist or is not a file. Skipping it." % tmp )
            continue
        cleanedList.append( tmp )
    return cleanedList

def _determineNCPus ( cpus_wanted, n_files ):
    """ determine the number of CPUs that are to be used.
    :param cpus_wanted: number of CPUs specified in parameter file
    :param n_files: number of files to be run on
    :returns: number of CPUs that are to be used
    """
    ncpusAll = runtime.nCPUs()
    # ncpus = parser.getint("parameters", "ncpus")
    ncpus = cpus_wanted
    if ncpus == 0 or ncpus < -1:
        logger.error ( "Weird number of ncpus given in ini file: %d" % ncpus )
        sys.exit()
    if ncpus == -1 or ncpus > ncpusAll: ncpus = ncpusAll
    ncpus = min ( n_files, ncpus )
    return ncpus

def testPoints(fileList, inDir, outputDir, parser, databaseVersion,
                 listOfExpRes, timeout, development, parameterFile):
    """
    Loop over all input files in fileList with testPoint, using ncpus CPUs
    defined in parser

    :param fileList: list of input files to be tested
    :param inDir: path to directory where input files are stored
    :param outputDir: path to directory where output is stored
    :param parser: ConfigParser storing information from parameter.ini file
    :param databaseVersion: Database version (printed to output files)
    :param listOfExpRes: list of ExpResult objects to be considered
    :param timeout: set a timeout for one model point (0 means no timeout)
    :param development: turn on development mode (e.g. no crash report)
    :param parameterFile: parameter file, for crash reports
    :returns: printer(s) output, if not run in parallel mode
    """
    if len( fileList ) == 0:
        logger.error ( "no files given." )
        return None

    cleanedList = _cleanList ( fileList, inDir )
    if len(cleanedList) == 1:
        return runSingleFile ( cleanedList[0], outputDir, parser, databaseVersion,
                               listOfExpRes, timeout, development, parameterFile )
    ncpus = _determineNCPus ( parser.getint("parameters", "ncpus"), len(cleanedList) )
    logger.info ("Running SModelS on %d cores" % ncpus )

    if ncpus == 1:
        return runSetOfFiles( cleanedList, outputDir, parser, databaseVersion,
                              listOfExpRes, timeout, development, parameterFile )

    ### now split up for every fork
    chunkedFiles = [cleanedList[x::ncpus] for x in range(ncpus)]
    children = []
    for (i,chunk) in enumerate(chunkedFiles):
        pid=os.fork()
        logger.debug("Forking: %s %s %s " % ( i,pid,os.getpid() ) )
        if pid == 0:
            logger.debug("chunk #%d: pid %d (parent %d)." %
                    ( i, os.getpid(), os.getppid() ) )
            logger.debug( " `-> %s" % " ".join ( chunk ) )
            runSetOfFiles(chunk, outputDir, parser, databaseVersion,
                            listOfExpRes, timeout, development, parameterFile)
            os._exit(0) ## not sys.exit(), return, nor continue
        if pid < 0:
            logger.error ( "fork did not succeed! Pid=%d" % pid )
            sys.exit()
        if pid > 0:
            children.append ( pid )
    for child in children:
        r = os.waitpid ( child, 0 )
        logger.debug ( "child %d terminated: %s" % (child,r) )
    logger.debug ( "all children terminated" )
    logger.debug ( "returning no output, because we are in parallel mode" )
    return None

def checkForSemicolon ( strng, section, var ):
    if ";" in strng:
        logger.warning ( "A semicolon (;) has been found in [%s] %s, in your config file. If this was meant as comment, then please a space before it." % ( section, var) )

def loadDatabase(parser, db):
    """
    Load database

    :parameter parser: ConfigParser with path to database
    :parameter db: binary database object. If None, then database is loaded,
                   according to databasePath. If True, then database is loaded,
                   and text mode is forced.
    :returns: database object, database version

    """
    try:
        dp = parser.get ( "path", "databasePath" )
        logger.error ( "``[path] databasePath'' in ini file is deprecated; " \
           "use ``[database] path'' instead. (See e.g. smodels/etc/parameters_default.ini)" )
        parser.set ( "database", "path", dp )
    except Exception as e:
        ## path.databasePath not set. This is good.
        pass
    try:
        database = db
        # logger.error ( "database=db: %s" % database )
        if database in [ None, True ]:
            databasePath = parser.get( "database", "path" )
            checkForSemicolon ( databasePath, "database", "path" )
            discard_zeroes = True
            try:
                discard_zeroes = parser.getboolean( "database", "discardZeroes" )
            except Exception as e: ## too complicated to be more specific
                logger.info ( "database:discardZeroes is not given in config file. Defaulting to 'True'." )
            force_load=None
            if database == True: force_load="txt"
            if os.path.isfile ( databasePath ):
                force_load="pcl"
            database = Database( databasePath, force_load=force_load, \
                                 discard_zeroes = discard_zeroes )
        databaseVersion = database.databaseVersion
    except DatabaseNotFoundException:
        logger.error("Database not found in ``%s''" % os.path.realpath(databasePath))
        sys.exit()
    return database, databaseVersion

def loadDatabaseResults(parser, database):
    """
    Load database entries specified in parser

    :parameter parser: ConfigParser, containing analysis and txnames selection
    :parameter database: Database object
    :returns: List of experimental results

    """
    """ In case that a list of analyses or txnames are given, retrieve list """
    tmp = parser.get("database", "analyses").split(",")
    analyses = [ x.strip() for x in tmp ]
    tmp_tx = parser.get("database", "txnames").split(",")
    txnames = [ x.strip() for x in tmp_tx ]
    if parser.get("database", "dataselector") == "efficiencyMap":
        dataTypes = ['efficiencyMap']
        datasetIDs = ['all']
    elif parser.get("database", "dataselector") == "upperLimit":
        dataTypes = ['upperLimit']
        datasetIDs = ['all']
    else:
        dataTypes = ['all']
        tmp_dIDs = parser.get("database", "dataselector").split(",")
        datasetIDs = [ x.strip() for x in tmp_dIDs ]

    useSuperseded=False
    useNonValidated=False
    if parser.has_option("database","useSuperseded"):
        useSuperseded = parser.getboolean("database", "usesuperseded")
    if parser.has_option("database","useNonValidated"):
        useNonValidated = parser.getboolean("database", "usenonvalidated")
    if useSuperseded:
        logger.info('Including superseded results')
    if useNonValidated:
        logger.info('Including non-validated results')


    """ Load analyses """

    ret = database.getExpResults(analysisIDs=analyses, txnames=txnames,
                                 datasetIDs=datasetIDs, dataTypes=dataTypes,
                                 useSuperseded=useSuperseded, useNonValidated=useNonValidated)
    return ret

def getParameters(parameterFile):
    """
    Read parameter file, exit in case of errors

    :parameter parameterFile: Path to parameter File
    :returns: ConfigParser read from parameterFile

    """
    try:
        parser = ConfigParser( inline_comment_prefixes=( ';', ) )
    except Exception as e:
        parser = SafeConfigParser()
    ret=parser.read(parameterFile)
    if ret == []:
        logger.error ( "No such file or directory: '%s'" % parameterFile )
        sys.exit()
    return parser

def getAllInputFiles(inFile):
    """
    Given inFile, return list of all input files

    :parameter inFile: Path to input file or directory containing input files
    :returns: List of all input files

    """
    if os.path.isdir(inFile):
        fileList = os.listdir(inFile)
        return fileList, inFile
    fileList = [ os.path.basename ( inFile ) ]
    return fileList, os.path.dirname ( inFile )
