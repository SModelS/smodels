#!/usr/bin/env python3

"""
.. module:: modelTester
   :synopsis: Functions to test (a set of) points, handling decomposition,
              result and coverage checks, parallelisation.

.. moduleauthor:: Ursula Laa <ursula.laa@lpsc.in2p3.fr>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""


from smodels.base.model import Model
from smodels.base.physicsUnits import GeV, fb, TeV
from smodels.base import runtime
from smodels.decomposition import decomposer
from smodels.experiment.exceptions import DatabaseNotFoundException
from smodels.experiment.databaseObj import Database
from smodels.matching import theoryPrediction
from smodels.matching.theoryPrediction import theoryPredictionsFor,TheoryPredictionsCombiner
from smodels.matching.exceptions import SModelSMatcherError as SModelSError
from smodels.share.models.SMparticles import SMList
from smodels.tools.particlesLoader import load
from smodels.tools import crashReport, timeOut
from smodels.tools.printers.masterPrinter import MPrinter
from smodels.tools.printerTools import printScanSummary
from smodels.base.smodelsLogging import logger
from smodels.tools import ioObjects
from smodels.tools import coverage


from collections import OrderedDict
import os
import sys
import time
import gc
try:
    from ConfigParser import SafeConfigParser, NoSectionError, NoOptionError
except ImportError:
    from configparser import ConfigParser, NoSectionError, NoOptionError
import logging


def testPoint(inputFile, outputDir, parser, database):
    """
    Test model point defined in input file (running decomposition, check
    results, test coverage)

    :parameter inputFile: path to input file
    :parameter outputDir: path to directory where output is be stored
    :parameter parser: ConfigParser storing information from parameters file
    :parameter database: Database holding the list of experiment results

    :return: dictionary with input filename as key and the MasterPrinter object as value
    """


    """ Set BSM model, if necessary """
    if parser.has_option("particles","model"):
        runtime.modelFile = parser.get( "particles", "model" )
    else:
        logger.debug(f'Model file has not been defined. Using input file {inputFile} to read quantum numbers.')
        runtime.modelFile = inputFile

    """Get run parameters and options from the parser"""
    sigmacut = parser.getfloat("parameters", "sigmacut") * fb
    minmassgap = parser.getfloat("parameters", "minmassgap") * GeV
    if parser.has_option("parameters","minmassgapISR"):
        minmassgapISR = parser.getfloat("parameters", "minmassgapISR") * GeV
    else:
        minmassgapISR = 1.0*GeV

    """Setup output printers"""
    masterPrinter = MPrinter()
    masterPrinter.setPrinterOptions(parser)
    masterPrinter.setOutPutFiles(os.path.join(
        outputDir, os.path.basename(inputFile)))

    """ Add list of analyses loaded to printer"""
    masterPrinter.addObj(database)

    """Check input file for errors"""
    inputStatus = ioObjects.FileStatus()
    if parser.getboolean("options", "checkInput"):
        inputStatus.checkFile(inputFile)
    """Initialize output status and exit if there were errors in the input"""
    printParameters = []
    if parser.has_section('parameters'):
        printParameters += list(parser.items('parameters'))
    if parser.has_section('particles'):
        printParameters += list(parser.items('particles'))
    if parser.has_section('options'):
        printParameters += list(parser.items('options'))

    printParameters = OrderedDict(printParameters)
    outputStatus = ioObjects.OutputStatus(inputStatus.status, inputFile,
                                          printParameters,
                                          database.databaseVersion)
    masterPrinter.addObj(outputStatus)
    if outputStatus.status < 0:
        return {os.path.basename(inputFile): masterPrinter}

    """
    Load the input model
    ====================
    """
    try:
        """
        Load the input model and  update it with the information from the input file
        """

        BSMList = load()
        model = Model(BSMparticles=BSMList, SMparticles=SMList)
        promptWidth = None
        stableWidth = None
        ignorePromptQNumbers = []
        if parser.has_option("particles", "promptWidth"):
            promptWidth = parser.getfloat("particles", "promptWidth")*GeV
        if parser.has_option("particles", "stableWidth"):
            stableWidth = parser.getfloat("particles", "stableWidth")*GeV
        if parser.has_option("particles", "ignorePromptQNumbers"):
            ignorePromptQNumbers = parser.get("particles", "ignorePromptQNumbers")
            ignorePromptQNumbers = ignorePromptQNumbers.replace(" ","")
            ignorePromptQNumbers = ignorePromptQNumbers.split(",")
        model.updateParticles(inputFile=inputFile,
                              promptWidth=promptWidth,
                              stableWidth=stableWidth,
                              ignorePromptQNumbers=ignorePromptQNumbers)
    except SModelSError as e:
        logger.error(f"Exception {e} {type(e)}")
        """ Update status to fail, print error message and exit """
        outputStatus.updateStatus(-1)
        return {os.path.basename(inputFile): masterPrinter}

    """
    Decompose input model
    =====================
    """

    try:

        """ Decompose the input model, store the output elements in smstoplist """
        sigmacut = parser.getfloat("parameters", "sigmacut") * fb
        smstoplist = decomposer.decompose(model, sigmacut,
                                          massCompress=parser.getboolean(
                                              "options", "doCompress"),
                                          invisibleCompress=parser.getboolean(
                                              "options", "doInvisible"),
                                          minmassgap=minmassgap,
                                          minmassgapISR=minmassgapISR)
    except SModelSError as e:
        print(f"Exception {e} {type(e)}")
        """ Update status to fail, print error message and exit """
        outputStatus.updateStatus(-1)
        return {os.path.basename(inputFile): masterPrinter}

    """ Print Decomposition output.
        If no topologies with sigma > sigmacut are found, update status, write
        output file, stop running """
    if not smstoplist:
        outputStatus.updateStatus(-3)
        return {os.path.basename(inputFile): masterPrinter}

    masterPrinter.addObj(smstoplist)

    """
    Compute theory predictions
    ====================================================
    """

    """ Get theory prediction for each analysis and print basic output """
    allPredictions = []
    combineResults = False
    useBest = True
    try:
        combineResults = parser.getboolean("options", "combineSRs")
    except (NoSectionError, NoOptionError):
        pass
    try:
        allSRs = parser.getboolean("options", "reportAllSRs")
        if allSRs:  # If set print out all SRs and skip combination
            useBest = False
    except (NoSectionError, NoOptionError):
        pass

    if parser.has_section ( "experimentalFeatures" ):
        featuresDict = dict(parser.items("experimentalFeatures"))
        setExperimentalFeatures( featuresDict )

    allPredictions = theoryPredictionsFor(database, smstoplist,
                                          useBestDataset=useBest,
                                          combinedResults=combineResults)

    """Compute chi-square and likelihood"""
    if parser.getboolean("options", "computeStatistics"):
        for theoPred in allPredictions:
            theoPred.computeStatistics()

    """ Define theory predictions list that collects all theoryPrediction objects which satisfy max condition."""
    maxcond = parser.getfloat("parameters", "maxcond")
    theoryPredictions = theoryPrediction.TheoryPredictionList(
        allPredictions, maxcond)

    if len(theoryPredictions) != 0:
        outputStatus.updateStatus(1)
        theoryPredictions._theoryPredictions = [tp for tp in theoryPredictions._theoryPredictions if not "CR" in os.path.basename(tp.dataset.path)] # Do not print CRs "results"
        masterPrinter.addObj(theoryPredictions)
    else:
        outputStatus.updateStatus(0)  # no results after enforcing maxcond

    if parser.getboolean("options", "testCoverage"):
        """ Testing coverage of model point, add results to the output file """
        if parser.has_option("options", "coverageSqrts"):
            sqrts = parser.getfloat("options", "coverageSqrts")*TeV
        else:
            sqrts = None
        uncovered = coverage.Uncovered(
            smstoplist, sigmacut=sigmacut, sqrts=sqrts)
        masterPrinter.addObj(uncovered)

    if parser.has_option("options", "combineAnas"):
        """ Combine analyses """

        combineAnas = parser.get("options", "combineAnas").replace(" ","").split(",")
        if combineAnas:
            if combineResults is True:
                logger.warning("Combining analyses with signal region combination (combineSRs=True) might significantly reduce CPU performance.")
            combiner = TheoryPredictionsCombiner.selectResultsFrom(theoryPredictions,
                                                                   combineAnas)
            # Only compute combination if at least one result was selected
            if combiner is not None:
                combiner.computeStatistics()
                masterPrinter.addObj(combiner)

    return {os.path.basename(inputFile): masterPrinter}


def runSingleFile(inputFile, outputDir, parser, database,
                  timeout, development, parameterFile):
    """
    Call testPoint on inputFile, write crash report in case of problems

    :parameter inputFile: path to input file
    :parameter outputDir: path to directory where output is be stored
    :parameter parser: ConfigParser storing information from parameter.ini file
    :parameter datbase: Database holding the list of selected results
    :parameter crashReport: if True, write crash report in case of problems
    :param timeout: set a timeout for one model point (0 means no timeout)
    :returns: output of printers
    """

    try:
        with timeOut.Timeout(timeout):
            res = testPoint(inputFile, outputDir, parser, database)
            for fname,mprinter in res.items():
                res[fname] = mprinter.flush()
            return res
    except Exception as e:
        crashReportFacility = crashReport.CrashReport()

        if development:
            print(crashReport.createStackTrace())
            raise e
        else:
            print(crashReport.createStackTrace())
            crashReportFacility.createCrashReportFile(inputFile, parameterFile)
            print(crashReportFacility.createUnknownErrorMessage())
    return {inputFile: None}


def runSetOfFiles(inputFiles, outputDir, parser, database,
                  timeout, development, parameterFile, 
                  return_dict ):
    """
    Loop over all input files in inputFiles with testPoint

    :parameter inputFiles: list of input files to be tested
    :parameter outputDir: path to directory where output is be stored
    :parameter parser: ConfigParser storing information from parameter.ini file
    :parameter database: Database with selected experimental results
    :parameter development: turn on development mode (e.g. no crash report)
    :parameter parameterFile: parameter file, for crash reports
    :returns: printers output
    """

    for inputFile in inputFiles:
        tmp=runSingleFile(inputFile, outputDir, parser, database,
                          timeout, development, parameterFile)
        return_dict.update ( tmp )

def _cleanList(fileList, inDir):
    """ clean up list of files """
    cleanedList = []
    for f in fileList:
        tmp = os.path.join(inDir, f)
        if not os.path.isfile(tmp):
            logger.info(
                f"{tmp} does not exist or is not a file. Skipping it.")
            continue
        cleanedList.append(tmp)
    return cleanedList


def _determineNCPus(cpus_wanted, n_files):
    """ determine the number of CPUs that are to be used.
    :param cpus_wanted: number of CPUs specified in parameter file
    :param n_files: number of files to be run on
    :returns: number of CPUs that are to be used
    """

    ncpusAll = runtime.nCPUs()
    # ncpus = parser.getint("parameters", "ncpus")
    ncpus = cpus_wanted
    if ncpus <= 0:
        ncpus = ncpusAll + ncpus
    ncpus = min(n_files, ncpus)
    if ncpus < 1:
        ncpus = 1
    return ncpus


def testPoints(fileList, inDir, outputDir, parser, database,
               timeout, development, parameterFile):
    """
    Loop over all input files in fileList with testPoint, using ncpus CPUs
    defined in parser

    :param fileList: list of input files to be tested
    :param inDir: path to directory where input files are stored
    :param outputDir: path to directory where output is stored
    :param parser: ConfigParser storing information from parameter.ini file
    :param database: Database with selected experimental results
    :param timeout: set a timeout for one model point (0 means no timeout)
    :param development: turn on development mode (e.g. no crash report)
    :param parameterFile: parameter file, for crash reports
    :returns: printer(s) output, if not run in parallel mode
    """

    t0 = time.time()
    if len(fileList) == 0:
        logger.error("no files given.")
        return None

    cleanedList = _cleanList(fileList, inDir)
    ncpus = _determineNCPus(parser.getint(
        "parameters", "ncpus"), len(cleanedList))
    nFiles = len(cleanedList)

    if nFiles == 0:
        logger.error("No valid input files found")
        return None
    elif nFiles == 1:
        logger.info("Running SModelS for a single file")
        runSingleFile(cleanedList[0], outputDir, parser,
                      database, timeout,
                      development, parameterFile)
    else:

        if ncpus == 1:
            logger.info("Running SModelS for %i files with a single process. Messages will be redirected to smodels.log"
                        % (nFiles))

            for hdlr in logger.handlers[:]:
                logger.removeHandler(hdlr)
                hdlr.close()
            fileLog = logging.FileHandler('./smodels.log')
            logger.addHandler(fileLog)

            # Run a single process:
            outputDict = {}
            runSetOfFiles(cleanedList, outputDir, parser,
                            database, timeout,
                            development, parameterFile, 
                            outputDict )
        else:
            logger.info("Running SModelS for %i files with %i processes. Messages will be redirected to smodels.log"
                        % (nFiles, ncpus))

            for hdlr in logger.handlers[:]:
                logger.removeHandler(hdlr)
                hdlr.close()
            fileLog = logging.FileHandler('./smodels.log')
            logger.addHandler(fileLog)

            # Launch multiple processes.
            # Split list of files
            chunkedFiles = [cleanedList[x::ncpus] for x in range(ncpus)]
            children = []
            from multiprocessing import Process, Manager
            manager = Manager()
            outputDict = manager.dict()
            for chunkFile in chunkedFiles:
                args = ( chunkFile, outputDir, parser, database, timeout, 
                         development, parameterFile, outputDict )
                p = Process ( target=runSetOfFiles, args = args )
                p.start()
                              
                children.append(p)
            ctr = 0
            nsteps = 10
            for p in children:
                p.join()
                ctr+=1
                if ctr % nsteps == 10:
                    t=(time.time()-t0)/60.
                    logger.info ( f"{ctr} of {len(children)} processes done in {t:.2f} min" )
            """
            iprint, nprint = 5, 5  # Define when to start printing and the percentage step
            # Check process progress until they are all finished
            while True:
                done = sum([p.ready() for p in children])
                fracDone = 100*float(done)/len(children)
                if fracDone >= iprint:
                    while fracDone >= iprint:
                        iprint += nprint
                    logger.info('%i%% of processes done in %1.2f min' %
                                (iprint-nprint, (time.time()-t0)/60.))
                if done == len(children):
                    break
                time.sleep(2)

            logger.debug("All children terminated")
            """

        # Collect output to build global summary:
        scanSummaryFile = os.path.join(outputDir, 'summary.txt')
        logger.info(f"A summary of the scan results can be found in {scanSummaryFile}")
        printScanSummary(outputDict, scanSummaryFile)
        # Remove summary log from logger
        logger.removeHandler(fileLog)
        fileLog.close()

    logger.info(f"Done in {(time.time() - t0) / 60.0:3.2f} min")

    return None


def checkForSemicolon(strng, section, var):
    if ";" in strng:
        logger.warning(
            f"A semicolon(;) has been found in [{section}] {var}, in your config file. If this was meant as comment, then please a space before it.")


def loadDatabase(parser, db):
    """
    Load database

    :parameter parser: ConfigParser with path to database
    :parameter db: binary database object. If None, then database is loaded,
                   according to databasePath. If True, then database is loaded,
                   and text mode is forced.
    :returns: database object

    """
    try:
        dp = parser.get("path", "databasePath")
        logger.error("``[path] databasePath'' in ini file is deprecated; "
                     "use ``[database] path'' instead.(See e.g. smodels/etc/parameters_default.ini)")
        parser.set("database", "path", dp)
    except (NoSectionError, NoOptionError):
        # path.databasePath not set. This is good.
        pass
    try:
        database = db
        # logger.error("database=db: %s" % database)
        if database in [None, True]:
            databasePath = parser.get("database", "path")
            checkForSemicolon(databasePath, "database", "path")
            force_load = None
            if database is True:
                force_load = "txt"
            if os.path.isfile(databasePath):
                force_load = "pcl"
            database = Database(databasePath, force_load=force_load)
    except DatabaseNotFoundException:
        logger.error(f"Database not found in ``{os.path.realpath(databasePath)}''")
        sys.exit()
    return database


def loadDatabaseResults(parser, database):
    """
    Restrict the (active) database results to the ones specified in parser

    :parameter parser: ConfigParser, containing analysis and txnames selection
    :parameter database: Database object

    """
    """ In case that a list of analyses or txnames are given, retrieve list """
    tmp = parser.get("database", "analyses").split(",")
    analyses = [x.strip() for x in tmp]
    tmp_tx = parser.get("database", "txnames").split(",")
    txnames = [x.strip() for x in tmp_tx]
    if parser.get("database", "dataselector") == "efficiencyMap":
        dataTypes = ['efficiencyMap']
        datasetIDs = ['all']
    elif parser.get("database", "dataselector") == "upperLimit":
        dataTypes = ['upperLimit']
        datasetIDs = ['all']
    else:
        dataTypes = ['all']
        tmp_dIDs = parser.get("database", "dataselector").split(",")
        datasetIDs = [x.strip() for x in tmp_dIDs]

    useNonValidated = False
    if parser.has_option("database", "useNonValidated"):
        useNonValidated = parser.getboolean("database", "usenonvalidated")
    if useNonValidated:
        logger.info('Including non-validated results')

    """ Load analyses """

    database.selectExpResults(analysisIDs=analyses, txnames=txnames,
                                 datasetIDs=datasetIDs, dataTypes=dataTypes,
                                 useNonValidated=useNonValidated)


def getParameters(parameterFile):
    """
    Read parameter file, exit in case of errors

    :parameter parameterFile: Path to parameter File
    :returns: ConfigParser read from parameterFile

    """
    try:
        parser = ConfigParser(inline_comment_prefixes=(';',))
    except NameError:
        parser = SafeConfigParser()
    ret = parser.read(parameterFile)
    if ret == []:
        logger.error(f"No such file or directory: '{parameterFile}'")
        sys.exit()

    if parser.has_section ( "experimentalFeatures" ):
        featuresDict = dict(parser.items("experimentalFeatures"))
        setExperimentalFeatures(featuresDict)
    
    try:
        runtime.modelFile = parser.get("particles", "model")
    except:
        pass
    try:
        pyhfbackend = parser.get("options","pyhfbackend")
        from smodels.statistics import pyhfInterface
        r = pyhfInterface.setBackend ( pyhfbackend )
    except:
        pass
    return parser


def setExperimentalFeatures(featuresDict):
    """ set the experimental features flats, if experimentalFeatures:* = True """

    for feature in featuresDict.keys():
        if not feature in runtime._experimental:
            logger.warning ( f"'{feature}' is not a known experimental feature. will ignore." )
            continue
        flag = False
        if featuresDict[feature].lower() in [ "true", "1", "yes" ]:
           flag = True
        runtime._experimental[feature]=flag

def getAllInputFiles(inFile):
    """
    Given inFile, return list of all input files

    :parameter inFile: Path to input file or directory containing input files
    :returns: List of all input files, and the directory name

    """
    if os.path.isdir(inFile):
        fileList = os.listdir(inFile)
        return fileList, inFile
    fileList = [os.path.basename(inFile)]
    return fileList, os.path.dirname(inFile)


def getCombiner(inputFile,parameterFile):
    """
    Facility for running SModelS, computing the theory predictions and returning the combination of analyses
    (defined in the parameterFile). Useful for plotting likelihoods!.
    Extracts and returns the TheoryPredictionsCombiner object from the master printer, if the object is found. Return None otherwise.

    :param inputFile: path to the input SLHA file
    :param parameterFile: path to parameters.ini file

    :return: TheoryPredictionsCombiner object generated by running SModelS.
    """

    # Get parameters
    parser = getParameters(parameterFile)

    # Load database and results
    database = loadDatabase(parser, None)
    loadDatabaseResults(parser, database)

    # Run SModelS for a single file and get the printer
    outputDir = './'
    output = testPoint(inputFile, outputDir, parser,
                       database)
    mprinter = list(output.values())[0]
    # Try to exctract the TheoryPredictionsCombiner object from one of the printers.
    combiner = None
    for p in mprinter.Printers.values():
        if combiner is not None:
            break
        for obj in p.toPrint:
            if isinstance(obj,TheoryPredictionsCombiner):
                combiner = obj
                break
    if combiner is None:
        logger.info("Combiner not found for input file %s with parameters from %s. Is combineAnas defined correctly? (At least one printer must be defined).")


    return combiner
