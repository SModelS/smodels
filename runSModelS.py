#!/usr/bin/env python

"""
.. module:: runSModelS
   :synopsis: Main code for running SModelS.
   
"""

from __future__ import print_function
import os, sys
import logging
from ConfigParser import SafeConfigParser
from smodels.experiment.databaseObj import Database
from smodels.installation import installDirectory
from smodels.theory import slhaDecomposer
from smodels.theory import lheDecomposer
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.tools.physicsUnits import GeV, fb
from smodels.tools import ioObjects
from smodels.tools import coverage
from smodels.tools import crashReport, timeOut
import smodels.tools.printer as prt
from smodels.experiment.exceptions import DatabaseNotFoundException

log = logging.getLogger(__name__)
currentFile = ""


def runSingleFile ( inFile, inputFile, outputDir ):
    return

def main(inFile, parameterFile, outputDir, verbosity = 'info', db=None ):
    """
    Provides a command line interface to basic SModelS functionalities.
    
    :param inputFile: input file name (either a SLHA or LHE file)
                      or directory name (path to directory containing input files)
    :param parameterFile: File containing the input parameters (default =
                          /etc/parameters_default.ini)
    :param outputDir: Output directory to write a summary of results to
    :param db: supply a smodels.experiment.databaseObj.Database object, so
            the database doesn't have to be loaded anymore. Will
            render a few parameters in the parameter file irrelevant.
            If None, load the database as described in parameterFile,
            If True, force loading the text database.
    
    """

    """
    Read and check parameter file
    =========================
    """
    global currentFile
    parser = SafeConfigParser()
    ret=parser.read(parameterFile)
    if ret == []:
        log.error ( "No such file or directory: '%s'" % parameterFile )
        sys.exit()

    """ Minimum value of cross-section for an element to be considered eligible
        for decomposition.  Too small sigmacut leads to too large decomposition
        time.  """
    sigmacut = parser.getfloat("parameters", "sigmacut") * fb

    """ Minimum value for considering two states non-degenerate (only used for
        mass compression) """
    minmassgap = parser.getfloat("parameters", "minmassgap") * GeV

    inputType = parser.get("options", "inputType").lower()
    if inputType != 'slha' and inputType != 'lhe':
        log.error("Unknown input type (must be SLHA or LHE): %s" % inputType)
        return

    """ Check database location and load database"""
    try:
        databasePath = parser.get("path", "databasePath")
        database = db
        if database in [ None, True ]:
            force_load=None
            if database == True: force_load="txt"
            database = Database(databasePath, force_load=force_load, verbosity=verbosity )
        databaseVersion = database.databaseVersion
    except DatabaseNotFoundException:
        log.error("Database not found in %s" % os.path.realpath(databasePath))
        databaseVersion = None
        return

    """ Get list of input files to be tested """
    if os.path.isdir(inFile):
        fileList = os.listdir(inFile)
    else: fileList = [inFile]

    """ Create output directory if missing """
    if not os.path.isdir(outputDir): os.mkdir(outputDir)

    """ loop over input files and run SModelS """
    for inputFile in fileList:
        runSingleFile ( inFile, inputFile, outputDir )
        if len(fileList) > 1: inputFile = os.path.join(inFile, inputFile)
        print("Now testing %s" %inputFile)
        currentFile = inputFile
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
        if outputStatus.status < 0: continue

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
            continue

        """ Print Decomposition output.
            If no topologies with sigma > sigmacut are found, update status, write
            output file, stop running """
        if not smstoplist:
            outputStatus.updateStatus(-3)
            continue

        outLevel = 0
        if parser.getboolean("stdout", "printDecomp"):
            outLevel = 1
            outLevel += parser.getboolean("stdout", "addElmentInfo")
        stdoutPrinter.addObj(smstoplist,outLevel)


        """
        Load analysis database
        ======================
        """

        """ In case that a list of analyses or txnames are given, retrieve list """
        analyses = parser.get("database", "analyses").split(",")
        txnames = parser.get("database", "txnames").split(",")
        if parser.get("database", "dataselector") == "efficiencyMap":
            dataTypes = ['efficiencyMap']
            datasetIDs = ['all']
        elif parser.get("database", "dataselector") == "upperLimit":
            dataTypes = ['upperLimit']
            datasetIDs = ['all']
        else:
            dataTypes = ['all']
            datasetIDs = parser.get("database", "dataselector").split(",")
        '''if parser.get("database", "datasets") == "None": datasetIDs = [None]
        else: datasetIDs = parser.get("database", "datasets").split(",")'''

        """ Load analyses """        

        listOfExpRes = database.getExpResults(analysisIDs=analyses, txnames=txnames, 
                            datasetIDs=datasetIDs, dataTypes=dataTypes)

        """ Print list of analyses loaded """
        outLevel = 0
        if parser.getboolean("stdout", "printAnalyses"):
            outLevel = 1
            outLevel += parser.getboolean("stdout", "addAnaInfo")          
        for expResult in listOfExpRes: stdoutPrinter.addObj(expResult,outLevel)
    
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


if __name__ == "__main__":
    import argparse
    """ Set default input and output files """
    parameterFile = "%s/etc/parameters_default.ini" % installDirectory()
    outputDir = "results"

    """ Get the name of input SLHA file and parameter file """
    ap = argparse.ArgumentParser()
    ap.add_argument('-f', '--filename', 
            help='name of SLHA or LHE input file, necessary input, if directory '
                 'is given, loop all files in the directory', required=True)
    ap.add_argument('-p', '--parameterFile', 
            help='name of parameter file, optional argument, if not set, use '
                           'all parameters from etc/parameters_default.ini', 
                           default=parameterFile)
    ap.add_argument('-o', '--outputDir', 
            help='name of output directory, optional argument, default is: ' +
                           outputDir, default=outputDir)
    ap.add_argument('--development', help='enable development output', 
            action='store_true')
    ap.add_argument('-t', '--force_txt', help='force loading the text database',
            action='store_true')
    ap.add_argument('--run-crashreport', 
            help='parse crash report file and use its contents for a SModelS run',
                           action='store_true')
    ap.add_argument('-v','--verbose', help='verbosity level. '
            'accepted values are: debug, info, warning, error.',
                           default = "info", type = str )
    ap.add_argument('--timeout', help='define a limit on the running time (in secs).'
                           ' If not set, run without a time limit', 
                           default = 0, type = int)
    
    
    args = ap.parse_args()

    db=None
    if args.force_txt: db=True
    
    if args.run_crashreport:
        args.filename, args.parameterFile = crashReport.readCrashReportFile(
                args.filename)
        with timeOut.Timeout(args.timeout):
            main(args.filename, args.parameterFile, args.outputDir, args.verbose, db )
        
    else:
        try:
            with timeOut.Timeout(args.timeout):
                main(args.filename, args.parameterFile, args.outputDir, args.verbose, db )
        except Exception:
            crashReportFacility = crashReport.CrashReport()
             
            if args.development:
                print(crashReport.createStackTrace())
            else:
                print(crashReport.createStackTrace())
                crashReportFacility.createCrashReportFile(currentFile, 
                                args.parameterFile)
                print(crashReportFacility.createUnknownErrorMessage())
