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
from smodels.tools.physicsUnits import GeV, fb
from smodels.tools.modelTester import runAllFiles
from smodels.tools import crashReport, timeOut
import smodels.tools.printer as prt
from smodels.experiment.exceptions import DatabaseNotFoundException

log = logging.getLogger(__name__)

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

    """ Setup output printers """
    stdoutPrinter = prt.TxTPrinter(output = 'stdout')

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

    runAllFiles ( fileList, inFile, outputDir, parser, databaseVersion, listOfExpRes )


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
    
    if args.run_crashreport: #FIXME overall crash report for the loop (before reading any input file)
        args.filename, args.parameterFile = crashReport.readCrashReportFile(
                args.filename)
        with timeOut.Timeout(args.timeout):
            main(args.filename, args.parameterFile, args.outputDir, args.verbose, db )
        
    else:
        try:
            with timeOut.Timeout(args.timeout):
                main( args.filename, args.parameterFile, args.outputDir, 
                      args.verbose, db )
        except Exception,e:
            crashReportFacility = crashReport.CrashReport()
             
            if args.development:
                print(crashReport.createStackTrace())
            else:
                print(crashReport.createStackTrace())
                crashReportFacility.createCrashReportFile( e.inputFile, 
                                args.parameterFile )
                print(crashReportFacility.createUnknownErrorMessage())
