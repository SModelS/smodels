#!/usr/bin/env python

"""
.. module:: runSModelS
   :synopsis: Main code for running SModelS.
   
"""

from __future__ import print_function
import os, sys
import logging
import argparse
from ConfigParser import SafeConfigParser
from smodels.experiment.databaseObjects import Database
from smodels.installation import installDirectory
from smodels.theory import slhaDecomposer
from smodels.theory import lheDecomposer
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.tools.physicsUnits import GeV
from smodels.tools.physicsUnits import fb
from smodels.tools import ioObjects
from smodels.tools import missingTopologies
from smodels.tools import crashReport
import smodels.tools.printer as prt
from smodels.experiment.exceptions import DatabaseNotFoundException

log = logging.getLogger(__name__)

def main(inputFile, parameterFile, outputFile):
    """
    Provides a command line interface to basic SModelS functionalities.
    
    :param inputFile: input file name (either a SLHA or LHE file)
    :param parameterFile: File containing the input parameters (default = /etc/parameters_default.ini)
    :param outputFile: Output file to write a summary of results
    
    """

    """
    Read and check input file
    =========================
    """
    parser = SafeConfigParser()
    parser.read(parameterFile)

    """ Minimum value of cross-section for an element to be considered eligible for decomposition.
        Too small sigmacut leads to too large decomposition time. """
    sigmacut = parser.getfloat("parameters", "sigmacut") * fb

    """ Minimum value for considering two states non-degenerate (only used for mass compression) """
    minmassgap = parser.getfloat("parameters", "minmassgap") * GeV

    if os.path.exists(outputFile):
        log.warning("Removing old output file " + outputFile)
    outfile = open(outputFile, 'w')
    outfile.close()


    inputType = parser.get("options", "inputType").lower()
    if inputType != 'slha' and inputType != 'lhe':
        log.error("Unknown input type (must be SLHA or LHE): %s" % inputType)
        return

    """ Check input file for errors """
    inputStatus = ioObjects.FileStatus()
    if parser.getboolean("options", "checkInput"):
        inputStatus.checkFile(inputType, inputFile, sigmacut)

    """ Check database location """
    try:
        databasePath = parser.get("path", "databasePath")
        database = Database(databasePath)
        databaseVersion = database.databaseVersion
    except DatabaseNotFoundException:
        log.error("Database not found in %s" % os.path.realpath(databasePath))
        databaseVersion = None
        return

    """ Initialize output status and exit if there were errors in the input """
    outputStatus = ioObjects.OutputStatus(inputStatus.status, inputFile, dict(parser.items("parameters")), databaseVersion, outputFile)
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
            smstoplist = slhaDecomposer.decompose(inputFile, sigmacut, doCompress=parser.getboolean("options", "doCompress"),
                         doInvisible=parser.getboolean("options", "doInvisible"), minmassgap=minmassgap)
        else:
            smstoplist = lheDecomposer.decompose(inputFile, doCompress=parser.getboolean("options", "doCompress"),
                         doInvisible=parser.getboolean("options", "doInvisible"), minmassgap=minmassgap)
    except:
        """ Update status to fail, print error message and exit """
        outputStatus.updateStatus(-1)
        return

    """ Print Decomposition output.
        If no topologies with sigma > sigmacut are found, update status, write output file, stop running """
    if not smstoplist:
        outputStatus.updateStatus(-3)
        return

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

    listOfExpRes = database.getExpResults(analysisIDs=analyses, txnames=txnames, datasetIDs=datasetIDs, dataTypes=dataTypes)

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

    """ If there is no best result, this means that there are no matching experimental results for the point """
    if not allPredictions:
        """ no experimental constraints found """
        outputStatus.updateStatus(0)
    else:
        outputStatus.updateStatus(1)

    stdoutPrinter.addObj(outputStatus)
    summaryPrinter.addObj(outputStatus)
    
    
    """ Define result list that collects all theoryPrediction objects."""
    maxcond = parser.getfloat("parameters", "maxcond")
    results = ioObjects.ResultList(allPredictions,maxcond)
    if not parser.getboolean("file", "expandedSummary"):
        results.useBestResult()

    outLevel = 0
    if not results.isEmpty():
        outLevel = 1
        outLevel += parser.getboolean("file", "addConstraintInfo")
    summaryPrinter.addObj(results,outLevel)
    if parser.getboolean("stdout", "printResults"):
        stdoutPrinter.addObj(results,outLevel)
    
    sqrts = max([xsec.info.sqrts for xsec in smstoplist.getTotalWeight()])
    if parser.getboolean("options", "findMissingTopos"):
        """ Look for missing topologies, add them to the output file """
        missingtopos = missingTopologies.MissingTopoList(sqrts)
        missingtopos.findMissingTopos(smstoplist, listOfExpRes, minmassgap, parser.getboolean("options", "doCompress"),
                      doInvisible=parser.getboolean("options", "doInvisible"))        
        summaryPrinter.addObj(missingtopos)
        stdoutPrinter.addObj(missingtopos,2)
         
        
    stdoutPrinter.close()
    summaryPrinter.close()


if __name__ == "__main__":
    """ Set default input and output files """
    parameterFile = "%s/etc/parameters_default.ini" % installDirectory()
    outputFile = "summary.txt"

    """ Get the name of input SLHA file and parameter file """
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-f', '--filename', help='name of SLHA or LHE input file, necessary input', required=True)
    argparser.add_argument('-p', '--parameterFile', help='name of parameter file, optional argument, if not set, use '
                           'all parameters from etc/parameters_default.ini', default=parameterFile)
    argparser.add_argument('-o', '--outputFile', help='name of output file, optional argument, default is: ' +
                           outputFile, default=outputFile)
    argparser.add_argument('--development', help='enable development output', action='store_true')
    argparser.add_argument('--run-crashreport', help='parse crash report file and use its contents for a SModelS run',
                           action='store_true')
    args = argparser.parse_args()
    
    if args.run_crashreport:
        args.filename, args.parameterFile = crashReport.readCrashReportFile(args.filename)
        main(args.filename, args.parameterFile, args.outputFile)
        
    else:
        try:
            main(args.filename, args.parameterFile, args.outputFile)
        except Exception:
            crashReportFacility = crashReport.CrashReport()
             
            if args.development:
                print(crashReport.createStackTrace())
            else:
                print(crashReport.createStackTrace())
                crashReportFacility.createCrashReportFile(args.filename, args.parameterFile)
                print(crashReportFacility.createUnknownErrorMessage())
