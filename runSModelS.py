#!/usr/bin/env python

"""
.. module:: Main code for running smodels
"""

import sys, os, logging
import argparse
from ConfigParser import SafeConfigParser
from smodels.tools.physicsUnits import GeV, fb
from smodels.tools import ioObjects, missingTopologies
from smodels.experiment import smsHelpers, smsAnalysisFactory
from smodels.theory import slhaDecomposer, lheDecomposer
from smodels.theory.theoryPrediction import theoryPredictionFor
from smodels.installation import installDirectory
log = logging.getLogger(__name__)


def main(inputFile, parameterFile, outputFile):
    """
    Provides a command line interface to basic SModelS functionalities.
    
    :param inputFile: input file name (either a SLHA or LHE file)
    :param parameterFile: File containing the input parameters (default = /etc/parameters_default.ini)
    :param outputFile: Output file to write a summary of results
    """

#-------------------------
#Read and check input file
#-------------------------
    #use ConfigParser to read input parameters
    parser = SafeConfigParser()
    parser.read(parameterFile)
    
    #sigmacut = minimum value of cross-section for an element to be considered eligible for decomposition.
    #Too small sigmacut leads to too large decomposition time.
    sigmacut = parser.getfloat("parameters","sigmacut")*fb
    #minmassgap = minimum value for considering two states non-degenerate (only used for mass compression)
    minmassgap = parser.getfloat("parameters","minmassgap")*GeV

    if os.path.exists(outputFile): #remove old output file
        log.warning("Removing old output file in "+outputFile)
    outfile = open(outputFile,'w')
    outfile.close()

    inputType = parser.get("options","inputType").lower()
    if inputType != 'slha' and inputType != 'lhe':
        log.error("Unknown input type (must be SLHA or LHE): %s" % inputType)
        sys.exit()
    #Check input file for errors
    inputStatus = ioObjects.FileStatus()
    if parser.getboolean("options","checkInput"): inputStatus.checkFile(inputType,inputFile,sigmacut)        

    #check database address    
    try:
        smsHelpers.base = parser.get("path","databasePath")
        databaseVersion = smsHelpers.databaseVersion()
    except:
        log.error("Database not found in %s" % os.path.realpath(smsHelpers.base))
        databaseVersion  = None
        sys.exit()

    #initialize output status and exit if there were errors in the input
    outputStatus = ioObjects.OutputStatus(inputStatus.status,inputFile,dict(parser.items("parameters")), databaseVersion, outputFile)
#---------------------------------------------------------------------------------------------------------


#-------------------------
#Decompose input file
#-------------------------
    try:
        # Decompose input SLHA file, store the output elements in smstoplist
        if inputType == 'slha':
            smstoplist = slhaDecomposer.decompose(inputFile, sigmacut, doCompress=parser.getboolean("options","doCompress"),
                         doInvisible=parser.getboolean("options","doInvisible"), minmassgap=minmassgap)
        else:
            smstoplist = lheDecomposer.decompose(inputFile, doCompress=parser.getboolean("options","doCompress"),
                         doInvisible=parser.getboolean("options","doInvisible"), minmassgap=minmassgap)   
    except:
        #Update status to fail, print error message and exit:
        outputStatus.updateStatus(-1)
        
#Print Decomposition output:
    # If no topologies with sigma > sigmacut are found, update status, write output file, stop running
    if not smstoplist:
        outputStatus.updateStatus(-3)
            
    outLevel= 0
    if parser.getboolean("stdout","printDecomp"):
        outLevel = 1
        outLevel += parser.getboolean("stdout","addElmentInfo")
    smstoplist.printout(outputLevel=outLevel)
#------------------------

#-------------------------
#Load analysis database
#-------------------------    
    #in case that a list of analyses or txnames are given, retrieve list
    analyses = parser.get("database","analyses")
    if "," in analyses: analyses = analyses.split(",")    
    txnames = parser.get("database","txnames")
    if "," in txnames: txnames = txnames.split(",")
    # Load analyses
    listofanalyses = smsAnalysisFactory.load(analyses, txnames)

#Print list of analyses loaded:
    if parser.getboolean("stdout","printAnalyses"):
        outLevel = 1
        outLevel += parser.getboolean("stdout","addAnaInfo")
        print("=======================\n == List of Analyses   ====\n ================")
        for analysis in listofanalyses:
            analysis.printout(outputLevel=outLevel)    
#------------------------------


#----------------------------------------------------
#Compute theory predictions and anlalyses constraints
#----------------------------------------------------
    
    #define result list that collects all theoryPrediction objects
    #variables set to define printing options
    results = ioObjects.ResultList(bestresultonly = not parser.getboolean("file","expandedSummary"), 
                                   describeTopo = parser.getboolean("file","addConstraintInfo"))
    
    #Get theory prediction for each analysis and print basic output
    for analysis in listofanalyses:        
        theorypredictions = theoryPredictionFor(analysis, smstoplist)
        if not theorypredictions:  continue
        if parser.getboolean("stdout","printResults"):
            print "================================================================================"
            theorypredictions.printout()
        print "................................................................................"

        # Create a list of results, to determine the best result
        for theoryprediction in theorypredictions:
            results.addResult(theoryprediction, maxcond = parser.getfloat("parameters","maxcond"))

    # If there is no best result, this means that there are no matching experimental results for the point.
    if results.isEmpty():
        # no experimental constraints found
        outputStatus.updateStatus(0)
    else:
        outputStatus.updateStatus(1)

    # write output file
    outputStatus.printout("file", outputFile)
    # add experimental constraints if found
    if outputStatus.status == 1: results.printout("file", outputFile)

    sqrts = max([xsec.info.sqrts for xsec in smstoplist.getTotalWeight()])
    if parser.getboolean("options","findMissingTopos"):
        #look for missing topologies, add them to the output file
        missingtopos = missingTopologies.MissingTopoList(sqrts)
        missingtopos.findMissingTopos(smstoplist, listofanalyses, minmassgap, parser.getboolean("options","doCompress"),
                         doInvisible=parser.getboolean("options","doInvisible"))
        missingtopos.printout("file", outputFile)
# ---------------------------------------------------------


if __name__ == "__main__":
    
    #set default input and output files
    parameterFile = "%s/etc/parameters_default.ini" % installDirectory()
    outputFile = "summary.txt"
    
    #get the name of input slha file (and parameter file)
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-f', '--filename', help = 'name of SLHA or LHE input file, necessary input', required = True)
    argparser.add_argument('-p', '--parameterFile', 
                            help = 'name of parameter file, optional argument, if not set, use all parameters from etc/parameters_default.ini',
                            default = parameterFile)
    argparser.add_argument('-o', '--outputFile', help = 'name of output file, optional argument, default is: '+outputFile, 
                           default = outputFile)
    args = argparser.parse_args() # pylint: disable-msg=C0103

    # execute run function
    main(args.filename, args.parameterFile, args.outputFile)
