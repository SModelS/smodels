#!/usr/bin/env python

from __future__ import print_function

"""
.. module:: printExample
   :synopsis: Basic file example for exempligyin how to use the printers SModelS.
   
   This file must be run under the installation folder.

"""

""" Import basic functions (this file must be executed in the installation folder) """

import sys
from smodels.theory import slhaDecomposer
from smodels.tools.physicsUnits import fb, GeV
from smodels.tools.printer import printout
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.databaseObjects import Database
from smodels.tools import printer, ioObjects
from smodels.tools import modpyslha as pyslha

#Set the address of the database folder
database = Database("../smodels-database/")


stdoutPrinter = printer.TxTPrinter(output = 'stdout')
summaryPrinter = printer.SummaryPrinter(output = 'file', filename = 'summary_print.txt')
pythonPrinter = printer.PyPrinter(output = 'file', filename = 'sms_output.py')
printerList = printer.MPrinter(stdoutPrinter,summaryPrinter,pythonPrinter)

def main():
    """
    Main program. Displays basic use case.

    """
    
    #Path to input file name (either a SLHA or LHE file)
#     slhafile = 'inputFiles/slha/gluino_squarks.slha'
    slhafile = 'inputFiles/slha/lightEWinos.slha'
#     slhafile = 'inputFiles/slha/compression.slha'
#     lhefile = 'inputFiles/lhe/gluino_squarks.lhe'

    #Add some file information to the printer
    pyslhaData = pyslha.readSLHAFile(slhafile)
    printerList.addObj(pyslhaData)
    
    
    #Set main options for decomposition:
    sigmacut = 0.3 * fb
    mingap = 5. * GeV

    """ Decompose model (use slhaDecomposer for SLHA input or lheDecomposer for LHE input) """
    smstoplist = slhaDecomposer.decompose(slhafile, sigmacut, doCompress=True, doInvisible=True, minmassgap=mingap)
    # smstoplist = lheDecomposer.decompose(lhefile, doCompress=True,doInvisible=True, minmassgap=mingap)

    #Add the decomposition result to the printers
    printerList.addObj(smstoplist)

    listOfExpRes = database.getExpResults()

    # Compute the theory predictions for each analysis
    for expResult in listOfExpRes:        
        predictions = theoryPredictionsFor(expResult, smstoplist)
        printerList.addObj(predictions)

    
    #Add additional information:
    outputStatus = ioObjects.OutputStatus([1,None], slhafile,
                                           {'sigmacut' : sigmacut, 'mingap' : mingap},
                                           'databaseVersion','outputfile')
    printerList.addObj(outputStatus)
    printerList.close()
    


if __name__ == '__main__':
    main()
