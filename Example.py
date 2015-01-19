#!/usr/bin/env python

"""
.. module:: Example
   :synopsis: Basic main file example for using SModelS.

"""

#Import basic functions (this file must be run under the installation folder)
from __future__ import print_function
import sys
from smodels.theory import slhaDecomposer
from smodels.theory import lheDecomposer
from smodels.tools.physicsUnits import fb, GeV
from smodels.experiment import smsAnalysisFactory
from smodels.theory.theoryPrediction import theoryPredictionFor
from smodels.experiment.databaseObjects import DataBase

#Set the address of the database folder
database = DataBase("/home/lessa/smodels-database/")


def main():
    """
    Main program. Displays basic use case.

    """

    #Path to input file name (either a SLHA or LHE file)
    slhafile = 'inputFiles/slha/gluino_squarks.slha'
#     lhefile = 'inputFiles/lhe/gluino_squarks.lhe'

    #Set main options for decomposition:
    sigmacut = 0.3 * fb
    mingap = 5. * GeV

    #Decompose model (use slhaDecomposer for SLHA input or lheDecomposer for LHE input):
    smstoplist = slhaDecomposer.decompose(slhafile, sigmacut, doCompress=True,doInvisible=True, minmassgap=mingap)
#     smstoplist = lheDecomposer.decompose(lhefile, doCompress=True,doInvisible=True, minmassgap=mingap)

    # Print decomposition summary. Set outputLevel=0 (no output), 1 (simple output), 2 (extended output)
    smstoplist.printout(outputLevel=1)

    # Load all analyses from database
    listofanalyses = database.getAnalyses()

    # Compute the theory predictions for each analysis
    analysesPredictions = [theoryPredictionFor(analysis, smstoplist) for analysis in listofanalyses]

    #Access information for each theory prediction/analysis
    for analysisPred in analysesPredictions:
        if not analysisPred: continue  #Skip non-applicable analyses
        #If the analysis prediction contains more than one theory prediction (cluster), loop over predictions:        
        for theoryPrediction in analysisPred:
            print("------------------------")
            print("Analysis name = ",theoryPrediction.analysis.label)   #Analysis name
            print("Prediction Mass = ",theoryPrediction.mass)    #Value for average cluster mass (average mass of the elements in cluster)
            print("Signal Cross-Section = ",theoryPrediction.value)   #Value for the cluster signal cross-section
            print("Condition Violation = ",theoryPrediction.conditions)  #Condition violation values
            
            #Get upper limit for the respective prediction:
            print("Analysis UL = ",theoryPrediction.analysis.getUpperLimitFor(theoryPrediction.mass)) 
    
    


if __name__ == '__main__':
    main()
    sys.exit()
