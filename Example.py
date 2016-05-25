#!/usr/bin/env python

from __future__ import print_function

"""
.. module:: Example
   :synopsis: Basic main file example for using SModelS.
   
   This file must be run under the installation folder.

"""

""" Import basic functions (this file must be executed in the installation folder) """

import sys
from smodels.theory import slhaDecomposer
from smodels.theory import lheDecomposer
from smodels.tools.physicsUnits import fb, GeV
from smodels.tools.printer import printout
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.databaseObj import Database

#Set the address of the database folder
database = Database("./smodels-database/")

def main():
    """
    Main program. Displays basic use case.

    """
    
    #Path to input file name (either a SLHA or LHE file)
#     slhafile = 'inputFiles/slha/gluino_squarks.slha'
    slhafile = 'inputFiles/slha/lightEWinos.slha'
#     lhefile = 'inputFiles/lhe/gluino_squarks.lhe'

    #Set main options for decomposition:
    sigmacut = 0.3 * fb
    mingap = 5. * GeV

    """ Decompose model (use slhaDecomposer for SLHA input or lheDecomposer for LHE input) """
    smstoplist = slhaDecomposer.decompose(slhafile, sigmacut, doCompress=True, doInvisible=True, minmassgap=mingap)
    # smstoplist = lheDecomposer.decompose(lhefile, doCompress=True,doInvisible=True, minmassgap=mingap)
    
    # Print decomposition summary. Set outputLevel=0 (no output), 1 (simple output), 2 (extended output)    
    printout(smstoplist,outputLevel=2)
    
    
    # Load all analyses from database
    listOfExpRes = database.getExpResults(dataTypes=['efficiencyMap'])

    # Compute the theory predictions for each analysis
    for expResult in listOfExpRes:
        predictions = theoryPredictionsFor(expResult, smstoplist)
        if not predictions: continue
        print('\n',expResult.getValuesFor('id')[0])
        for theoryPrediction in predictions:
            dataset = theoryPrediction.dataset
            datasetID = dataset.getValuesFor('dataId')[0]            
            mass = theoryPrediction.mass
            txnames = [str(txname) for txname in theoryPrediction.txnames]
            PIDs =  theoryPrediction.PIDs         
            print("------------------------")
            print("Dataset = ",datasetID)   #Analysis name
            print("TxNames = ",txnames)   
            print("Prediction Mass = ",mass)    #Value for average cluster mass (average mass of the elements in cluster)
            print("Prediction PIDs = ",PIDs)    #Value for average cluster mass (average mass of the elements in cluster)
            print("Theory Prediction = ",theoryPrediction.value)   #Value for the cluster signal cross-section
            print("Condition Violation = ",theoryPrediction.conditions)  #Condition violation values
              
            #Get upper limit for the respective prediction:
            if expResult.getValuesFor('dataType')[0] == 'upperLimit':
                ul = expResult.getUpperLimitFor(txname=theoryPrediction.txnames[0],mass=mass)                     
            elif expResult.getValuesFor('dataType')[0] == 'efficiencyMap':
                ul = expResult.getUpperLimitFor(dataID=datasetID)
            else: print('weird:',expResult.getValuesFor('dataType'))
            print("Theory Prediction UL = ",ul)
            print("R = ",theoryPrediction.value[0].value/ul)
      
    


if __name__ == '__main__':
    main()
