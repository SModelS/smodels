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
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.databaseObjects import DataBase

#Set the address of the database folder
database = DataBase("../smodels-database/")


def main():
    """
    Main program. Displays basic use case.

    """

    #Path to input file name (either a SLHA or LHE file)
#     slhafile = 'inputFiles/slha/gluino_squarks.slha'
    slhafile = 'inputFiles/slha/lightEWinos.slha'
#     slhafile = 'inputFiles/slha/compression.slha'
#     lhefile = 'inputFiles/lhe/gluino_squarks.lhe'

    #Set main options for decomposition:
    sigmacut = 0.3 * fb
    mingap = 5. * GeV

    """ Decompose model (use slhaDecomposer for SLHA input or lheDecomposer for LHE input) """
    smstoplist = slhaDecomposer.decompose(slhafile, sigmacut, doCompress=True, doInvisible=True, minmassgap=mingap)
    # smstoplist = lheDecomposer.decompose(lhefile, doCompress=True,doInvisible=True, minmassgap=mingap)

    # Print decomposition summary. Set outputLevel=0 (no output), 1 (simple output), 2 (extended output)
    smstoplist.printout(outputLevel=0)

    # Load all analyses from database
    listOfExpRes = [database.getExpResults(datasetIDs=[None],analysisIDs=['ATLAS-SUSY-2013-11'])]
#     listOfExpRes = database.getExpResults(datasetIDs=[None])

    # Compute the theory predictions for each analysis
    for expResult in listOfExpRes:
        print('\n',expResult)
        predictions = theoryPredictionsFor(expResult, smstoplist)
        if not predictions: continue
        dataset = predictions.dataset
        datasetID = dataset.getValuesFor('dataid')
        for theoryPrediction in predictions:
            mass = theoryPrediction.mass
            txname = theoryPrediction.txname            
            print("------------------------")
            print("TxName = ",txname)   #Analysis name
            print("Prediction Mass = ",mass)    #Value for average cluster mass (average mass of the elements in cluster)
            print("Theory Prediction = ",theoryPrediction.value)   #Value for the cluster signal cross-section
            print("Condition Violation = ",theoryPrediction.conditions)  #Condition violation values
              
            #Get upper limit for the respective prediction:
            if expResult.getValuesFor('datatype') == 'upper-limit':
                print("Theory Prediction UL = ",expResult.getUpperLimitFor(txname=txname,mass=mass))
            elif expResult.getValuesFor('datatype') == 'efficiency-map':
                print("Theory Prediction UL = ",expResult.getUpperLimitFor(dataID=datasetID))
            else: print('weird:',expResult.getValuesFor('type'))
      
    


if __name__ == '__main__':
    main()
