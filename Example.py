#!/usr/bin/env python

from __future__ import print_function

"""
.. module:: Example
   :synopsis: Basic main file example for using SModelS.
   
   This file must be run under the installation folder.

"""

""" Import basic functions (this file must be executed in the installation folder) """

from smodels.theory import slhaDecomposer
from smodels.theory import lheDecomposer
from smodels.tools.physicsUnits import fb, GeV
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.databaseObj import Database

#Set the address of the database folder
database = Database("../smodels-database/")

def main():
    """
    Main program. Displays basic use case.

    """
    
    #Path to input file name (either a SLHA or LHE file)
    slhafile = 'inputFiles/slha/lightEWinos.slha'
#     lhefile = 'inputFiles/lhe/gluino_squarks.lhe'

    #Set main options for decomposition:
    sigmacut = 0.3 * fb
    mingap = 5. * GeV

    """ Decompose model (use slhaDecomposer for SLHA input or lheDecomposer for LHE input) """
    smstoplist = slhaDecomposer.decompose(slhafile, sigmacut, doCompress=True, doInvisible=True, minmassgap=mingap)
    # smstoplist = lheDecomposer.decompose(lhefile, doCompress=True,doInvisible=True, minmassgap=mingap)
    
    # Print basic information from decomposition:
    print("\n\033[40;37m Decomposition Results: \033[0m")
    print("\t \033[31m Total number of topologies: %i \033[0m" %len(smstoplist))
    for top in smstoplist:
        print("\t Topology: ",top)
        print("\t\t Total number of elements =",len(top.elementList))
        print("\t\t Total cross-section =",top.getTotalWeight())
            
    
    
    # Load all analyses from database
    listOfExpRes = database.getExpResults()
    #Count number of loaded UL and EM experimental results:
    nUL, nEM = 0, 0
    for exp in listOfExpRes:
        expType = exp.getValuesFor('dataType')[0]
        if expType == 'upperLimit':
            nUL += 1
        elif  expType == 'efficiencyMap':
            nEM += 1
    #Print basic information about the loaded database:
    print("\nLoaded Database with %i UL results and %i EM results" %(nUL,nEM))

    # Compute the theory predictions for each experimental result and print the results:
    print("\n\033[40;37mTheory Predictions and Constraints:\033[0m")
    rmax = 0.
    for expResult in listOfExpRes:
        predictions = theoryPredictionsFor(expResult, smstoplist)
        if not predictions: continue #Skip if there are no constraints from this result
        print('\n \033[31m %s \033[0m' %expResult.globalInfo.id)
        for theoryPrediction in predictions:
            dataset = theoryPrediction.dataset
            datasetID = dataset.dataInfo.dataId            
            mass = theoryPrediction.mass
            txnames = [str(txname) for txname in theoryPrediction.txnames]
            PIDs =  theoryPrediction.PIDs         
            print("------------------------")
            print("Dataset = ",datasetID)   #Analysis name
            print("TxNames = ",txnames)   
            print("Prediction Mass = ",mass)    #Value for average cluster mass (average mass of the elements in cluster)
            print("Prediction PIDs = ",PIDs)    #Value for average cluster mass (average mass of the elements in cluster)
            print("Theory Prediction = ",theoryPrediction.xsection)   #Signal cross-section
            print("Condition Violation = ",theoryPrediction.conditions)  #Condition violation values
              
            #Get upper limit for the respective prediction:
            if expResult.getValuesFor('dataType')[0] == 'upperLimit':
                ul = expResult.getUpperLimitFor(txname=theoryPrediction.txnames[0],mass=mass)                     
            elif expResult.getValuesFor('dataType')[0] == 'efficiencyMap':
                ul = expResult.getUpperLimitFor(dataID=datasetID)
            else: print('weird:',expResult.getValuesFor('dataType'))
            print("Theory Prediction UL = ",ul)
            r = theoryPrediction.xsection.value/ul
            print("R = ",r)
            rmax = max(r,rmax)
            
    print("\nThe largest r-value (theory/upper limit ratio) is ",rmax)
    if rmax > 1.:
        print("\033[31m (The input model is likely excluded) \033[0m")
    else:
        print("\033[32m (The input model is not excluded by the simplified model results) \033[0m")
      
    


if __name__ == '__main__':
    main()
