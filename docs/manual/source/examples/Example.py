#!/usr/bin/env python3

from __future__ import print_function
"""
.. module:: Example
   :synopsis: Basic main file example for using SModelS.   
   This file must be run under the installation folder.
"""
""" Import basic functions (this file must be executed in the installation folder) """


from smodels.tools import runtime
from smodels.theory import decomposer
from smodels.tools.physicsUnits import fb, GeV, TeV
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.databaseObj import Database
from smodels.tools import coverage
from smodels.tools.smodelsLogging import setLogLevel
from smodels.particlesLoader import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model
setLogLevel("info")

# Set the path to the database
database = Database("official")

def main():
    """
    Main program. Displays basic use case.
    """
    #Define your model (list of BSM particles)
    runtime.modelFile = 'smodels.share.models.mssm'
#     runtime.modelFile = 'mssmQNumbers.slha'
    model = Model(BSMparticles=BSMList, SMparticles=SMList)

    # Path to input file (either a SLHA or LHE file)
#     lhefile = 'inputFiles/lhe/gluino_squarks.lhe'
    slhafile = 'inputFiles/slha/lightEWinos.slha'
#     model.updateParticles(inputFile=lhefile)
    model.updateParticles(inputFile=slhafile)


    # Set main options for decomposition
    sigmacut = 0.01*fb
    mingap = 5.*GeV

    # Decompose model
    toplist = decomposer.decompose(model, sigmacut, doCompress=True, doInvisible=True, minmassgap=mingap)
    
    # Access basic information from decomposition, using the topology list and topology objects:
    print("\n Decomposition Results: " )
    print("\t  Total number of topologies: %i " %len(toplist) )
    nel = sum([len(top.elementList) for top in toplist])
    print("\t  Total number of elements = %i " %nel )
    #Print information about the m-th topology:
    m = 2
    if len(toplist) > m:
        top = toplist[m]
        print( "\t\t %i-th topology  = " %m,top,"with total cross section =",top.getTotalWeight() )
        #Print information about the n-th element in the m-th topology:
        n = 0
        el = top.elementList[n]
        print( "\t\t %i-th element from %i-th topology  = " %(n,m),el, end="" )
        print( "\n\t\t\twith final states =",el.getFinalStates(),"\n\t\t\twith cross section =",el.weight,"\n\t\t\tand masses = ",el.mass )

    # Load the experimental results to be used.
    # In this case, all results are employed.
    listOfExpRes = database.getExpResults()

    # Print basic information about the results loaded.
    # Count the number of loaded UL and EM experimental results:
    nUL, nEM = 0, 0
    for exp in listOfExpRes:
        expType = exp.datasets[0].dataInfo.dataType
        if expType == 'upperLimit':
            nUL += 1
        elif  expType == 'efficiencyMap':
            nEM += 1
    print("\n Loaded Database with %i UL results and %i EM results " %(nUL,nEM) )


    # Compute the theory predictions for each experimental result and print them:
    print("\n Theory Predictions and Constraints:")
    rmax = 0.
    bestResult = None
    for expResult in listOfExpRes:
        predictions = theoryPredictionsFor(expResult, toplist, combinedResults=False, marginalize=False)
        if not predictions: continue # Skip if there are no constraints from this result
        print('\n %s ' %expResult.globalInfo.id)
        for theoryPrediction in predictions:
            dataset = theoryPrediction.dataset
            datasetID = theoryPrediction.dataId()
            mass = theoryPrediction.mass
            txnames = [str(txname) for txname in theoryPrediction.txnames]
            PIDs =  theoryPrediction.PIDs         
            print("------------------------" )
            print("Dataset = ",datasetID )   #Analysis name
            print("TxNames = ",txnames )  
            print("Prediction Mass = ",mass )   #Value for average cluster mass (average mass of the elements in cluster)
            print("Prediction PIDs = ",PIDs )   #Value for average cluster mass (average mass of the elements in cluster)
            print("Theory Prediction = ",theoryPrediction.xsection )  #Signal cross section
            print("Condition Violation = ",theoryPrediction.conditions ) #Condition violation values
              
            # Get the corresponding upper limit:
            print("UL for theory prediction = ",theoryPrediction.upperLimit )

            # Compute the r-value
            r = theoryPrediction.getRValue()
            print("r = ",r )
            #Compute likelihhod and chi^2 for EM-type results:
            if dataset.getType() == 'efficiencyMap':
                theoryPrediction.computeStatistics()
                print('Chi2, likelihood=', theoryPrediction.chi2, theoryPrediction.likelihood )
            if r > rmax:
                rmax = r
                bestResult = expResult.globalInfo.id
            
    # Print the most constraining experimental result
    print("\nThe largest r-value (theory/upper limit ratio) is ",rmax )
    if rmax > 1.:
        print("(The input model is likely excluded by %s)" %bestResult )
    else:
        print("(The input model is not excluded by the simplified model results)" )
      
    #Find out missing topologies for sqrts=8*TeV:
    uncovered = coverage.Uncovered(toplist,sqrts=8.*TeV)
    #First sort coverage groups by label
    groups = sorted(uncovered.groups[:], key = lambda g: g.label)
    #Print uncovered cross-sections:
    for group in groups:
        print("\nTotal cross-section for %s (fb): %10.3E\n" %(group.description,group.getTotalXSec()))

    missingTopos = uncovered.getGroup('missing (prompt)')
    #Print some of the missing topologies:
    if missingTopos.generalElements:
        print('Missing topologies (up to 3):' )
        for genEl in missingTopos.generalElements[:3]:
            print('Element:', genEl)
            print('\tcross-section (fb):', genEl.missingX)
    else:
        print("No missing topologies found\n")
    
    missingDisplaced = uncovered.getGroup('missing (displaced)')
    #Print elements with displaced decays:
    if missingDisplaced.generalElements:
        print('\nElements with displaced vertices (up to 2):' )    
        for genEl in missingDisplaced.generalElements[:2]:
            print('Element:', genEl)
            print('\tcross-section (fb):', genEl.missingX)
    else:
        print("\nNo displaced decays")
      
        
if __name__ == '__main__':
    main()
