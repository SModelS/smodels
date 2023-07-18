#!/usr/bin/env python3

from __future__ import print_function
"""
.. module:: Example
   :synopsis: Basic main file example for using SModelS.
   This file must be run under the installation folder.
"""
""" Import basic functions (this file must be executed in the installation folder) """

from smodels.base import runtime
# Define your model (list of BSM particles)
runtime.modelFile = 'smodels.share.models.mssm'
# runtime.modelFile = 'mssmQNumbers.slha'

from smodels.decomposition import decomposer
from smodels.base.physicsUnits import fb, GeV, TeV
from smodels.matching.theoryPrediction import theoryPredictionsFor,TheoryPredictionsCombiner
from smodels.experiment.databaseObj import Database
from smodels.tools import coverage
from smodels.base.smodelsLogging import setLogLevel
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model
import time
setLogLevel("info")

# Set the path to the database
import os


def main(inputFile='./inputFiles/slha/lightEWinos.slha', sigmacut=0.05*fb,
         database = 'official'):
    """
    Main program. Displays basic use case.
    """

    # Set the path to the database
    database = Database(database)

    model = Model(BSMparticles=BSMList, SMparticles=SMList)
    # Path to input file (either a SLHA or LHE file)
#     lhefile = 'inputFiles/lhe/gluino_squarks.lhe'
    slhafile = os.path.abspath(inputFile)
#     model.updateParticles(inputFile=lhefile)
    model.updateParticles(inputFile=slhafile)

    # Set main options for decomposition
    sigmacut = sigmacut
    mingap = 5.*GeV

    t0 = time.time()
    # Decompose model
    topDict = decomposer.decompose(model, sigmacut,
                                   massCompress=True, invisibleCompress=True,
                                   minmassgap=mingap)

    # Access basic information from decomposition, using the topology list and topology objects:
    print("\n Decomposition done in %1.2fm" %((time.time()-t0)/60.))
    print("\n Decomposition Results: ")
    print("\t  Total number of topologies: %i " % len(topDict))
    nSMS = len(topDict.getSMSList())
    print("\t  Total number of SMS = %i " % nSMS)
    # Print information about the m-th topology:
    m = 2
    if len(topDict) > m:
        cName = sorted(topDict.keys())[m]
        smsList = topDict[cName]
        print("\t\t %i topology  = " % cName)
        # Print information about the n-th element in the m-th topology:
        n = 0
        sms = smsList[n]
        print("\t\t %i-th SMS  = " % (n), sms, end="")
        print("\n\t\t\twith final states =", sms.getFinalStates(), "\n\t\t\twith cross section =", sms.weightList, "\n\t\t\tand masses = ", sms.mass)

    # Load the experimental results to be used.
    # In this case, all results are employed.
    listOfExpRes = database.expResultList

    t0 = time.time()
    # Print basic information about the results loaded.
    # Count the number of loaded UL and EM experimental results:
    nUL, nEM = 0, 0
    for exp in listOfExpRes:
        expType = exp.datasets[0].dataInfo.dataType
        if expType == 'upperLimit':
            nUL += 1
        elif expType == 'efficiencyMap':
            nEM += 1
    print("\n Loaded Database with %i UL results and %i EM results " % (nUL, nEM))

    # Compute the theory predictions for each experimental result and print them:
    print("\n Theory Predictions and Constraints:")
    rmax = 0.
    bestResult = None
    allPredictions = theoryPredictionsFor(database, topDict, combinedResults=False)
    for theoryPrediction in allPredictions:
        print('\n %s ' % theoryPrediction.analysisId())
        dataset = theoryPrediction.dataset
        datasetID = theoryPrediction.dataId()
        txnames = sorted([str(txname) for txname in theoryPrediction.txnames])
        print("------------------------")
        print("Dataset = ", datasetID)  # Analysis name
        print("TxNames = ", txnames)
        print("Theory Prediction = ", theoryPrediction.xsection)  # Signal cross section
        print("Condition Violation = ", theoryPrediction.conditions)  # Condition violation values

        # Get the corresponding upper limit:
        print("UL for theory prediction = ", theoryPrediction.upperLimit)

        # Compute the r-value
        r = theoryPrediction.getRValue()
        print("r = %1.3E" % r)
        # Compute likelihoods for EM-type results:
        if dataset.getType() == 'efficiencyMap':
            theoryPrediction.computeStatistics()
            print('L_BSM, L_SM, L_max = %1.3E, %1.3E, %1.3E' % (theoryPrediction.likelihood(),
                    theoryPrediction.lsm(), theoryPrediction.lmax()))
        if r > rmax:
            rmax = r
            bestResult = theoryPrediction.analysisId()

    # Print the most constraining experimental result
    print("\nThe largest r-value (theory/upper limit ratio) is %1.3E" % rmax)
    if rmax > 1.:
        print("(The input model is likely excluded by %s)" % bestResult)
    else:
        print("(The input model is not excluded by the simplified model results)")

    print("\n Theory Predictions done in %1.2fm" %((time.time()-t0)/60.))
    t0 = time.time()
    # Select a few results results for combination:
    combineAnas = ['ATLAS-SUSY-2013-11', 'CMS-SUS-13-013']
    selectedTheoryPreds = []
    for tp in allPredictions:
        expID = tp.analysisId()
        if expID not in combineAnas:
            continue
        if tp.likelihood() is None:
            continue
        selectedTheoryPreds.append(tp)
    # Make sure each analysis appears only once:
    expIDs = [tp.analysisId() for tp in selectedTheoryPreds]
    if len(expIDs) != len(set(expIDs)):
        print("\nDuplicated results when trying to combine analyses. Combination will be skipped.")
    # Only compute combination if at least two results were selected
    elif len(selectedTheoryPreds) > 1:
        combiner = TheoryPredictionsCombiner(selectedTheoryPreds)
        combiner.computeStatistics()
        llhd = combiner.likelihood()
        lmax = combiner.lmax()
        lsm = combiner.lsm()
        print("\n\nCombined analyses:", combiner.analysisId())
        print("Combined r value: %1.3E" % combiner.getRValue())
        print("Combined r value (expected): %1.3E" % combiner.getRValue(expected=True))
        print("Likelihoods: L, L_max, L_SM = %10.3E, %10.3E, %10.3E\n" % (llhd, lmax, lsm))

    print("\n Combination of analyses done in %1.2fm" %((time.time()-t0)/60.))
    t0 = time.time()
    # Find out missing topologies for sqrts=13*TeV:
    uncovered = coverage.Uncovered(topDict, sqrts=13.*TeV)
    print("\n Coverage done in %1.2fm" %((time.time()-t0)/60.))
    # First sort coverage groups by label
    groups = sorted(uncovered.groups[:], key=lambda g: g.label)
    # Print uncovered cross-sections:
    for group in groups:
        print("\nTotal cross-section for %s (fb): %10.3E\n" % (group.description, group.getTotalXSec()))

    missingTopos = uncovered.getGroup('missing (prompt)')
    # Print some of the missing topologies:
    if missingTopos.finalStateSMS:
        print('Missing topologies (up to 3):')
        for genEl in missingTopos.finalStateSMS[:3]:
            print('Element:', genEl)
            print('\tcross-section (fb):', genEl.missingX)
    else:
        print("No missing topologies found\n")

    missingDisplaced = uncovered.getGroup('missing (displaced)')
    # Print elements with displaced decays:
    if missingDisplaced.finalStateSMS:
        print('\nElements with displaced vertices (up to 2):')
        for genEl in missingDisplaced.finalStateSMS[:2]:
            print('Element:', genEl)
            print('\tcross-section (fb):', genEl.missingX)
    else:
        print("\nNo displaced decays")

    return topDict,allPredictions,uncovered

if __name__ == '__main__':

    topDict,allPredictions,uncovered = main()
