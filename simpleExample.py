#!/usr/bin/env python

"""SModelS basic use case.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from tools import set_path
import tools.loggingConfiguration
tools.loggingConfiguration.configure()
from theory import slhaDecomposer, lheDecomposer
from tools.physicsUnits import addunit
from experiment import smsAnalysisFactory
from theory.theoryPrediction import theoryPredictionFor
import logging
#logger = logging.getLogger(__name__)

def main():
    """
    Main program. Displays basic use case.
    
    """

#Decompose model (SLHA or LHE input):    
    slhafile = "inputFiles/slha/andrePT4.slha"
    lhefile = "inputFiles/lhe/ued_1.lhe"
    mingap = addunit(5.,'GeV')
    sigmacut = addunit(0.1,'fb')
    smsTopList = slhaDecomposer.decompose(slhafile, sigmacut, doCompress=True
                                          , doInvisible=True, minmassgap=mingap)
#     smsTopList = lheDecomposer.decompose(lhefile, doCompress=True, doInvisible=True, minmassgap=mingap)
#Print decomposition summary
    smsTopList.printout()
    
#Load analyses
    listOfAnalyses = smsAnalysisFactory.load()
#Get theory prediction for each analysis and print basic output       
    for ana in listOfAnalyses:
        preds = theoryPredictionFor(ana, smsTopList)
        if not preds: continue
        print ana.label
        for pred in preds:
            print 'mass=', pred.mass
            print 'theory prediction=', pred.value
            print 'theory conditions:'
            if not pred.conditions:
                print pred.conditions
            else:
                for cond in pred.conditions:
                    print pred.conditions[cond]
            print 'experimental limit=',ana.getUpperLimitFor(pred.mass)
            print '\n'
    

if __name__ == '__main__':
    main()
