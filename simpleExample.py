#!/usr/bin/env python

"""
.. module:: simpleExample
   :synopsis: Basic use case for the SModelS framework.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from __future__ import print_function
from theory import slhaDecomposer
from theory import lheDecomposer
from tools.physicsUnits import addunit
from experiment import smsAnalysisFactory
from theory.theoryPrediction import theoryPredictionFor
import logging

logger = logging.getLogger(__name__) # pylint: disable-msg=C0103

slhafile = 'inputFiles/slha/andrePT4.slha'
lhefile = 'inputFiles/lhe/ued_1.lhe'


def main():
    """
    Main program. Displays basic use case.
    
    """
    # Decompose model (SLHA or LHE input):
    sigmacut = addunit(0.1,'fb')
    mingap = addunit(5.,'GeV')
    
    smsTopList = slhaDecomposer.decompose(slhafile, sigmacut, doCompress=True,
                                          doInvisible=True, minmassgap=mingap)
    # smsTopList = lheDecomposer.decompose(lhefile, doCompress=True,
    #                                      doInvisible=True, minmassgap=mingap)
    
    # Print decomposition summary
    smsTopList.printout()
    
    # Load analyses
    listOfAnalyses = smsAnalysisFactory.load()
    
    # Get theory prediction for each analysis and print basic output       
    for analysis in listOfAnalyses:
        theoryPredictions = theoryPredictionFor(analysis, smsTopList)
        if not theoryPredictions:
            continue
        print(analysis.label)
        theoryPredictions.printout()
    

if __name__ == '__main__':
    main()
