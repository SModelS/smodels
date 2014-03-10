#!/usr/bin/env python

"""
SModelS basic use case.

"""

from theory import slhaDecomposer
from tools.PhysicsUnits import addunit
from experiment import smsanalysisFactory
from theory.theoryPrediction import theoryPredictionFor
import logging

logger = logging.getLogger(__name__)

def main():
    """
    Main program. Displays basic use case.
    
    """
    # useXsec = CrossSection.XSectionInfo()
    # useXsec.sqrts = addunit(8,'TeV')
    # useXsec.order = 2
    # useXsec.label = 'tev8'
    # UseXSecs = [useXsec]
    
    listOfAnalyses = smsanalysisFactory.load()
    slhafile = "slha/andrePT4.slha"
    # lhefile = "lhe/ued_1.lhe"
    # lhefile = "lhe/TChiChipmSlepL_1.lhe"
    # nevts = 10000
    doCompress = True
    doInvisible = True
    minmassgap = addunit(5.,'GeV')
    sigmacut = addunit(0.1,'fb')
    smsTopList = slhaDecomposer.decompose(slhafile, sigmacut, doCompress,
                                          doInvisible, minmassgap)
    # smsTopList = lheDecomposer.decompose(lhefile, None, None, doCompress,
    #                                      doInvisible, minmassgap)

    smsTopList.printout()
    
    for ana in listOfAnalyses:    
        preds = theoryPredictionFor(ana, smsTopList)
        if not preds:
            continue
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
            print '\n'
    

if __name__ == '__main__':
    main()
