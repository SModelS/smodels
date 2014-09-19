#!/usr/bin/env python

"""
.. module:: simpleExample
   :synopsis: Basic use case for the SModelS framework.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from __future__ import print_function
import setPath  # # set to python path for smodels
# from smodels import SModelS
from smodels.theory import slhaDecomposer
from smodels.theory import lheDecomposer
from smodels.tools.physicsUnits import fb, GeV
from smodels.experiment import smsAnalysisFactory
from smodels.theory.theoryPrediction import theoryPredictionFor
# import logging


def main():
    """
    Main program. Displays basic use case.

    """

#     slhafile = 'doublet_123f981ec.slha'
    slhafile = 'inputFiles/slha/andrePT4.slha'
#     lhefile = 'inputFiles/lhe/jory_dqdqbar_mgo_inf_msq300_mlsp100.lhe'

    # Decompose model (SLHA or LHE input):
    sigmacut = 0.001 * fb
    mingap = 5. * GeV

    smstoplist = slhaDecomposer.decompose(slhafile, sigmacut, doCompress=True,
                                        doInvisible=True, minmassgap=mingap)
#     smstoplist = lheDecomposer.decompose(lhefile, doCompress=True,
#                                         doInvisible=True, minmassgap=mingap)

    # Print decomposition summary
    smstoplist.printout()
    for top in smstoplist:
      for el in top.elementList:
        print(el,el.weight)

    # Load analyses
    listofanalyses = smsAnalysisFactory.load()

    # Get theory prediction for each analysis and print basic output
    for analysis in listofanalyses:
        theorypredictions = theoryPredictionFor(analysis, smstoplist)
        if not theorypredictions:
            continue
        for cluster in theorypredictions:
            resdic = {}
            resdic['mass'] = cluster.mass
            resdic['conditions'] = cluster.conditions
            resdic['theoryPrediction'] = str(cluster.value)
            resdic['analysisName'] = cluster.analysis.label
            resdic['UpperLimit'] = cluster.analysis.getUpperLimitFor(cluster.mass)
            print(resdic)
        theorypredictions.printout()


if __name__ == '__main__':
    main()
