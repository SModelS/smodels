"""
.. module:: limitGetter
   :synopsis: Access the proper experimental limits to given analysis objects.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""

from . import smsResults
from tools.physicsUnits import rmvunit
import copy
import sys
import logging

logger = logging.getLogger(__name__)

def limit(analysis, addTheoryPredictions=[]):
    """
    Get limit from an analysis object.

    :param addTheoryPredictions: list of theory predictions to add, e.g.,
    [ '7 TeV (NLL)', '7 TeV (LO)' ]
    :type addTheoryPredictions: [String]
    """
    sqrts = rmvunit(analysis.sqrts,"TeV")
    # lead=analysis.Top.leadingElement()
    ret = []
    for (constraint, condition) in analysis.results.items():
        if len(addTheoryPredictions) > 0:
            if not analysis.computeTheoryPredictions() or \
                    len(analysis.ResultList) == 0:
                continue
            theoRes = analysis.ResultList[0]
        tx = analysis.plots[constraint][0]
        for ana in analysis.plots[constraint][1]:
            for (index, element) in enumerate(analysis.Top.elements()):
                for (mi, masses1) in enumerate(element.B[0].masses):
                    # masses1=element.B[0].masses[0] # ,lead.B[1].masses[0]
                    masses2 = element.B[1].masses[mi] # ,lead.B[1].masses[0]
                    ul = smsResults.getSmartUpperLimit(ana, tx, masses1,
                                                       masses2)
                    tmp = {"ul": ul, "analysis": ana, "tx": tx, "m1": masses1,
                           "m2": masses2, "sqrts": sqrts}
                    if len(addTheoryPredictions) > 0:
                        theory = theoRes.prediction()
                        tmp["theory"] = theory
                        allexcl = False
                        for t in addTheoryPredictions:
                            excl = rmvunit(theory[t], "fb") > rmvunit(ul, "fb")
                            tmp["excl_%s" % t] = excl
                            allexcl = allexcl or excl
                        tmp["excluded"] = allexcl
                    ret.append(tmp)
    return ret

        
def getPlotLimit(inmass, analysis, complain=False):
    """
    Get upper limit on sigma*BR for a specific array of masses from plot.
    
    :param inmass: Array of masses in SModelS graph.
    :param analysis: experiment.analysis.ULanalysis.
    
    """
    massArray = copy.deepcopy(inmass)

    # Skip empty mass arrays
    if len(massArray) < 2: 
        logger.info('M = ' + str(massArray))
        sys.exit()
        if complain:
            logger.error("Length of mass-array < 2.")
        return False

    # Make sure the two branches have equal masses
    if massArray[0] != massArray[1]:
        if complain:
            logger.error("Masses differ between branches.")
        return False

    masslist = [rmvunit(mass,'GeV') for mass in massArray[0]]

    # Run label
    run = analysis.run
    # If run has not been defined, use latest run
    if run == "":
        run = None
     
    analysis, cmsLabel = analysis.label.split(':')    
    limit = smsResults.getSmartUpperLimit(analysis, cmsLabel, masslist, run,
                                          debug=complain) 
    return limit
