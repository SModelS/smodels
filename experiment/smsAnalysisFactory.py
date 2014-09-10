#!/usr/bin/env python

"""
.. module:: experiment.smsAnalysisFactory
   :synopsis: Create a list of analysis objects from a results database.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from __future__ import print_function
from smodels.experiment import smsResults, smsHelpers
from smodels.tools.physicsUnits import rmvunit
from smodels.theory.particleNames import elementsInStr
from smodels.theory import analysis
from smodels.theory import element
import logging

logger = logging.getLogger(__name__)


def load(analyses=None, topologies=None, sqrts=[7, 8], usePrivate=None):
    """
    Create an analysis objects from the info given in the SMS results database.
)
    :param analyses: If a list of analyses is passed, create only objects for
                     these analyses (the database naming convention is used).
    :param topologies: If a list of topologies is passed, only these topologies
                       are considered.
    :param sqrts: Array of center-of-mass energies of the analyses that are to
                  be considered.
    :param usePrivate: Include analyses flagged as private.
    :returns: list of analyses
    
    """
    # Enable supplying a single analysis/topology
    if isinstance(topologies, str):
        topologies = [topologies]
    if isinstance(analyses, str):
        analyses = [analyses]
    if isinstance(sqrts, int):
        sqrts = [sqrts]
    if isinstance(sqrts, float):
        sqrts = [int(sqrts)]

    listOfAnalyses = []

    if analyses == None:
        analyses = smsResults.getAllResults().keys()
    for ana in analyses:
        if smsResults.isPrivate(ana) and not usePrivate:
            logger.warning("Skipping private analysis %s.",str(ana))
            continue
        logger.debug("Building analysis %s.", str(ana))
        ss = rmvunit(smsResults.getSqrts(ana), "TeV")
        if ss == None:
            logger.debug("SS: " + str(ss) + str(ana))
            continue
        ss = int(ss)
        if not ss in sqrts:
            continue
        for tx in smsResults.getTopologies(ana):
            if topologies != None and tx not in topologies:
                continue
            logger.debug(str(tx))
            newAnalysis = analysis.ULanalysis()
            newAnalysis.sqrts = smsResults.getSqrts(ana)
            stopo = _getRealTopo(tx)
            newAnalysis.label = ana + ":" + tx
            # "2012"
            newAnalysis.run = smsHelpers.getRun(ana)
            constraint = smsResults.getConstraints(ana, topology=stopo)
            cond = smsResults.getConditions(ana, topology=stopo)
            if not constraint or constraint == "Not yet assigned":
                logger.debug("Constraint for %s %s (%s) does not exist.",
                             str(ana), str(tx), str(stopo))
                continue

            if cond:
                cond = cond.replace("'", "").replace(" ", "").split(';')
            newAnalysis.constraint = constraint
            newAnalysis.conditions = cond
            newAnalysis.elementsEff = _getElementsEffs(constraint,cond)
            # Add analysis to list of analyses:
            listOfAnalyses.append(newAnalysis)

    return listOfAnalyses


def _getRealTopo(tx):
    """
    Get real topology, e.g., T3w025 -> T3w, etc.
    
    """
    ret = tx
    ret.replace("050", "").replace("x1C180", "").replace("025", "")
    if ret.find("x") > -1:
        ret = ret[:ret.find("x")]
    return ret


def _getElementsEffs(constraint,conditions):
    """
    Generate a dictionary of elements with their simple efficiencies as values.    
    Efficiencies are = 1. if the element appears in the constraint or conditions.
    
    """
    # Get element strings appearing in constraint
    elStrings = elementsInStr(constraint)
    if conditions:
        for cond in conditions: elStrings += elementsInStr(cond)
    elementsEff = {}
    elStrings = set(elStrings)
    for elstr in elStrings:
        el = element.Element(elstr)
        elementsEff[el] = 1.
    return elementsEff



def _getArray(constraint):
    """
    Get number of vertices, branches and insertions from a constraint string.
    
    This maps, e.g.,
    2*([[['L'],['L']],[['L'],['nu']]] + [[['L'],['L']],[['nu'],['L']]])
    to
    [[['L'],['L']],[['L'],['nu']]]
    
    """
    c = constraint.replace(" ", "")
    c = c.replace("'", "")
    elStrings = []
    while "[[[" in c:
        el = c[c.find("[[["):c.find("]]]") + 3]
        c = c.replace(el, "", 1)
        elStrings.append(el)
    return elStrings


if __name__ == "__main__":
    load()
    print("List of analyses/results: ")
    _listOfAnalyses = load()
    for (ct, _ana) in enumerate(_listOfAnalyses):
        print(ct, _ana.label, _ana.Top.vertnumb, _ana.Top.vertparts)
