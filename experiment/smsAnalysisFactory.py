#!/usr/bin/env python

"""
.. module:: smsAnalysisFactory
   :synopsis: Creates a list of Analysis objects from a results database.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from __future__ import print_function
from . import smsResults
from tools.physicsUnits import rmvunit
from theory import analysis
from theory import element
import logging

logger = logging.getLogger(__name__)

def getRealTopo(tx):
    """
    T3w025 -> T3w, etc.
    
    """
    ret = tx
    ret.replace("050","").replace("x1C180","").replace("025","")
    if ret.find("x") > -1:
        ret = ret[:ret.find("x")]
    return ret


def getArray(constraint):
    """
    Extracts the number of vertices, branches and insertions from the
    constraint string. It maps, e.g.,
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


def getElementsEffs(constraint):
    """
    Generates a dictionary of elements with their simple efficiencies as values
    for an upper limit-type of analysis. Equals the element multiplicative
    factor appearing in the constraint.
    
    """      
    # Get element strings appearing in constraint
    elStrings = getArray(constraint)
    cons = constraint.replace(" ", "")
    cons = cons.replace("'", "")
    elementsEff = {}
    for el in elStrings:
        newcons = cons.replace(el, "1.", 1)
        for el2 in elStrings:
            if el2 != el:
                newcons = newcons.replace(el2, "0.", 1)
        elEff = eval(newcons)
        elementsEff[element.Element(el)] = elEff        
    return elementsEff        
          

def load(anas=None, topos=None, sqrts=[7, 8]):
    """
    Creates the Analysis objects from the info given in the SMS results
    database.

    :param anas: if given as a list, then we create only objects for these
    analyses (the database naming convention is used).
    :param topos: if given as a list, then only these topos are considered.
    :param sqrts: array of center-of-mass energies of the analyses that are to
    be considered.
    :returns: list of Analyses
    
    """   
    # lets make sure the user can also supply a single topo/ana without having
    # to code an array
    if isinstance(topos, str):
        topos = [topos]
    if isinstance(anas, str):
        anas = [anas]
    if isinstance(sqrts, int):
        sqrts = [sqrts]
    if isinstance(sqrts, float):
        sqrts = [int(sqrts)]

    listOfAnalyses = []

        
    if anas == None:
        anas = smsResults.getAllResults().keys()
    for ana in anas:
        logger.debug("Building ana " + str(ana))
        ss = rmvunit(smsResults.getSqrts(ana), "TeV")
        if ss == None:
            logger.debug("SS = " + str(ss) + str(ana))
            continue
        ss = int(ss)
        if not ss in sqrts:
            continue
        for tx in smsResults.getTopologies(ana):
            if topos != None and tx not in topos:
                continue
            logger.debug(str(tx))                        
            newAnalysis = analysis.ULanalysis()
            newAnalysis.sqrts = smsResults.getSqrts(ana)
            stopo = getRealTopo (tx)
            newAnalysis.label = ana + ":" + tx
            # "2012"
            newAnalysis.run = smsResults.getRun(ana)
            constraint = smsResults.getConstraints(ana, topo=stopo)
            cond = smsResults.getConditions(ana, topo=stopo)
            if not constraint or constraint == "Not yet assigned":
                logger.debug("dont have a constraint for " + str(ana) + \
                            str(tx) + "(" + str(stopo) + ")")
                continue
            
            if cond:
                cond = cond.replace("'", "").replace(" ", "").split(';')
            newAnalysis.constraint = constraint
            newAnalysis.conditions = cond
            newAnalysis.elementsEff = getElementsEffs(constraint)
            # Add ana to list of analyses:
            listOfAnalyses.append(newAnalysis)

    return listOfAnalyses


if __name__ == "__main__":
    load()
    print("List of analyses/results: ")
    listOfAnalyses = load()
    for (ct, ana) in enumerate(listOfAnalyses):
        # .label, ana.sqrts
        print(ct, ana.label, ana.Top.vertnumb, ana.Top.vertparts)
