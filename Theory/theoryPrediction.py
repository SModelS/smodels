#!/usr/bin/env python

"""
.. module:: TheoryPrediction
        :synopsis: Classes encapsulating the results of the computation of reference cross sections and related methods
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

import copy, time
import clusterTools
import crossSection, element, analysis
import logging
from auxiliaryFunctions import elementsInStr, Csim, Cgtr
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class TheoryPrediction(object):
    """Main class to store the results of the theory prediction for a given analysis."""
    
    def __init__(self):
        self.analysis = None
        self.value = None
        self.conditions = None
        self.mass = None


def theoryPredictionFor(Analysis,SMSTopList,maxMassDist=0.2):
    """Main method to compute theory predictions. Collects the elements and efficiencies, combine the masses
    (if needed) and compute the conditions (if any). Returns a list of TheoryPrediction objects"""

#Select elements constrained by analysis and apply efficiencies
    t1 = time.time()                        
    elements = getElementsFrom(SMSTopList,Analysis)
    print 'getElementsFrom done in',time.time()-t1,'s'   
    if len(elements) == 0: return None      
#Combine masses
    t1 = time.time()
    print 'nels=',len(elements)    
    clusters = combineElements(elements,Analysis,maxDist=maxMassDist)
    print 'combineElements done in',time.time()-t1,'s'
#Collect results and evaluate conditions:
    predictions = []
    for cluster in clusters:
        theoPrediction = TheoryPrediction()
        theoPrediction.value = cluster.getTotalXSec()
        theoPrediction.conditions = evalConditions(cluster,Analysis)
        theoPrediction.mass = cluster.getAvgMass()
        predictions.append(theoPrediction)

    if len(predictions) == 0: return None
    return predictions


def getElementsFrom(SMSTopList,Analysis):
    """ Loops over all elements in SMSTopList and returns the elements which are constrained by the analysis (have
    efficiency != 0). The elements weights are multiplied by their respective efficiency and the cross-sections not
    matching the analysis sqrts are removed."""
    
    elements = []
    for el in SMSTopList.getElements():
        eff = Analysis.getEfficiencyFor(el)
        if eff == 0.: continue
        element = el.copy()
        element.weight = crossSection.XSectionList()
        for xsec in el.weight:
            if xsec.info.sqrts == Analysis.sqrts: element.weight.add(copy.deepcopy(xsec*eff))
        if len(element.weight) > 0: elements.append(element)
            
    return elements


def combineElements(elements,Analysis,maxDist):
    """ Combines elements according to the Analysis type. If Analysis = upper limit type, group elements
    into mass clusters. If Analysis = signal region type, group all elements into a single cluster"""
    
    if type(Analysis) == type(analysis.SRanalysis()):
        clusters = [clusterTools.groupAll(elements)]
    elif type(Analysis) == type(analysis.ULanalysis()):
        clusters = clusterTools.clusterElements(elements,Analysis,maxDist)
    return clusters


def evalConditions(cluster,Analysis):
    """If analysis type = upper limit, evaluates the analysis conditions inside an element cluster.
    If  analysis type = signal region, returns None"""
        
    if type(Analysis) == type(analysis.SRanalysis()): return None
    elif type(Analysis) == type(analysis.ULanalysis()):
        if not Analysis.conditions: return Analysis.conditions
        conditions = {}
#Loop over conditions        
        for cond in Analysis.conditions:
            condElements = [element.Element(el_str) for el_str in elementsInStr(cond)]  #Get elements appearing in conditions            
            newcond = cond
            for iel,el in enumerate(condElements):
                newcond = newcond.replace(str(el),"condElements["+str(iel)+"].weight")
                for el1 in cluster.elements:   
                    if el1.particlesMatch(el):                        
                        el.weight.combineWith(el1.weight)
                   
            conditions[cond] = eval(newcond)
            
        return conditions
    

