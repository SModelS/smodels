"""
.. module:: theory.theoryPrediction
   :synopsis: Provides a class to encapsulate the results of the computation of
              reference cross sections and related functions.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

import copy
from smodels.theory import clusterTools, crossSection, element
from smodels.theory.particleNames import elementsInStr
from smodels.theory.auxiliaryFunctions import cSim, cGtr  # pylint: disable=W0611
from smodels.theory.analysis import SRanalysis
from smodels.theory.analysis import ULanalysis
from smodels.theory.printer import Printer
import logging

logger = logging.getLogger(__name__)


class TheoryPrediction(object):
    """
    An instance of this class represents the results of the theory prediction
    for an analysis.
    
    """
    def __init__(self):
        self.analysis = None
        self.value = None
        self.conditions = None
        self.mass = None
        
    def getmaxCondition(self):      
        maxcond = 0.
        for value in self.conditions.values():
            if value == 'N/A': return value
            if value == None: continue
            maxcond = max(maxcond,value)
        return maxcond        


class TheoryPredictionList(Printer):
    """
    An instance of this class represents the a collection of theory prediction
    objects.
    
    """
    def __init__(self, theoryPredictions):
        self._theoryPredictions = theoryPredictions


    def __iter__(self):
        for theoryPrediction in self._theoryPredictions:
            yield theoryPrediction


    def formatData(self):
        """
        Select data preparation method through dynamic binding.
        
        """
        return Printer.formatTheoryPredictionData(self)


def theoryPredictionFor(analysis, smsTopList, maxMassDist=0.2):
    """
    Compute theory predictions.
    
    Collect elements and efficiencies, combine the masses (if needed) and
    compute the conditions (if existing).
    
    :returns: list of TheoryPrediction objects
    
    """
    # Select elements constrained by analysis and apply efficiencies
    elements = _getElementsFrom(smsTopList, analysis)
    if len(elements) == 0:
        return None
    # Combine masses
    clusters = _combineElements(elements, analysis, maxDist=maxMassDist)
    # Collect results and evaluate conditions
    predictions = []
    for cluster in clusters:
        theoryPrediction = TheoryPrediction()
        theoryPrediction.analysis = analysis
        theoryPrediction.value = _evalConstraint(cluster,analysis)
        theoryPrediction.conditions = _evalConditions(cluster, analysis)
        theoryPrediction.mass = cluster.getAvgMass()
        predictions.append(theoryPrediction)

    if len(predictions) == 0:
        return None
    return TheoryPredictionList(predictions)


def _getElementsFrom(smsTopList, analysis):
    """
    Get elements, that are constrained by the analysis.
    
    Loop over all elements in smsTopList and returns the elements which are
    constrained by the analysis (have efficiency != 0). The elements weights
    are multiplied by their respective efficiency and the cross-sections not
    matching the analysis sqrts are removed.
    
    """
    elements = []
    for el in smsTopList.getElements():
        eff = analysis.getEfficiencyFor(el)
        if eff == 0.:
            continue
        element = el.copy()
        element.weight = crossSection.XSectionList()
        for xsec in el.weight:
            if xsec.info.sqrts == analysis.sqrts:
                element.weight.add(copy.deepcopy(xsec * eff))
        if len(element.weight) > 0:
            elements.append(element)

    return elements


def _combineElements(elements, analysis, maxDist):
    """
    Combine elements according to the analysis type.
    
    If analysis == upper limit type, group elements into mass clusters. If
    analysis == signal region type, group all elements into a single cluster.
    
    """
    if type(analysis) == type(SRanalysis()):
        clusters = [clusterTools.groupAll(elements)]
    elif type(analysis) == type(ULanalysis()):
        clusters = clusterTools.clusterElements(elements, analysis, maxDist)
    return clusters


def _evalConstraint(cluster, analysis):
    """
    Evaluate the analysis constraint inside an element cluster.
    
    If analysis type == upper limit, evaluates the analysis constraint inside
    an element cluster.
    
    :retunrs: total cluster cross-section, if analysis type == signal region
    
    """    
    
    if type(analysis) == type(SRanalysis()):
        return cluster.getTotalXSec()
    elif type(analysis) == type(ULanalysis()):
        if not analysis.constraint:
            return analysis.constraint
        
        exprvalue = _evalExpression(analysis.constraint,cluster,analysis)
        return exprvalue
    

def _evalConditions(cluster, analysis):
    """
    Evaluate the analysis conditions inside an element cluster.
    
    If analysis type == upper limit, evaluates the analysis conditions inside
    an element cluster.
    
    :retunrs: None, if analysis type == signal region
    
    """    
    
    if type(analysis) == type(SRanalysis()):
        return None
    elif type(analysis) == type(ULanalysis()):
        if not analysis.conditions:
            return analysis.conditions
        conditions = {}
        # Loop over conditions
        for cond in analysis.conditions:
            exprvalue = _evalExpression(cond,cluster,analysis)
            if type(exprvalue) == type(crossSection.XSectionList()):
                conditions[cond] = exprvalue[0].value
            else:
                conditions[cond] = exprvalue
                
        return conditions    
        
        
def _evalExpression(stringExpr,cluster,analysis):
    """
    Auxiliary method to evaluate a string expression using the weights of the elements in the cluster.
    
    :return: cross-section value for the expression or expression value if it is not numerical (None,string,...)
    """

#Generate elements appearing in the string expression with zero cross-sections:
    elements = []
    for elStr in elementsInStr(stringExpr):
        el = element.Element(elStr)      
        elements.append(el)

#Replace elements in strings by their weights and add weights from cluster to the elements list:
    expr = stringExpr[:].replace("'","").replace(" ","") 
    for iel, el in enumerate(elements):        
        expr = expr.replace(str(el), "elements["+ str(iel) +"].weight")        
        for el1 in cluster.elements:                    
            if el1.particlesMatch(el):
                el.weight.combineWith(el1.weight)
                el.combineMotherElements(el1) ## keep track of all mothers

    if expr.find("Cgtr") >= 0 or expr.find("Csim") >= 0:
        expr = expr.replace("Cgtr", "cGtr")
        expr = expr.replace("Csim", "cSim")
        logger.warning(analysis.label + " using deprecated functions "
                               "'Cgtr'/'Csim'. Auto-replacing with 'cGtr'/'cSim'.")
    exprvalue = eval(expr)
    if type(exprvalue) == type(crossSection.XSectionList()):
        if len(exprvalue) != 1:
            logger.error("Evaluation of expression "+expr+" returned multiple values.")
        return exprvalue
    else:
        return exprvalue
