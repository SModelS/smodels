"""
.. module:: theory.theoryPrediction
   :synopsis: Provides a class to encapsulate the results of the computation of
              reference cross sections and related functions.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

from smodels.theory import clusterTools, crossSection, element
from smodels.theory.particleNames import elementsInStr
from smodels.theory.auxiliaryFunctions import cSim, cGtr  #DO NOT REMOVE
from smodels.theory.printer import Printer
import logging

logger = logging.getLogger(__name__)

class TheoryPrediction(object):
    """
    An instance of this class represents the results of the theory prediction
    for an analysis.
    
    :ivar analysis: holds the analysis (ULanalysis or EManalysis object)
                    to which the prediction refers to
    :ivar value: value of the theory prediction 
                (relevant cross-section to be compared with the experimental limits).
                It is a XSection object.
    :ivar conditions: list of values for the analysis conditions
                      (only for upper limit-type analysis, e.g. analysis=ULanalysis)
    :ivar mass: mass of the cluster to which the theory prediction refers to
                (only for upper limit-type analysis, e.g. analysis=ULanalysis)    
    """
    def __init__(self):
        self.analysis = None
        self.value = None
        self.conditions = None
        self.mass = None
        
    def getmaxCondition(self):
        """
        Returns the maximum value from the list conditions
        
        :returns: maximum condition value (float)
        """
            
        maxcond = 0.
        for value in self.conditions.values():
            if value == 'N/A': return value
            if value == None: continue
            maxcond = max(maxcond,value)
        return maxcond        

class TheoryPredictionList(Printer):
    """
    An instance of this class represents a collection of theory prediction
    objects.
    
    :ivar _theoryPredictions: list of TheoryPrediction objects    
    """
    def __init__(self, theoryPredictions):
        """
        Initializes the list.
        
        :parameter theoryPredictions: list of TheoryPrediction objects
        """
        self._theoryPredictions = theoryPredictions


    def __iter__(self):      
        for theoryPrediction in self._theoryPredictions:
            yield theoryPrediction

    def __getitem__(self, index):
        return self._theoryPredictions[index]

    def __len__(self):
        return len(self._theoryPredictions)
    
    def __add__(self,theoPred):
        if isinstance(theoPred,TheoryPrediction):
            self._theoryPredictions.append(theoPred)

    def formatData(self,outputLevel):
        """
        Select data preparation method through dynamic binding.
        :param outputLevel: general control for the output depth to be printed 
                            (0 = no output, 1 = basic output, 2 = detailed output,...
        
        """
        return Printer.formatTheoryPredictionData(self,outputLevel)

def theoryPredictionFor(expResult, smsTopList, maxMassDist=0.2):
    """
    Compute theory predictions for the given experimental result, using the list of elements
    in smsTopList.
    For each Txname appearing in expResult, it collects the elements and efficiencies, 
    combine the masses (if needed) and compute the conditions (if existing).
    
    :parameter expResult: expResult to be considered (ExpResult object)
    :parameter smsTopList: list of topologies containing elements (TopologyList object)
    :parameter maxMassDist: maximum mass distance for clustering elements (float)
    :returns:  a TheoryPredictionList object containing a list of TheoryPrediction objects    
    """
    
    predictionList = TheoryPredictionList()
    #Loop over the txnames in expResult:
    for txname in expResult.txnames:
        # Select elements belonging to TxName and apply efficiencies
        elements = _getElementsFrom(smsTopList, txname)
        if len(elements) == 0: continue        
        #Remove the cross-sections which do not match the experimental analysis:
        for el in elements:
            for xsec in el.weight:
                if xsec.info.sqrts != expResult.info.sqrts:
                    el.weight.delete(xsec)
 
        # Combine masses
        clusters = _combineElements(elements, txname, maxDist=maxMassDist)
        # Collect results and evaluate conditions
        for cluster in clusters:
            theoryPrediction = TheoryPrediction()
            theoryPrediction.txname = txname.txname
            theoryPrediction.value = _evalConstraint(cluster,txname)            
            theoryPrediction.conditions = _evalConditions(cluster, txname)
            theoryPrediction.mass = cluster.getAvgMass()
            predictionList.add(theoryPrediction)

    if len(predictionList) == 0: return None
    else: predictionList


def _getElementsFrom(smsTopList, txname):
    """
    Get elements, that belong to the TxName (appear in constraint).    
    Loop over all elements in smsTopList and returns a copy of the elements belonging
    to TxName (have efficiency != 0). The copied elements
    have their weights multiplied by their respective efficiencies.
    
    :parameter txname: TxName to be considered (TxName object)
    :parameter smsTopList: list of topologies containing elements (TopologyList object)
    :returns: list of elements (Element objects)
    """
    elements = []
    for el in smsTopList.getElements():
        eff = txname.getEfficiencyFor(el)
        if eff == 0.: continue
        element = el.copy()
        element.weight *= eff
    return elements


def _combineElements(elements, txname, maxDist):
    """
    Combine elements according to the txname type.    
    If txname == upper limit type, group elements into mass clusters. If
    txname == efficiency map type, group all elements into a single cluster.
    
    :parameter elements: list of elements (Element objects)
    :parameter analysis: analysis to be considered (ULanalysis or EManalysis object)
    :returns: list of element clusters (ElementCluster objects)
    """
    if txname.txnameData.type == 'efficiencyMap':
        clusters = [clusterTools.groupAll(elements)]
    elif txname.txnameData.type == 'upperLimits':        
        clusters = clusterTools.clusterElements(elements, txname.txnameData, maxDist)
    return clusters


def _evalConstraint(cluster, txname):
    """
    Evaluate the txname constraint inside an element cluster.      
    If txname type == upper limit, sum all the elements' weights
    according to the analysis constraint.
    If txname type == efficiency map, sum all the elements' weights.
    
    :parameter cluster: cluster of elements (ElementCluster object)
    :parameter txname: TxName object holding the constraint
    :returns: cluster cross-section
    """    
    
    if txname.txnameData.type == 'efficiencyMap':
        return cluster.getTotalXSec()
    elif txname.txnameData.type == 'upperLimits':
        if not txname.constraint:
            return txname.constraint
        
        exprvalue = _evalExpression(txname.constraint,cluster)
        return exprvalue
    

def _evalConditions(cluster, txname):
    """        
    If txname type == upper limit (ULanalysis), evaluates the analysis conditions inside
    an element cluster.
    If txname type == efficiency map (EManalysis), returns None
    
    :parameter cluster: cluster of elements (ElementCluster object)
    :parameter txname: TxName object holding the conditions
    :returns: list of condition values (floats) if analysis type == upper limit. None, otherwise.    
    """    
    
    if txname.txnameData.type == 'efficiencyMap':
        return None
    elif txname.txnameData.type == 'upperLimits':        
        if not txname.fuzzycondition:
            return txname.fuzzycondition
        conditions = {}
        # Loop over conditions
        for cond in txname.fuzzycondition:
            exprvalue = _evalExpression(cond,cluster)
            if type(exprvalue) == type(crossSection.XSectionList()):
                conditions[cond] = exprvalue[0].value
            else:
                conditions[cond] = exprvalue
                
        return conditions    
        
        
def _evalExpression(stringExpr,cluster):
    """
    Auxiliary method to evaluate a string expression using the weights of the elements in the cluster.
    Replaces the elements in stringExpr (in bracket notation) by their weights and evaluate the 
    expression.
    e.g. computes the total weight of string expressions such as "[[[e^+]],[[e^-]]]+[[[mu^+]],[[mu^-]]]"
    or ratios of weights of string expressions such as "[[[e^+]],[[e^-]]]/[[[mu^+]],[[mu^-]]]"
    and so on...    
    
    :parameter stringExpr: string containing the expression to be evaluated
    :parameter cluster: cluster of elements (ElementCluster object)
    :returns: value for the expression. Can be a XSectionList object, a float or not numerical (None,string,...)
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
    exprvalue = eval(expr)
    if type(exprvalue) == type(crossSection.XSectionList()):
        if len(exprvalue) != 1:
            logger.error("Evaluation of expression "+expr+" returned multiple values.")
        return exprvalue
    else:
        return exprvalue
