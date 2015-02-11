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
import logging,sys

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

def theoryPredictionsFor(expResult, smsTopList, maxMassDist=0.2):
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
    
    predictionList = []
    # Select elements belonging to expResult and apply efficiencies
    elements = _getElementsFrom(smsTopList, expResult)
    if len(elements) == 0: return None
    for el in elements:
        for xsec in el.weight:
            if xsec.info.sqrts != expResult.info.sqrts:
                el.weight.delete(xsec)    

    # Combine elements according to their respective constraints and masses
    # (For efficiencyMap analysis group all elements)
    clusters = _combineElements(elements, expResult, maxDist=maxMassDist)
    # Collect results and evaluate conditions
    for cluster in clusters:
        theoryPrediction = TheoryPrediction()
        theoryPrediction.expResult = expResult
        theoryPrediction.txname = cluster.txname  # = None for efficiency map results
        theoryPrediction.value = _evalConstraint(cluster)            
        theoryPrediction.conditions = _evalConditions(cluster)
        theoryPrediction.mass = cluster.getAvgMass()
        predictionList.append(theoryPrediction)

    if len(predictionList) == 0: return None
    else: return TheoryPredictionList(predictionList)


def _getElementsFrom(smsTopList, expResult):
    """
    Get elements that belong to any of the TxNames in expResult 
    (appear in any of constraints in the result).    
    Loop over all elements in smsTopList and returns a copy of the elements belonging
    to any of the constraints (TxNames, i.e. have efficiency != 0). The copied elements
    have their weights multiplied by their respective efficiencies.
    
    :parameter expResult: Experimental result to be considered (ExpResult object)
    :parameter smsTopList: list of topologies containing elements (TopologyList object)
    :returns: list of elements (Element objects)
    """
    
    elements = []
    for el in smsTopList.getElements():
        for txname in expResult.txnames:   
            eff = txname.getEfficiencyFor(el)       
            if eff == 0.: continue
            element = el.copy()
            element.weight *= eff
            elements.append(element)  
    return elements


def _combineElements(elements, expResult, maxDist):
    """
    Combine elements according to the experimental result type.    
    If expResult == upper limit type, first group elements with different TxNames
    and then into mass clusters.
    If expResult == efficiency map type, group all elements into a single cluster.
    
    :parameter elements: list of elements (Element objects)
    :parameter expResult: experimental result to be considered (ExpResult object)
    :returns: list of element clusters (ElementCluster objects)
    """
    
    clusters = []
    if expResult.txnames[0].txnameData.type == 'efficiencyMap':
        cluster = clusterTools.groupAll(elements)
        cluster.txname = None
        clusters.append(cluster)
    elif expResult.txnames[0].txnameData.type == 'upperLimits':        
        for txname in expResult.txnames:
            txnameEls = []
            for element in elements:
                for el in txname._elements:                
                    if element.particlesMatch(el):
                        txnameEls.append(element)
                        break
            txnameClusters = clusterTools.clusterElements(txnameEls, txname.txnameData, maxDist)
            for cluster in txnameClusters: cluster.txname = txname
            clusters += txnameClusters
    return clusters


def _evalConstraint(cluster):
    """
    Evaluate the constraint inside an element cluster.      
    If the cluster refers to a specific TxName, sum all the elements' weights
    according to the analysis constraint.
    If cluster.txname = None (efficiency map), sum all the elements' weights.
    
    :parameter cluster: cluster of elements (ElementCluster object)
    :returns: cluster cross-section
    """    
    
    if cluster.txname is None:
        return cluster.getTotalXSec()
    elif cluster.txname.txnameData.type == 'upperLimits':
        if not cluster.txname.constraint or cluster.txname.constraint == "not yet assigned":
            return cluster.txname.constraint
        exprvalue = _evalExpression(cluster.txname.constraint,cluster)
        return exprvalue
    else:
        logger.error("Unknown cluster type")
        sys.exit()
    

def _evalConditions(cluster):
    """
    Evaluate the conditions (if any) inside an element cluster.   
    If the cluster refers to a specific TxName, evaluates the analysis conditions inside
    an element cluster.
    If cluster.txname = None (efficiency map), returns None
    
    :parameter cluster: cluster of elements (ElementCluster object)    
    :returns: list of condition values (floats) if analysis type == upper limit. None, otherwise.    
    """    
    
    if cluster.txname is None:
        return None
    elif cluster.txname.txnameData.type == 'upperLimits':  
        if not cluster.txname.fuzzycondition or cluster.txname.fuzzycondition == 'None'\
         or cluster.txname.fuzzycondition == "not yet assigned":
            return cluster.txname.fuzzycondition
        conditions = {}
        # Loop over conditions
        for cond in cluster.txname.fuzzycondition:
            exprvalue = _evalExpression(cond,cluster)
            if type(exprvalue) == type(crossSection.XSectionList()):
                conditions[cond] = exprvalue[0].value
            else:
                conditions[cond] = exprvalue                
        return conditions    
    else:
        logger.error("Unknown cluster type")
        sys.exit()
        
        
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
