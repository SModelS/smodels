"""
.. module:: theory.theoryPrediction
   :synopsis: Provides a class to encapsulate the results of the computation of
              reference cross sections and related functions.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

from smodels.theory import clusterTools, crossSection, element
from smodels.theory.particleNames import elementsInStr
from smodels.theory.auxiliaryFunctions import cSim, cGtr  #DO NOT REMOVE
import logging,sys
from smodels.tools.physicsUnits import TeV

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

        if not self.conditions: return 0.        
        maxcond = 0.
        for value in self.conditions.values():
            if value == 'N/A': return value
            if value == None: continue
            maxcond = max(maxcond,value)
        return maxcond        

class TheoryPredictionList(object):
    """
    An instance of this class represents a collection of theory prediction
    objects.
    
    :ivar _theoryPredictions: list of TheoryPrediction objects   
     
    """
    def __init__(self, theoryPredictions=None):
        """
        Initializes the list.
        
        :parameter theoryPredictions: list of TheoryPrediction objects
        """        
        self._theoryPredictions = []
        if theoryPredictions and isinstance(theoryPredictions,list):
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
    
    dataSetResults = []
    #Compute predictions for each data set (for UL analyses there is one single set)
    for dataset in expResult.datasets:
        predList = _getDataSetPredictions(dataset,smsTopList,maxMassDist)
        if predList: dataSetResults.append(predList)
    if not dataSetResults: return None
    
    #For results with more than one dataset, select the best data set 
    #according to the expect upper limit
    bestResults = _getBestResults(dataSetResults)
    bestResults.expResult = expResult
    for op in bestResults:
        op.expResult = expResult
    return bestResults


def _getBestResults(dataSetResults):
    """
    Returns the best result according to the expected upper limit
    :param dataSetResults: list of TheoryPredictionList objects
    :return: best result (TheoryPredictionList object)
    """
        
    #In the case of UL analyses or efficiency-maps with a single signal region
    #return the single result:        
    if len(dataSetResults) == 1:
        return dataSetResults[0]
    
    #For efficiency-map analyses with multipler signal regions,
    #select the best one according to the expected upper limit:
    bestExpectedR = 0.
    for predList in dataSetResults:
        dataset = predList.dataset
        if dataset.dataInfo.datatype != 'efficiency-map':
            logger.error("Multiple data sets should only exist for efficiency map results!")
            sys.exit()
        if len(predList) != 1:
            logger.error("Multiple clusters should only exist for upper limit results!")
            sys.exit()
        pred = predList[0]
        if len(pred.value) != 1:
            logger.error("Signal region prediction should correspond to a single cross-section!")
            sys.exit()
        xsec = pred.value[0]        
        expectedR = xsec.value/dataset.getUpperLimit(0.05,True)
        if expectedR > bestExpectedR:
            bestExpectedR = expectedR
            bestPredList = predList
    
    return bestPredList
    

def _getDataSetPredictions(dataset,smsTopList,maxMassDist):   
    """
    Compute theory predictions for a given data set.
    For upper-limit results returns the list of theory predictions for the experimental result.
    For efficiency-map results returns the list of theory predictions for the signal region.
    Uses the list of elements in smsTopList.
    For each Txname appearing in dataset, it collects the elements and efficiencies, 
    combine the masses (if needed) and compute the conditions (if existing).
    
    :parameter dataset: Data Set to be considered (DataSet object)
    :parameter smsTopList: list of topologies containing elements (TopologyList object)
    :parameter maxMassDist: maximum mass distance for clustering elements (float)
    :returns:  a TheoryPredictionList object containing a list of TheoryPrediction objects    
    """     
    
    predictionList = TheoryPredictionList()
    predictionList.dataset = dataset
    # Select elements belonging to expResult and apply efficiencies
    elements = _getElementsFrom(smsTopList, dataset)
    if len(elements) == 0: return None
    if (dataset.globalInfo.sqrts/TeV).normalize()._unit:
            ID = dataset.globalInfo.id
            logger.error("Sqrts defined with wrong units for %s" %(ID) )
            return False

    for el in elements:
        for xsec in el.weight:          
            if xsec.info.sqrts != dataset.globalInfo.sqrts:
                el.weight.delete(xsec)    

    # Combine elements according to their respective constraints and masses
    # (For efficiencyMap analysis group all elements)
    clusters = _combineElements(elements, dataset, maxDist=maxMassDist)
    
    # Collect results and evaluate conditions    
    for cluster in clusters:
        theoryPrediction = TheoryPrediction()
        theoryPrediction.txname = cluster.txname  # = None for efficiency map results
        theoryPrediction.value = _evalConstraint(cluster)            
        theoryPrediction.conditions = _evalConditions(cluster)
        theoryPrediction.mass = cluster.getAvgMass()
        theoryPrediction.PIDs = cluster.getPIDs()
        predictionList + theoryPrediction
        

    

    if len(predictionList) == 0: return None
    else: return predictionList


def _getElementsFrom(smsTopList, dataset):
    """
    Get elements that belong to any of the TxNames in dataset 
    (appear in any of constraints in the result).    
    Loop over all elements in smsTopList and returns a copy of the elements belonging
    to any of the constraints (TxNames, i.e. have efficiency != 0). The copied elements
    have their weights multiplied by their respective efficiencies.
    
    :parameter dataset:  Data Set to be considered (DataSet object)
    :parameter smsTopList: list of topologies containing elements (TopologyList object)
    :returns: list of elements (Element objects)
    
    """
    
    elements = []
    for el in smsTopList.getElements():             
        for txname in dataset.txnameList:
            hasEl = txname.hasElementAs(el)
            if not hasEl: continue   
            eff = txname.getEfficiencyFor(hasEl)            
            if eff == 0.: continue
            hasEl.weight *= eff
            elements.append(hasEl)
            break  
    return elements


def _combineElements(elements, dataset, maxDist):
    """
    Combine elements according to the data set type.    
    If expResult == upper limit type, first group elements with different TxNames
    and then into mass clusters.
    If expResult == efficiency map type, group all elements into a single cluster.
    
    :parameter elements: list of elements (Element objects)
    :parameter expResult: Data Set to be considered (DataSet object)
    :returns: list of element clusters (ElementCluster objects)
    
    """
    
    clusters = []
    if dataset.dataInfo.datatype == 'efficiency-map':
        cluster = clusterTools.groupAll(elements)        
        clusters.append(cluster)
    elif dataset.dataInfo.datatype == 'upper-limit':        
        for txname in dataset.txnameList:
            txnameEls = []
            for element in elements:
                #Check if element really belongs to txname:
                if not txname.getEfficiencyFor(element): continue
                else: txnameEls.append(element)
            txnameClusters = clusterTools.clusterElements(txnameEls, txname, maxDist)            
            clusters += txnameClusters
    else:
        logger.warning("Unkown data type: %s. Data will be ignored." 
                       % dataset.dataInfo.datatype)
                
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
        if not cluster.txname.fuzzycondition \
         or cluster.txname.fuzzycondition == "not yet assigned":
            return cluster.txname.fuzzycondition
        conditions = {}
        #Make sure conditions is always a list
        if isinstance(cluster.txname.fuzzycondition,str): 
            cluster.txname.fuzzycondition = [cluster.txname.fuzzycondition]
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

#Get cross-section info from cluster (to generate zero cross-section values):
    infoList = cluster.elements[0].weight.getInfo()    
#Generate elements appearing in the string expression with zero cross-sections:
    elements = []
    for elStr in elementsInStr(stringExpr):
        el = element.Element(elStr)
        el.weight = crossSection.XSectionList(infoList)
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
