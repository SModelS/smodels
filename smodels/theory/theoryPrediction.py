"""
.. module:: theoryPrediction
   :synopsis: Provides a class to encapsulate the results of the computation of
              reference cross sections and related functions.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. autofunction:: _getElementsFrom    
"""

from smodels.theory import clusterTools, crossSection, element
from smodels.theory.particleNames import elementsInStr
from smodels.theory.auxiliaryFunctions import cSim, cGtr  #DO NOT REMOVE
import sys
from smodels.tools.physicsUnits import TeV,fb
from smodels.theory.exceptions import SModelSTheoryError as SModelSError

from smodels.tools.smodelsLogging import logger


class TheoryPrediction(object):
    """
    An instance of this class represents the results of the theory prediction
    for an analysis.
    
    :ivar analysis: holds the analysis (ULanalysis or EManalysis object)
                    to which the prediction refers
    :ivar xsection: xsection of the theory prediction 
                (relevant cross section to be compared with the experimental limits).
                It is a XSection object.
    :ivar conditions: list of values for the analysis conditions
                      (only for upper limit-type analysis, e.g. analysis=ULanalysis)
    :ivar mass: mass of the cluster to which the theory prediction refers
                (only for upper limit-type analysis, e.g. analysis=ULanalysis)   
                 
    """
    def __init__(self):
        self.analysis = None
        self.xsection = None
        self.conditions = None
        self.mass = None        

    def computeStatistics(self):
        """
        Compute the likelihood, chi-square and expected upper limit for this theory prediction.
        The resulting values are stored as the likelihood, chi2 and expectedUL attributes.        
        """
        if not hasattr(self, "dataset") or self.dataset.dataInfo.dataType == 'upperLimit':
            self.likelihood = None
            self.chi2 = None
            self.expectedUL = None
            return
        
        lumi = self.dataset.globalInfo.lumi  
        nsig = (self.xsection.value*lumi).asNumber()
        llhd = self.dataset.likelihood(nsig)
        chi2 = self.dataset.chi2(nsig)
        expectedUL = self.dataset.getSRUpperLimit(alpha = 0.05, expected = True)
        
        self.likelihood =  llhd
        self.chi2 =  chi2
        self.expectedUL = expectedUL
        
    def getmaxCondition(self):
        """
        Returns the maximum xsection from the list conditions
        
        :returns: maximum condition xsection (float)        
        """

        if not self.conditions: return 0.        
        values = [ 0. ]
        for value in self.conditions.values():
            if value == 'N/A': return value
            if value == None: continue
            values.append ( float(value) )
        return max(values)
        
    
    def __str__(self):
        return "%s:%s" % ( self.analysis, self.xsection )



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
    
    def __add__(self,theoPredList):
        if isinstance(theoPredList,TheoryPredictionList):
            res = TheoryPredictionList()
            res._theoryPredictions = self._theoryPredictions + theoPredList._theoryPredictions
            return res
        else:
            return None
        
    def __radd__(self, theoPredList):
        if theoPredList == 0:
            return self
        else:
            return self.__add__(theoPredList)        

def theoryPredictionsFor(expResult, smsTopList, maxMassDist=0.2, useBestDataset=True):
    """
    Compute theory predictions for the given experimental result, using the list of elements
    in smsTopList.
    For each Txname appearing in expResult, it collects the elements and efficiencies, 
    combine the masses (if needed) and compute the conditions (if existing).
    
    :parameter expResult: expResult to be considered (ExpResult object)
    :parameter smsTopList: list of topologies containing elements (TopologyList object)
    :parameter maxMassDist: maximum mass distance for clustering elements (float)
    :parameter useBestDataset: If True, uses only the best dataset (signal region).
               If False, returns predictions for all datasets.
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
    if useBestDataset:
        bestResults = _getBestResults(dataSetResults)
        bestResults.expResult = expResult
        for theoPred in bestResults:
            theoPred.expResult = expResult    
        return bestResults
    else:        
        allResults = sum(dataSetResults)
        for theoPred in allResults: theoPred.expResult = expResult       
        return allResults

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
    bestXsec = 0.*fb
    for predList in dataSetResults:
        if len(predList) != 1:
            logger.error("Multiple clusters should only exist for upper limit results!")
            raise SModelSError()
        dataset = predList[0].dataset
        if dataset.dataInfo.dataType != 'efficiencyMap':
            logger.error("Multiple data sets should only exist for efficiency map results!")
            raise SModelSError()                    
        pred = predList[0]
        xsec = pred.xsection        
        expectedR = ( xsec.value/dataset.getSRUpperLimit(0.05,True,False) ).asNumber()
        if expectedR > bestExpectedR or (expectedR == bestExpectedR and xsec.value > bestXsec):
            bestExpectedR = expectedR
            bestPredList = predList
            bestXsec = xsec.value
    
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
    # Select elements belonging to expResult and apply efficiencies
    elements = _getElementsFrom(smsTopList, dataset)    
    
    #Check dataset sqrts format:
    if (dataset.globalInfo.sqrts/TeV).normalize()._unit:
            ID = dataset.globalInfo.id
            logger.error("Sqrts defined with wrong units for %s" %(ID) )
            return False
            
    #Remove unwanted cross sections
    newelements = []
    for el in elements:
        el.weight = el.weight.getXsecsFor(dataset.globalInfo.sqrts)
        if not el.weight: continue
        newelements.append(el)
    elements = newelements
    if len(elements) == 0: return None
       

    # Combine elements according to their respective constraints and masses
    # (For efficiencyMap analysis group all elements)
    clusters = _combineElements(elements, dataset, maxDist=maxMassDist)


    # Collect results and evaluate conditions    
    for cluster in clusters:
        theoryPrediction = TheoryPrediction()
        theoryPrediction.dataset = dataset
        theoryPrediction.txnames = cluster.txnames  
        theoryPrediction.xsection = _evalConstraint(cluster)
        theoryPrediction.conditions = _evalConditions(cluster)
        theoryPrediction.cluster = cluster
        theoryPrediction.mass = cluster.getAvgMass()
        theoryPrediction.PIDs = cluster.getPIDs()
        theoryPrediction.IDs = cluster.getIDs()
        predictionList._theoryPredictions.append(theoryPrediction)        
        

    if len(predictionList) == 0: return None
    else: return predictionList

def _getElementsFrom(smsTopList, dataset):
    """
    Get elements that belong to any of the TxNames in dataset 
    (appear in any of constraints in the result).    
    Loop over all elements in smsTopList and returns a copy of the elements belonging
    to any of the constraints (i.e. have efficiency != 0). The copied elements
    have their weights multiplied by their respective efficiencies.
    
    :parameter dataset:  Data Set to be considered (DataSet object)
    :parameter smsTopList: list of topologies containing elements (TopologyList object)
    :returns: list of elements (Element objects)    
    """
    elements = []  
    newMissingelements = []
    for txname in dataset.txnameList:
        for top in smsTopList:      
            allNewElements = []
            hasTop = txname._topologyList.hasTopology(top)
            if not hasTop: continue                 
            for i,el in enumerate(top.getElements()):       
                txnameEl = txname.hasElementAs(el)  #Check if element appears in txname                
                if not txnameEl: continue  
                newElements, factors = txname.getNewElementsAndFactors(txnameEl)

                for i,newEl in enumerate(newElements):                                
                    if not factors[i]: continue                    
                    newElement, eff = txname.checkDecayTypes(newEl)    
                    if eff: 
                        newElement.tested = True                                        
                        newElement.txname = txname                         
                        newElement.weight *= factors[i]*eff                                     
                        newElement.eff = factors[i]*eff  
                    else: newElement.weight *= factors[i]     
      
                    if newElement.covered:
                        el.covered = True     
                                
                    if newElement.tested:                                               
                        el.tested = True
                        elements.append(newElement) #Save element with correct branch ordering that are tested by experiment
                    else: 
                        if not newElement in newMissingelements: 
                            newMissingelements.append(newElement)

    # collect new elements that were reweighted and them  
    newMissingelements = [missingEl for missingEl in newMissingelements if not missingEl in elements]
    smsTopList.extraMissingElements = newMissingelements  
                      
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
    
    if dataset.dataInfo.dataType == 'efficiencyMap':        
        cluster = clusterTools.groupAll(elements)  
        clusters.append(cluster)
    elif dataset.dataInfo.dataType == 'upperLimit':
        txnames = list(set([el.txname for el in elements]))        
        for txname in txnames:
            txnameEls = []
            for element in elements:
                if not element.txname == txname:
                    continue
                else: txnameEls.append(element)
                    
            txnameClusters = clusterTools.clusterElements(txnameEls, maxDist)         
            clusters += txnameClusters
    else:
        logger.warning("Unkown data type: %s. Data will be ignored." 
                       % dataset.dataInfo.dataType)

    return clusters


def _evalConstraint(cluster):
    """
    Evaluate the constraint inside an element cluster.      
    If the cluster refers to a specific TxName, sum all the elements' weights
    according to the analysis constraint.
    For efficiency map results, sum all the elements' weights.
    
    :parameter cluster: cluster of elements (ElementCluster object)
    :returns: cluster cross section
    """

    if cluster.getDataType() == 'efficiencyMap':
        return cluster.getTotalXSec()
    elif cluster.getDataType() == 'upperLimit':
        if len(cluster.txnames) != 1:
            logger.error("An upper limit cluster should never contain more than one TxName")
            raise SModelSError()
        txname = cluster.txnames[0]
        if not txname.constraint or txname.constraint == "not yet assigned":
            return txname.constraint
        exprvalue = _evalExpression(txname.constraint,cluster)

        return exprvalue
    else:
        logger.error("Unknown data type %s" %(str(cluster.getDataType())))
        raise SModelSError()
    

def _evalConditions(cluster):
    """
    Evaluate the conditions (if any) inside an element cluster.
    
    :parameter cluster: cluster of elements (ElementCluster object)    
    :returns: list of condition values (floats) if analysis type == upper limit. None, otherwise.    
    """
    conditionVals = {}
    for txname in cluster.txnames:
        if not txname.condition or txname.condition == "not yet assigned":
            continue        
        #Make sure conditions is always a list
        if isinstance(txname.condition,str):
            conditions =  [txname.condition]
        elif isinstance(txname.condition,list):
            conditions = txname.condition
        else:
            logger.error("Conditions should be a list or a string")
            raise SModelSError()
            
        # Loop over conditions
               
        for cond in conditions:
            exprvalue = _evalExpression(cond,cluster)            
            if type(exprvalue) == type(crossSection.XSection()):
                conditionVals[cond] = exprvalue.value
            else:
                conditions[cond] = exprvalue
    
    if not conditionVals:
        return None
    else:            
        return conditionVals    
        
        
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
    :returns: xsection for the expression. Can be a XSection object, a float or not numerical (None,string,...)
    
    """
#Get cross section info from cluster (to generate zero cross section values):
    infoList = cluster.elements[0].weight.getInfo()    

#Generate elements appearing in the string expression with zero cross sections:
 
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
        return exprvalue[0] #Return XSection object
    
    return exprvalue
