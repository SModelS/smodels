"""
.. module:: theoryPrediction
   :synopsis: Provides a class to encapsulate the results of the computation of
              reference cross sections and related functions.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. autofunction:: _getElementsFrom
"""

from smodels.theory import clusterTools, crossSection, element
from smodels.theory.particleNames import elementsInStr
from smodels.theory.auxiliaryFunctions import cSim, cGtr
from smodels.tools.physicsUnits import TeV,fb
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.experiment.datasetObj import CombinedDataSet
from smodels.tools.smodelsLogging import logger

class TheoryPrediction(object):
    """
    An instance of this class represents the results of the theory prediction
    for an analysis.

    :ivar analysis: holds the analysis (ULanalysis or EManalysis object)
                    to which the prediction refers
    :ivar xsection: xsection of the theory prediction
                (relevant cross section to be compared with the experimental limits).
                For EM-type analyses, it corresponds to sigma*eff, for
                UL-type analyses, eff is considered to be 1.
                It is a XSection object.
    :ivar effectiveEff: the effective efficiency for EM-type analyses.
    :ivar conditions: list of values for the analysis conditions
                      (only for upper limit-type analysis, e.g. analysis=ULanalysis)
    :ivar mass: mass of the cluster to which the theory prediction refers
                (only for upper limit-type analysis, e.g. analysis=ULanalysis)

    """
    def __init__(self):
        self.xsection = None
        self.effectiveEff = None
        self.conditions = None
        self.mass = None

    def dataId(self):
        """
        Return ID of dataset
        """
        
        return self.dataset.getID()
    
    def analysisId(self):
        """
        Return experimental analysis ID
        """
        
        return self.dataset.globalInfo.id

    def dataType(self):
        """
        Return the type of dataset
        """
                
        return self.dataset.getType()
    

    def getUpperLimit(self, expected=False ):
        """
        Get the upper limit on sigma*eff.
        For UL-type results, use the UL map. For EM-Type returns
        the corresponding dataset (signal region) upper limit.
        For combined results, returns the upper limit on the
        total sigma*eff (for all signal regions/datasets). 
        
        :param expected: return expected Upper Limit, instead of observed.
        
        :return: upper limit (Unum object)
        """
        
        #First check if the upper-limit and expected upper-limit have already been computed.
        #If not, compute it and store them.
        if not hasattr(self, 'expectedUL') or not hasattr(self, 'upperLimit'):
            
            if self.dataType() == 'efficiencyMap':
                self.expectedUL = self.dataset.getSRUpperLimit(expected=True)
                self.upperLimit = self.dataset.getSRUpperLimit(expected=False)
            if self.dataType() == 'upperLimit':                
                self.expectedUL = self.dataset.getUpperLimitFor(mass=self.mass,
                                                                txnames=self.txnames,
                                                                expected=True)
                self.upperLimit = self.dataset.getUpperLimitFor(mass=self.mass,
                                                                txnames=self.txnames,
                                                                expected=False)
            if self.dataType() == 'combined':
                lumi = self.expResult.globalInfo.lumi
                #Create a list of signal events in each dataset/SR sorted according to datasetOrder
                srNsigDict = dict([[pred.dataset.getID(),(pred.xsection.value*lumi).asNumber()] for pred in self.datasetPredictions])
                srNsigs = [srNsigDict[dataID] if dataID in srNsigDict else 0. for dataID in self.dataset.globalInfo.datasetOrder]
                self.expectedUL = self.dataset.getCombinedUpperLimitFor(srNsigs,expected=True)
                self.upperLimit = self.dataset.getCombinedUpperLimitFor(srNsigs,expected=False)
                          
        #Return the expected or observed UL:
        if expected:
            return self.expectedUL
        else:
            return self.upperLimit

    def getRValue(self, expected = False):
        """
        Get the r value = theory prediction / experimental upper limit
        """
        
        upperLimit = self.getUpperLimit(expected)
        if upperLimit is None or upperLimit.asNumber(fb)==0.:                
            return None
        else:
            return (self.xsection.value/upperLimit).asNumber()

    def computeStatistics(self,marginalize=False):
        """
        Compute the likelihood, chi2 and expected upper limit for this theory prediction.
        The resulting values are stored as the likelihood and chi2
        attributes.
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        """
        
        
        if self.dataType()  == 'upperLimit':
            self.likelihood = None
            self.chi2 = None

        elif self.dataType() == 'efficiencyMap':            
            lumi = self.dataset.globalInfo.lumi
            nsig = (self.xsection.value*lumi).asNumber()
            llhd = self.dataset.likelihood(nsig,marginalize=marginalize)
            chi2 = self.dataset.chi2(nsig,marginalize=marginalize)    
            self.likelihood =  llhd
            self.chi2 =  chi2
            
        elif self.dataType() == 'combined':
            lumi = self.expResult.globalInfo.lumi
            #Create a list of signal events in each dataset/SR sorted according to datasetOrder
            srNsigDict = dict([[pred.dataset.getID(),(pred.xsection.value*lumi).asNumber()] for pred in self.datasetPredictions])
            srNsigs = [srNsigDict[dataID] if dataID in srNsigDict else 0. for dataID in self.dataset.globalInfo.datasetOrder]                
            self.likelihood = self.dataset.combinedLikelihood(srNsigs, marginalize=marginalize)
            self.chi2 = self.dataset.totalChi2(srNsigs, marginalize=marginalize)
            

    def getmaxCondition(self):
        """
        Returns the maximum xsection from the list conditions

        :returns: maximum condition xsection (float)
        """

        if not self.conditions: return 0.
        # maxcond = 0.
        values = [ 0. ]
        for value in self.conditions.values():
            if value == 'N/A': return value
            if value == None: continue
            #print ( "value=",value,type(value),float(value) )
            #maxcond = max(maxcond,float(value))
            values.append ( float(value) )
        return max(values)
        # return maxcond

    def __str__(self):
        ret = "%s:%s" % ( self.analysisId(), self.xsection )
        #self.computeStatistics()
        #ret += ", llhd=%f." % self.likelihood
        return ret

    def describe ( self ):
        if not hasattr ( self, "chi2" ):
            self.computeStatistics()
        # return a lengthy description
        ret =  "[theoryPrediction] analysis: %s\n" % self.analysisId()
        ret += "     prediction (sigma*eff): %s\n" % self.xsection
        ret += "         prediction (sigma): %s\n" % ( self.xsection.value / self.effectiveEff )
        ret += "       effective efficiency: %s\n" % self.effectiveEff
        ds = "None"
        if type (self.dataset) == list:
            ds = "multiple (%d combined)" % len(self.dataset)
        else:
            dataId = self.dataset.getID()
            ds = "%s (%s)" % ( dataId, self.dataset.folderName() )
        ret += "                   datasets: %s\n" % ds
        ret += "      obs limit (sigma*eff): %s\n" % self.getUpperLimit()
        ret += "      exp limit (sigma*eff): %s\n" % self.getUpperLimit( expected=True )
        ret += "          obs limit (sigma): %s\n" % (self.getUpperLimit() / self.effectiveEff )
        eul_ = self.getUpperLimit( expected=True )
        eul = eul_
        if type(eul)!=type(None):
            eul = eul / self.effectiveEff
         
        ret += "          exp limit (sigma): %s\n" % ( eul )
        ret += "                      obs r: %f\n" % ( self.xsection.value / self.getUpperLimit() )
        exp_r = None
        if type(eul_)!=type(None):
                exp_r = self.xsection.value / eul_
        ret += "                      exp r: %f\n" % ( exp_r )
        ret += "                       chi2: %s\n" % ( self.chi2 )
        return ret

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

    def append(self,theoryPred):
        self._theoryPredictions.append ( theoryPred )

    def __str__(self):
        if len ( self._theoryPredictions ) == 0:
            return "no predictions."
        ret = "%d predictions: " % len ( self._theoryPredictions )
        ret += ", ".join ( [ str(s) for s in self._theoryPredictions ] )
        return ret

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

def theoryPredictionsFor(expResult, smsTopList, maxMassDist=0.2,
                useBestDataset=True, combinedResults=True, marginalize=False ):
    """
    Compute theory predictions for the given experimental result, using the list of
    elements in smsTopList.
    For each Txname appearing in expResult, it collects the elements and
    efficiencies, combine the masses (if needed) and compute the conditions
    (if exist).

    :parameter expResult: expResult to be considered (ExpResult object)
    :parameter smsTopList: list of topologies containing elements
                           (TopologyList object)
    :parameter maxMassDist: maximum mass distance for clustering elements (float)
    :parameter useBestDataset: If True, uses only the best dataset (signal region).
               If False, returns predictions for all datasets.
    :parameter combinedResults: add theory predictions that result from
               combining datasets.
    :parameter marginalize: If true, marginalize nuisances. If false, profile them.
    :returns:  a TheoryPredictionList object containing a list of TheoryPrediction
               objects
    """

    dataSetResults = []
    #Compute predictions for each data set (for UL analyses there is one single set)
    for dataset in expResult.datasets:
        predList = _getDataSetPredictions(dataset,smsTopList,maxMassDist)
        if predList:
            dataSetResults.append(predList)
    if not dataSetResults:
        return None
    elif len(dataSetResults) == 1:
        result = dataSetResults[0]
        for theoPred in result:
            theoPred.expResult = expResult
            theoPred.upperLimit = theoPred.getUpperLimit()
        return result

    #For results with more than one dataset, return all dataset predictions
    #if useBestDataSet=False and combinedResults=False:
    if not useBestDataset and not combinedResults:
        allResults = sum(dataSetResults)
        for theoPred in allResults:
            theoPred.expResult = expResult
            theoPred.upperLimit = theoPred.getUpperLimit()
        return allResults
    
    #Else include best signal region results
    bestResults = TheoryPredictionList()
    bestResults.append(_getBestResult(dataSetResults))
    #If combinedResults = True, also include the combined result (when available):
    if combinedResults and len(dataSetResults) > 1:
        combinedDataSetResult = _getCombinedResultFor(dataSetResults,
                                                      expResult,marginalize)
        if combinedDataSetResult:
            bestResults.append(combinedDataSetResult)

    for theoPred in bestResults:
        theoPred.expResult = expResult
        theoPred.upperLimit = theoPred.getUpperLimit()
        
    return bestResults

def _getCombinedResultFor(dataSetResults,expResult,marginalize=False):
    """
    Compute the compbined result for all datasets, if covariance
    matrices are available. Return a TheoryPrediction object
    with the signal cross-section summed over all the signal regions
    and the respective upper limit.
    
    :param datasetPredictions: List of TheoryPrediction objects for each signal region
    :param expResult: ExpResult object corresponding to the experimental result
    :parameter marginalize: If true, marginalize nuisances. If false, profile them.    
    
    :return: TheoryPrediction object
    """
        
    if len(dataSetResults) == 1:
        return dataSetResults[0]
    elif not expResult.hasCovarianceMatrix():
        return None
    
    txnameList = []
    elementList = []
    totalXsec = None
    massList = []
    PIDList = []
    IDList = []
    datasetPredictions = []
    for predList in dataSetResults:
        if len(predList) != 1:
            raise SModelSError("Results with multiple datasets should have a single theory prediction (EM-type)!")
        pred = predList[0]
        datasetPredictions.append(pred)
        txnameList += pred.txnames
        elementList += pred.elements
        if not totalXsec:
            totalXsec = pred.xsection
        else:
            totalXsec += pred.xsection
        massList.append(pred.mass)
        for pidEntry in pred.PIDs:
            if not pidEntry in PIDList:
                PIDList.append(pidEntry)
        IDList += pred.IDs
        
    txnameList = list(set(txnameList))
    massList = list(set(massList))
    IDList = list(set(IDList))
    if len(massList) > 1:
        mass = None
    else:
        mass = massList[0]
    
    
    #Create a combinedDataSet object:
    combinedDataset = CombinedDataSet(expResult)
    combinedDataset._marginalize = marginalize
    #Create a theory preidction object for the combined datasets:
    theoryPrediction = TheoryPrediction()
    theoryPrediction.dataset = combinedDataset
    theoryPrediction.txnames = txnameList
    theoryPrediction.xsection = totalXsec
    theoryPrediction.datasetPredictions = datasetPredictions
    theoryPrediction.conditions = None
    theoryPrediction.elements = elementList
    theoryPrediction.mass = mass
    theoryPrediction.PIDs = PIDList
    theoryPrediction.IDs = IDList
    
    
    return theoryPrediction
    
def _getBestResult(dataSetResults):
    """
    Returns the best result according to the expected upper limit

    :param datasetPredictions: list of TheoryPredictionList objects
    :return: best result (TheoryPrediction object)
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
        if dataset.getType() != 'efficiencyMap':
            logger.error("Multiple data sets should only exist for efficiency map results!")
            raise SModelSError()
        pred = predList[0]
        xsec = pred.xsection
        expectedR = (xsec.value/dataset.getSRUpperLimit(0.05,True,False) ).asNumber()
        if expectedR > bestExpectedR or (expectedR == bestExpectedR and xsec.value > bestXsec):
            bestExpectedR = expectedR
            bestPred = pred
            bestXsec = xsec.value

    return bestPred

def _getDataSetPredictions(dataset,smsTopList,maxMassDist):
    """
    Compute theory predictions for a given data set.
    For upper-limit results returns the list of theory predictions for the
    experimental result.  For efficiency-map results returns the list of theory
    predictions for the signal region.  Uses the list of elements in
    smsTopList.
    For each Txname appearing in dataset, it collects the elements and efficiencies,
    combine the masses (if needed) and compute the conditions (if existing).

    :parameter dataset: Data Set to be considered (DataSet object)
    :parameter smsTopList: list of topologies containing elements 
                           (TopologyList object)
    :parameter maxMassDist: maximum mass distance for clustering elements (float)
    :returns:  a TheoryPredictionList object containing a list of TheoryPrediction 
               objects
    """

    predictionList = TheoryPredictionList()
    # Select elements belonging to expResult and apply efficiencies
    elements = _getElementsFrom(smsTopList, dataset)

    #Check dataset sqrts format:
    if (dataset.globalInfo.sqrts/TeV).normalize()._unit:
            ID = dataset.globalInfo.id
            logger.error( "Sqrt(s) defined with wrong units for %s" % (ID) )
            return False

    #Remove unwanted cross sections
    newelements = []
    for el in elements:
        el.weight = el.weight.getXsecsFor(dataset.globalInfo.sqrts)
        if not el.weight: continue
        newelements.append(el)
    elements = newelements
    if len(elements) == 0:
        return None

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
        theoryPrediction.elements = cluster.elements
        theoryPrediction.mass = cluster.getAvgMass()
        theoryPrediction.PIDs = cluster.getPIDs()
        theoryPrediction.IDs = cluster.getIDs()
        predictionList._theoryPredictions.append(theoryPrediction)

    if len(predictionList) == 0:
        return None
    else:
        return predictionList

def _getElementsFrom(smsTopList, dataset):
    """
    Get elements that belong to any of the TxNames in dataset
    (appear in any of constraints in the result).
    Loop over all elements in smsTopList and returns a copy of the elements belonging
    to any of the constraints (i.e. have efficiency != 0). The copied elements
    have their weights multiplied by their respective efficiencies.

    :parameter dataset:  Data Set to be considered (DataSet object)
    :parameter smsTopList: list of topologies containing elements 
                           (TopologyList object)
    :returns: list of elements (Element objects)
    """

    elements = []
    for txname in dataset.txnameList:
        for top in smsTopList:
            itop = txname._topologyList.index(top)  #Check if the topology appear in txname
            if itop is None: continue
            for el in top.getElements():
                newEl = txname.hasElementAs(el)  #Check if element appears in txname
                if not newEl: continue
                el.covered = True
                eff = txname.getEfficiencyFor(newEl.getMasses())
                if not eff:
                    continue
                el.tested = True
                newEl.eff = eff
                newEl.weight *= eff
                newEl.txname = txname
                elements.append(newEl) #Save element with correct branch ordering

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

    if dataset.getType() == 'efficiencyMap':
        cluster = clusterTools.groupAll(elements)
        clusters.append(cluster)
    elif dataset.getType() == 'upperLimit':
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
                       % dataset.getType())

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

#Get weights for elements appearing in stringExpr
    weightsDict = {}
    evalExpr = stringExpr.replace("'","").replace(" ","")
    for i,elStr in enumerate(elementsInStr(evalExpr)):
        el = element.Element(elStr)
        weightsDict['w%i'%i] = crossSection.XSectionList(infoList)
        for el1 in cluster.elements:
            if el1.particlesMatch(el):
                weightsDict['w%i'%i].combineWith(el1.weight)
                el.combineMotherElements(el1)
        evalExpr = evalExpr.replace(elStr,'w%i'%i)

    weightsDict.update({"Cgtr" : cGtr, "cGtr" : cGtr, "cSim" : cSim, "Csim" : cSim})
    exprvalue = eval(evalExpr, weightsDict)
    if type(exprvalue) == type(crossSection.XSectionList()):
        if len(exprvalue) != 1:
            logger.error("Evaluation of expression "+evalExpr+" returned multiple values.")
        return exprvalue[0] #Return XSection object
    return exprvalue
