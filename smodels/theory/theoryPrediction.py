"""
.. module:: theoryPrediction
   :synopsis: Provides a class to encapsulate the results of the computation of
              reference cross sections and related functions.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. autofunction:: _getElementsFrom
"""

from smodels.theory import clusterTools, crossSection, element
from smodels.theory.particleNames import elementsInStr
from smodels.theory.auxiliaryFunctions import cSim, cGtr # DO NOT REMOVE
import copy
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

    def dataType( self ):
        """
        Return the type of dataset
        """
                
        return self.dataset.getType()
    

    def getUpperLimit(self, expected=False ):
        """
        Get the upper limit on sigma*eff
        :param expected: return expected Upper Limit, instead of observed.
        """
               
        if expected:
            if not hasattr(self,'expectedUL'):
                self.expectedUL = self.expResult.getUpperLimitFor(mass=self.mass, \
                           dataID=self.dataset, expected=expected, txname = self.txnames[0])
            return self.expectedUL
        else:
            if not hasattr(self,'upperLimit'):
                self.upperLimit = self.expResult.getUpperLimitFor(mass=self.mass,
                                                                  dataID=self.dataId(),
                                                                  expected=expected,
                                                                  txname = self.txnames[0])
                return self.upperLimit

    def getRValue(self, expected = False ):
        """
        Get the r value = theory prediction / experimental upper limit
        """
        
        if expected:
            eul=self.getUpperLimit(expected=True)
            if type(eul)==type(None) or eul.asNumber(fb)==0.:
                
                return None
            else:
                return self.xsection.value/eul

        ul = self.getUpperLimit()
        if type(ul) == type(None) or ul.asNumber(fb) == 0.:
            return None
        
        return self.xsection.value/ul

    def computeStatistics(self,marginalize=False):
        """
        Compute the likelihood, chi2 and expected upper limit for this theory prediction.
        The resulting values are stored as the likelihood, chi2 and expectedUL
        attributes.
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        """
        if type ( self.dataset ) == list:
            ## a prediction for a combined result? special
            lumi = self.expResult.globalInfo.lumi
            pred = (self.xsection.value*lumi).asNumber() / self.effectiveEff
            nsig = [ pred * x for x in self.efficiencies ]
            self.likelihood = self.expResult.combinedLikelihood ( nsig, marginalize=marginalize )
            self.chi2 = self.expResult.totalChi2 ( nsig, marginalize=marginalize )
            # self.expectedUL = None
            return

        if self.dataset.dataInfo.dataType == 'upperLimit':
            ## FIXME need to write procedure for combined results!!
            self.likelihood = None
            self.chi2 = None
            self.expectedUL = None
            return

        lumi = self.dataset.globalInfo.lumi
        nsig = (self.xsection.value*lumi).asNumber()
        llhd = self.dataset.likelihood(nsig,marginalize=marginalize)
        chi2 = self.dataset.chi2(nsig,marginalize=marginalize)
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

    def isCombined ( self ):
        """ does this theory pred stem from combining several datasets?
            (implies the existence of a covariance matrix) """
        return type(self.dataset)==list

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
            dataId = self.dataset.dataInfo.dataId
            folderName = self.dataset.dataInfo.path
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

def theoryPredictionsFor( expResult, smsTopList, maxMassDist=0.2,
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
    :parameter marginalize: If true, marginalize nuisances. Ff false, profile them.
    :returns:  a TheoryPredictionList object containing a list of TheoryPrediction
               objects
    """

    preds = _sortPredictions( expResult, smsTopList, maxMassDist, combinedResults)
    for preds in preds.values():
        effs = [ pred.effectiveEff for pred in preds ]
        if sum( effs ) == 0.:
            logger.info( "all efficiencies of combination in %s are zero. will skip." % expResult.globalInfo.id )
            break

    dataSetResults = []
    #Compute predictions for each data set (for UL analyses there is one single set)
    for dataset in expResult.datasets:
        predList = _getDataSetPredictions(dataset,smsTopList,maxMassDist)
        if predList:
            dataSetResults.append(predList)
    if not dataSetResults:
        return None

    #For results with more than one dataset, return all dataset predictions
    #if useBestDataSet=False and combinedResults=False:
    if not useBestDataset and not combinedResults:
        allResults = sum(dataSetResults)
        for theoPred in allResults:
            theoPred.expResult = expResult
        return allResults
    
    #Else include best signal region results
    bestResults = _getBestResults(dataSetResults)
    #If combinedResults = True, also include the combined
    #result (when available):
    if combinedResults:
        combinedDataSetResult = _getCombinedResultFor(dataSetResults,expResult)
        if combinedDataSetResult:
            bestResults.append(combinedDataSetResult)

    for theoPred in bestResults:
        theoPred.expResult = expResult
        
    return bestResults


def _getCombinedResultsFor(dataSetResults,expResult):
    """
    Compute the compbined result for all datasets, if covariance
    matrices are available. Return a TheoryPrediction object
    with the signal cross-section summed over all the signal regions
    and the respective upper limit.
    """

def _mergePredictions ( preds, combinedUL, combinedEUL ):
    """ merge theory predictions, for the combined prediction. """
    if len(preds) == 0: return None
    ret=copy.deepcopy( preds[0] )
    ret.efficiencies = []
    eff, wtot = 0., 0.
    for pred in preds:
        w = pred.xsection.value.asNumber(fb)
        eff += pred.effectiveEff * w
        ret.efficiencies.append ( pred.effectiveEff )
        wtot += w
    eff = eff / wtot
    # print ( "combinedUL=",combinedUL )
    ret.xsection.value = ret.xsection.value / ret.effectiveEff * eff ## / preds[0].effectiveEff
    ret.combinedUL = None
    if type(combinedUL) != type(None):
        ret.combinedUL = combinedUL * eff
    ret.combinedExpectedUL = None
    if combinedEUL is not None:
        ret.combinedExpectedUL = combinedEUL * eff
    ret.effectiveEff = eff
    # ret.dataset = FIXME special
    ret.dataset = [ x.dataset for x in preds ] # we collect all datasets
    return ret

def _sortPredictions ( expResult, smsTopList, maxMassDist, combine ):
    """ returns dictionary of predictions, sorted by XSectionInfo.
        if combine is false, return empty dict.
    """
    preds={}
    if not hasattr ( expResult.globalInfo, "covariance" ) or \
       not hasattr ( expResult.globalInfo, "datasetOrder" ) or \
       not combine:
           return preds
    dsOrder = expResult.globalInfo.datasetOrder
    if type ( dsOrder ) == str:
        ## for debugging only, we allow a single dataset
        dsOrder = [ dsOrder ]
    for dsname in dsOrder:
        dataset=expResult.getDataset ( dsname )
        if dataset == None:
            txt = "In %s: dataset %s does not exist." % \
                  ( expResult.globalInfo.id, dsname )
            raise SModelSError ( txt )
        predList = _getDataSetPredictions(dataset,smsTopList,maxMassDist,True)
        if predList:
            for pred in predList:
                info = pred.xsection.info
                if not info in preds.keys():
                    preds[info]=[]
                preds[info].append ( pred )
        else:
            pass
            # logger.error ( "this is the culprit. we have no predlist. but we need one for combination. lets make one artificially." )
    return preds

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
        expectedR = (xsec.value/dataset.getSRUpperLimit(0.05,True,False) ).asNumber()
        if expectedR > bestExpectedR or (expectedR == bestExpectedR and xsec.value > bestXsec):
            bestExpectedR = expectedR
            bestPredList = predList
            bestXsec = xsec.value

    return bestPredList

def _getDataSetPredictions(dataset,smsTopList,maxMassDist,force_creation=False):
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
    :parameter force_creation: force creation of a prediction, even if efficiency
                               is zero. We need this to have consistent input for
                               combinations.
    :returns:  a TheoryPredictionList object containing a list of TheoryPrediction 
               objects
    """

    predictionList = TheoryPredictionList()
    # Select elements belonging to expResult and apply efficiencies
    elements = _getElementsFrom(smsTopList, dataset,force_creation)

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
        theoryPrediction.upperLimit = theoryPrediction.getUpperLimit()
        predictionList._theoryPredictions.append(theoryPrediction)

    if len(predictionList) == 0: return None
    else: return predictionList

def _getElementsFrom(smsTopList, dataset, force_creation):
    """
    Get elements that belong to any of the TxNames in dataset
    (appear in any of constraints in the result).
    Loop over all elements in smsTopList and returns a copy of the elements belonging
    to any of the constraints (i.e. have efficiency != 0). The copied elements
    have their weights multiplied by their respective efficiencies.

    :parameter dataset:  Data Set to be considered (DataSet object)
    :parameter smsTopList: list of topologies containing elements 
                           (TopologyList object)
    :parameter force_creation: force creation of element, even if eff == 0. 
                               This is needed for combined results.
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
                if not eff and not force_creation: continue
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

    exprvalue = eval(expr,globals(),
                     {'cGtr' : cGtr, 'cSim' : cSim,'Cgtr' : cGtr, 'Csim' : cSim})
    if type(exprvalue) == type(crossSection.XSectionList()):
        if len(exprvalue) != 1:
            logger.error("Evaluation of expression "+expr+" returned multiple values.")
        return exprvalue[0] #Return XSection object
    return exprvalue
