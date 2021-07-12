"""
.. module:: datasetObj
   :synopsis: Holds the classes and methods used to read and store the information in the
              data folders.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


import os,glob,json
import numpy as np
from smodels.experiment import txnameObj,infoObj
from smodels.tools.physicsUnits import fb
from smodels.tools.simplifiedLikelihoods import LikelihoodComputer, Data, UpperLimitComputer
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from smodels.theory.auxiliaryFunctions import getAttributesFrom,getValuesForObj
from smodels.tools.smodelsLogging import logger
from smodels.theory.auxiliaryFunctions import elementsInStr
from smodels.theory.element import Element

import itertools

# if on, will check for overlapping constraints
_complainAboutOverlappingConstraints = True

class DataSet(object):
    """
    Holds the information to a data set folder (TxName objects, dataInfo,...)
    """

    def __init__(self, path=None, info=None, createInfo=True,
                    discard_zeroes=True, databaseParticles = None):
        """
        :param discard_zeroes: discard txnames with zero-only results
        """

        self.path = path
        self.globalInfo = info
        self.txnameList = []

        if path and createInfo:
            logger.debug('Creating object based on data folder : %s' %self.path)

            #Get data folder info:
            if not os.path.isfile(os.path.join(path,"dataInfo.txt")):
                logger.error("dataInfo.txt file not found in " + path)
                raise TypeError
            self.dataInfo = infoObj.Info(os.path.join(path,"dataInfo.txt"))

            #Get list of TxName objects:
            for txtfile in glob.iglob(os.path.join(path,"*.txt")):
                try:
                    txname = txnameObj.TxName(txtfile,self.globalInfo,
                                            self.dataInfo, databaseParticles)
                    if discard_zeroes and txname.hasOnlyZeroes():
                        logger.debug ( "%s, %s has only zeroes. discard it." % \
                                         ( self.path, txname.txName ) )
                        continue
                    self.txnameList.append(txname)
                except TypeError: continue

            self.txnameList.sort()
            self.checkForRedundancy(databaseParticles)

    def isCombinableWith( self, other ):
        """
        Function that reports if two datasets are mutually uncorrelated = combinable.

        :param other: datasetObj to compare self with
        """
        id1, id2 = self.globalInfo.id, other.globalInfo.id
        if id1 == id2: ## we are always correlated with ourselves
            return False
        from smodels.tools.physicsUnits import TeV
        ds = abs (self.globalInfo.sqrts.asNumber(TeV) - other.globalInfo.sqrts.asNumber(TeV) )
        if ds > 1e-5: ## not the same
            return True
        def getCollaboration ( ds ):
            return "CMS" if "CMS" in ds.globalInfo.id else "ATLAS"
        coll1, coll2 = getCollaboration ( self ), getCollaboration ( other )
        if coll1 != coll2:
            return True

        if self.isGlobalFieldCombinableWith_ ( other ):
            return True
        if other.isGlobalFieldCombinableWith_ ( self ):
            return True
        if self.isLocalFieldCombinableWith_ ( other ):
            return True
        if other.isLocalFieldCombinableWith_ ( self ):
            return True
        if self.isCombMatrixCombinableWith_ ( other ):
            return True
        return False

    def isCombMatrixCombinableWith_ ( self, other ):
        """ check for combinability via the combinations matrix """
        if not hasattr ( self.globalInfo, "_combinationsmatrix" ):
            return False
        if self.globalInfo._combinationsmatrix == None:
            return False
        idSelf = self.globalInfo.id
        didSelf = self.dataInfo.dataId
        selflabel = f"{idSelf}:{didSelf}"
        idOther = other.globalInfo.id
        didOther = other.dataInfo.dataId
        otherlabel = f"{idOther}:{didOther}"
        for label, combs in self.globalInfo._combinationsmatrix.items():
            if label in [ idSelf, didSelf ]:
                ## match! with self! is "other" in combs?
                if idOther in combs or otherlabel in combs:
                    return True
            if label in [ idOther, didOther ]:
                ## match! with other! is "self" in combs?
                if idSelf in combs or selflabel in combs:
                    return True
        return False

    def isGlobalFieldCombinableWith_ ( self, other ):
        """ check for 'combinableWith' fields in globalInfo, check if <other> matches.
        this check is at analysis level (not at dataset level).

        :params other: a dataset to check against
        :returns: true, if pair is marked as combinable, else false
        """
        if not hasattr ( self.globalInfo, "combinableWith" ):
            return False
        tokens = self.globalInfo.combinableWith.split ( "," )
        idOther = other.globalInfo.id
        for t in tokens:
            if ":" in t:
                logger.error ( "combinableWith field in globalInfo is at the analysis level. You specified a dataset-level combination %s." % t )
                raise SModelSError()
        if idOther in tokens:
            return True
        return False

    def isLocalFieldCombinableWith_ ( self, other ):
        """ check for 'combinableWith' fields in globalInfo, check if <other> matches.
        this check is at dataset level (not at dataset level).

        :params other: a dataset to check against
        :returns: true, if pair is marked as combinable, else false
        """
        if not hasattr ( self.dataInfo, "combinableWith" ):
            return False
        tokens = self.dataInfo.combinableWith.split ( "," )
        for t in tokens:
            if not ":" in t:
                logger.error ( "combinableWith field in dataInfo is at the dataset level. You specified an analysis-level combination %s." % t )
                raise SModelSError()
        idOther = other.globalInfo.id
        didOther = other.dataInfo.dataId
        label = f"{idOther}:{didOther}"
        if label in tokens:
            return True
        return False

    def checkForRedundancy(self,databaseParticles):
        """ In case of efficiency maps, check if any txnames have overlapping
            constraints. This would result in double counting, so we dont
            allow it. """
        if self.getType() == "upperLimit":
            return False
        logger.debug ( "checking for redundancy" )
        datasetElements = []
        for tx in self.txnameList:
            if hasattr(tx, 'finalState'):
                finalState = tx.finalState
            else:
                finalState = ['MET','MET']
            if hasattr(tx, 'intermediateState'):
                intermediateState = tx.intermediateState
            else:
                intermediateState = None
            for el in elementsInStr(str(tx.constraint)):
                newEl = Element(el,finalState,intermediateState,
                        model=databaseParticles)
                datasetElements.append(newEl)
        combos = itertools.combinations(datasetElements, 2)
        for x,y in combos:
            if x == y and _complainAboutOverlappingConstraints:
                errmsg ="Constraints (%s) and (%s) appearing in dataset %s:%s overlap "\
                        "(may result in double counting)." % \
                        (x,y,self.getID(),self.globalInfo.id )
                logger.error( errmsg )
                raise SModelSError ( errmsg )

    def __ne__ ( self, other ):
        return not self.__eq__ ( other )

    def __str__ ( self ):
        if self.dataInfo.dataId:
            ret = "Dataset %s: %s" % (self.dataInfo.dataId, ", ".join ( map ( str, self.txnameList ) ) )
        else:
            ret = "Dataset: %s" % (", ".join ( map ( str, self.txnameList ) ) )
        return ret

    def __repr__(self):
        if self.dataInfo.dataId:
            return self.dataInfo.dataId
        else:
            return 'Dataset'

    def __eq__ ( self, other ):
        if type ( other ) != type ( self ):
            return False
        if self.dataInfo != other.dataInfo:
            return False
        if len(self.txnameList ) != len ( other.txnameList ):
            return False
        return True

    def getType(self):
        """
        Return the dataset type (EM/UL)
        """

        return self.dataInfo.dataType

    def getID(self):
        """
        Return the dataset ID
        """

        return self.dataInfo.dataId

    def getLumi(self):
        """
        Return the dataset luminosity. If not defined for the dataset, use
        the value defined in globalInfo.lumi.
        """

        if hasattr(self,'lumi'):
            return self.lumi
        else:
            return self.globalInfo.lumi

    def getTxName(self,txname):
        """
        get one specific txName object.
        """
        for tn in self.txnameList:
            if tn.txName == txname:
                return tn
        return None

    def getEfficiencyFor(self,txname,mass):
        """
        Convenience function. Get efficiency for mass
        assuming no lifetime rescaling. Same as self.getTxName(txname).getEfficiencyFor(m)
        """
        txname = self.getTxName(txname)
        if txname:
            return txname.getEfficiencyFor(mass)
        else:
            return None

    def getValuesFor(self,attribute):
        """
        Returns a list for the possible values appearing in the ExpResult
        for the required attribute (sqrts,id,constraint,...).
        If there is a single value, returns the value itself.

        :param attribute: name of a field in the database (string).
        :return: list of unique values for the attribute
        """

        return getValuesForObj(self,attribute)

    def likelihood(self, nsig, deltas_rel=0.2, marginalize=False, expected=False ):
        """
        Computes the likelihood to observe nobs events,
        given a predicted signal "nsig", assuming "deltas_rel"
        error on the signal efficiency.
        The values observedN, expectedBG, and bgError are part of dataInfo.

        :param nsig: predicted signal (float)
        :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param expected: Compute expected instead of observed likelihood
        :returns: likelihood to observe nobs events (float)
        """
        obs = self.dataInfo.observedN
        if expected:
            obs = self.dataInfo.expectedBG

        m = Data( obs, self.dataInfo.expectedBG, self.dataInfo.bgError**2,
                       deltas_rel=deltas_rel )
        computer = LikelihoodComputer(m)
        return computer.likelihood(nsig, marginalize=marginalize)

    def lmax( self, deltas_rel=0.2, marginalize=False, expected=False,
              allowNegativeSignals = False ):
        """
        Convenience function, computes the likelihood at nsig = observedN - expectedBG,
        assuming "deltas_rel" error on the signal efficiency.
        The values observedN, expectedBG, and bgError are part of dataInfo.

        :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :param expected: Compute expected instead of observed likelihood
        :param allowNegativeSignals: if False, then negative nsigs are replaced with 0.

        :returns: likelihood to observe nobs events (float)
        """
        obs = self.dataInfo.observedN
        if expected:
            obs = self.dataInfo.expectedBG

        m = Data( obs, self.dataInfo.expectedBG, self.dataInfo.bgError**2,
                       deltas_rel=deltas_rel )
        computer = LikelihoodComputer(m)
        return computer.lmax(marginalize=marginalize, nll=False,
                             allowNegativeSignals = allowNegativeSignals )

    def chi2(self, nsig, deltas_rel=0.2, marginalize=False):
        """
        Computes the chi2 for a given number of observed events "nobs",
        given number of signal events "nsig", and error on signal "deltas".
        nobs, expectedBG and bgError are part of dataInfo.
        :param nsig: predicted signal (float)
        :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
        :param marginalize: if true, marginalize nuisances. Else, profile them.
        :return: chi2 (float)
        """

        m = Data(self.dataInfo.observedN, self.dataInfo.expectedBG,
                    self.dataInfo.bgError**2,deltas_rel=deltas_rel)
        computer = LikelihoodComputer(m)
        ret = computer.chi2(nsig, marginalize=marginalize)

        return ret

    def folderName(self):
        """
        Name of the folder in text database.
        """
        return os.path.basename( self.path )

    def getAttributes(self, showPrivate=False):
        """
        Checks for all the fields/attributes it contains as well as the
        attributes of its objects if they belong to smodels.experiment.

        :param showPrivate: if True, also returns the protected fields (_field)
        :return: list of field names (strings)

        """

        attributes = getAttributesFrom(self)

        if not showPrivate:
            attributes = list(filter(lambda a: a[0] != '_', attributes))

        return attributes

    def getUpperLimitFor(self,element=None,expected = False, txnames = None
                         ,compute=False,alpha=0.05,deltas_rel=0.2):
        """
        Returns the upper limit for a given element (or mass) and txname. If
        the dataset hold an EM map result the upper limit is independent of
        the input txname or mass.
        For UL results if an Element object is given the corresponding upper limit
        will be rescaled according to the lifetimes of the element intermediate particles.
        On the other hand, if a mass is given, no rescaling will be applied.

        :param txname: TxName object or txname string (only for UL-type results)
        :param element: Element object or mass array with units (only for UL-type results)
        :param alpha: Can be used to change the C.L. value. The default value is 0.05
                      (= 95% C.L.) (only for  efficiency-map results)
        :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
        :param expected: Compute expected limit, i.e. Nobserved = NexpectedBG
                         (only for efficiency-map results)
        :param compute: If True, the upper limit will be computed
                        from expected and observed number of events.
                        If False, the value listed in the database will be used
                        instead.
        :return: upper limit (Unum object)
        """


        if self.getType() == 'efficiencyMap':
            upperLimit =  self.getSRUpperLimit(expected=expected,alpha=alpha,compute=compute,
                                               deltas_rel=deltas_rel)
            if (upperLimit/fb).normalize()._unit:
                logger.error("Upper limit defined with wrong units for %s and %s"
                              %(self.globalInfo.id,self.getID()))
                return False
            else:
                return upperLimit

        elif self.getType() == 'upperLimit':
            if not txnames or not element:
                logger.error("A TxName and mass array must be defined when \
                             computing ULs for upper-limit results.")
                return False
            elif isinstance(txnames,list):
                if len(txnames) != 1:
                    logger.error("txnames must be a TxName object, a string or a list with a single Txname object")
                    return False
                else:
                    txname = txnames[0]
            else:
                txname = txnames

            if not isinstance(txname, txnameObj.TxName) and \
            not isinstance(txname, str):
                logger.error("txname must be a TxName object or a string")
                return False

            if not isinstance(element, list) and not isinstance(element,Element):
                logger.error("Element must be an element object or a mass array")
                return False

            for tx in self.txnameList:
                if tx == txname or tx.txName == txname:
                    upperLimit = tx.getULFor(element,expected)

            return upperLimit

        else:
            logger.warning("Unkown data type: %s. Data will be ignored.",
                           self.getType())
            return None

    def getSRUpperLimit(self,alpha = 0.05, expected = False, compute = False, deltas_rel=0.2):
        """
        Computes the 95% upper limit on the signal*efficiency for a given dataset (signal region).
        Only to be used for efficiency map type results.

        :param alpha: Can be used to change the C.L. value. The default value is 0.05 (= 95% C.L.)
        :param expected: Compute expected limit ( i.e. Nobserved = NexpectedBG )
        :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
        :param compute: If True, the upper limit will be computed
                        from expected and observed number of events. If False, the value listed
                        in the database will be used instead.


        :return: upper limit value
        """

        if not self.getType() == 'efficiencyMap':
            logger.error("getSRUpperLimit can only be used for efficiency map results!")
            raise SModelSError()

        if not compute:
            if expected:
                try:
                    return self.dataInfo.expectedUpperLimit
                except AttributeError:
                    logger.info("expectedUpperLimit field not found. Returning None instead.")
                    return None
                    #return self.dataInfo.upperLimit
                    #logger.info("expectedUpperLimit field not found. Using observed UL instead.")
                    #return self.dataInfo.upperLimit
            else:
                return self.dataInfo.upperLimit

        Nobs = self.dataInfo.observedN  #Number of observed events
        if expected:
            Nobs = self.dataInfo.expectedBG
        Nexp = self.dataInfo.expectedBG  #Number of expected BG events
        bgError = self.dataInfo.bgError # error on BG

        m = Data(Nobs,Nexp,bgError**2,deltas_rel=deltas_rel)
        computer = UpperLimitComputer(cl=1.-alpha )
        maxSignalXsec = computer.ulSigma(m)
        if maxSignalXsec != None:
            maxSignalXsec = maxSignalXsec/self.getLumi()

        return maxSignalXsec

class CombinedDataSet(object):
    """
    Holds the information for a combined dataset (used for combining multiple datasets).
    """

    def __init__(self, expResult):

        self.path = expResult.path
        self.globalInfo = expResult.globalInfo
        self._datasets = expResult.datasets[:]
        self._marginalize = False
        self.sortDataSets()
        self.findType()

    def findType ( self ):
        """ find the type of the combined dataset """
        self.type = "bestSR" ## type of combined dataset, e.g. pyhf, or simplified
        if hasattr(self.globalInfo, "covariance"):
            self.type = "simplified"
        if hasattr(self.globalInfo, "jsonFiles" ):
            self.type = "pyhf"

    def __str__(self):
        ret = "Combined Dataset (%i datasets)" %len(self._datasets)
        return ret

    def sortDataSets(self):
        """
        Sort datasets according to globalInfo.datasetOrder.
        """
        if hasattr(self.globalInfo, "covariance"):
            datasets = self._datasets[:]
            if not hasattr(self.globalInfo, "datasetOrder" ):
                raise SModelSError("datasetOrder not given in globalInfo.txt for %s" % self.globalInfo.id )
            datasetOrder = self.globalInfo.datasetOrder
            if isinstance(datasetOrder,str):
                datasetOrder = [datasetOrder]

            if len(datasetOrder) != len(datasets):
                raise SModelSError("Number of datasets in the datasetOrder field does not match the number of datasets for %s"
                                   %self.globalInfo.id)
            for dataset in datasets:
                if not dataset.getID() in datasetOrder:
                    raise SModelSError("Dataset ID %s not found in datasetOrder" %dataset.getID())
                dsIndex = datasetOrder.index(dataset.getID())
                self._datasets[dsIndex] = dataset

    def getType(self):
        """
        Return the dataset type (combined)
        """

        return 'combined'

    def getID(self):
        """
        Return the ID for the combined dataset
        """

        return '(combined)'

    def getLumi(self):
        """
        Return the dataset luminosity. For CombinedDataSet always return
        the value defined in globalInfo.lumi.
        """

        return self.globalInfo.lumi

    def getDataSet(self,datasetID):
        """
        Returns the dataset with the corresponding dataset ID.
        If the dataset is not found, returns None.

        :param datasetID: dataset ID (string)

        :return: DataSet object if found, otherwise None.
        """

        for dataset in self._datasets:
            if datasetID == dataset.getID():
                return dataset

        return None
