"""
.. module:: datasetObj
   :synopsis: Holds the classes and methods used to read and store the information in the
              data folders.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


import os,glob
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

    def __init__(self, path=None, info=None, createInfo=True, discard_zeroes=True):
        """ :param discard_zeroes: discard txnames with zero-only results """

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
                    txname = txnameObj.TxName(txtfile,self.globalInfo,self.dataInfo)
                    if discard_zeroes and txname.hasOnlyZeroes():
                        logger.debug ( "%s, %s has only zeroes. discard it." % \
                                         ( self.path, txname.txName ) )
                        continue
                    self.txnameList.append(txname)
                except TypeError: continue

            self.txnameList.sort()
            self.checkForRedundancy()

    def checkForRedundancy(self):
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
            for el in elementsInStr(str(tx.constraint)):
                newEl = Element(el,finalState)
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
        given a predicted signal "nsig", assuming "deltas"
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
                    logger.info("expectedUpperLimit field not found. Using observed UL instead.")
                    return self.dataInfo.upperLimit
            else:
                return self.dataInfo.upperLimit

        Nobs = self.dataInfo.observedN  #Number of observed events
        if expected:
            Nobs = self.dataInfo.expectedBG
        Nexp = self.dataInfo.expectedBG  #Number of expected BG events
        bgError = self.dataInfo.bgError # error on BG        

        m = Data(Nobs,Nexp,bgError,detlas_rel=deltas_rel)
        computer = UpperLimitComputer(cl=1.-alpha )
        maxSignalXsec = computer.ulSigma(m)
        maxSignalXsec = maxSignalXsec/self.globalInfo.lumi

        return maxSignalXsec



class CombinedDataSet(object):
    """
    Holds the information for a combined dataset (used for combining multiple datasets).    
    """
    
    def __init__(self, expResult):
        
        self.globalInfo = expResult.globalInfo
        self._datasets = expResult.datasets[:]
        self._marginalize = False        
        self.sortDataSets()

    def __str__(self):
        ret = "Combined Dataset (%i datasets)" %len(self._datasets)
        return ret

                
                
    def sortDataSets(self):
        """
        Sort datasets according to globalInfo.datasetOrder.
        """
        
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
    
        
    def getCombinedUpperLimitFor(self, nsig, expected=False, deltas_rel=0.2):
        """
        Get combined upper limit.
        
        :param nsig: list of signal events in each signal region/dataset. The list
                        should obey the ordering in globalInfo.datasetOrder.
        :param expected: return expected, not observed value
        :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.        
        
        :returns: upper limit on sigma*eff
        """
        
        
        if not hasattr(self.globalInfo, "covariance" ):
            logger.error ( "no covariance matrix given in globalInfo.txt for %s" % self.globalInfo.id )
            raise SModelSError( "no covariance matrix given in globalInfo.txt for %s" % self.globalInfo.id )
        cov = self.globalInfo.covariance
        if type(cov) != list:
            raise SModelSError( "covariance field has wrong type: %s" % type(cov))
        if len(cov) < 1:
            raise SModelSError( "covariance matrix has length %d." % len(cov))

        computer = UpperLimitComputer(ntoys=10000)
        
        nobs = [x.dataInfo.observedN for x in self._datasets]
        bg = [x.dataInfo.expectedBG for x in self._datasets]
        no = nobs
        
        ret = computer.ulSigma(Data(observed=no, backgrounds=bg, covariance=cov, 
                                     third_moment=None, nsignal=nsig, deltas_rel=deltas_rel), 
                               marginalize=self._marginalize,
                               expected=expected)
        
        #Convert limit on total number of signal events to a limit on sigma*eff
        ret = ret/self.globalInfo.lumi
        
        return ret
        
    
    def combinedLikelihood(self, nsig, marginalize=False, deltas_rel=0.2):
        """
        Computes the (combined) likelihood to observe nobs events, given a
        predicted signal "nsig", with nsig being a vector with one entry per
        dataset.  nsig has to obey the datasetOrder. Deltas is the error on
        the signal.
        :param nsig: predicted signal (list)
        :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
        
        :returns: likelihood to observe nobs events (float)
        """
        
        if len(self._datasets) == 1:
            if isinstance(nsig,list):
                nsig = nsig[0]
            return self._datasets[0].likelihood(nsig,marginalize=marginalize)
        
        if not hasattr(self.globalInfo, "covariance" ):
            logger.error("Asked for combined likelihood, but no covariance error given." )
            return None

        nobs = [ x.dataInfo.observedN for x in self._datasets]
        bg = [ x.dataInfo.expectedBG for x in self._datasets]
        cov = self.globalInfo.covariance
            
        
        computer = LikelihoodComputer(Data(nobs, bg, cov, None, nsig, deltas_rel=deltas_rel))
        
        return computer.likelihood(nsig, marginalize=marginalize )

    def totalChi2(self, nsig, marginalize=False, deltas_rel=0.2):
        """
        Computes the total chi2 for a given number of observed events, given a
        predicted signal "nsig", with nsig being a vector with one entry per
        dataset. nsig has to obey the datasetOrder. Deltas is the error on
        the signal efficiency.
        :param nsig: predicted signal (list)
        :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.
        
        :returns: chi2 (float)
        """
        
        if len(self._datasets) == 1:
            if isinstance(nsig,list):
                nsig = nsig[0]            
            return self._datasets[0].chi2(nsig, marginalize=marginalize)
        
        if not hasattr(self.globalInfo, "covariance" ):
            logger.error("Asked for combined likelihood, but no covariance error given." )
            return None
        
        nobs = [x.dataInfo.observedN for x in self._datasets ]
        bg = [x.dataInfo.expectedBG for x in self._datasets ]
        cov = self.globalInfo.covariance
        
        computer = LikelihoodComputer(Data(nobs, bg, cov, deltas_rel=deltas_rel))
        
        return computer.chi2(nsig, marginalize=marginalize)
