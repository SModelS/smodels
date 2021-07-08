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
        self.bestCB = None # To store the index of the best combination
        self.findType()

    def findType ( self ):
        """ find the type of the combined dataset """
        self.type = None ## type of combined dataset, e.g. pyhf, or simplified
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

    def lmax ( self, nsig, marginalize, deltas_rel, nll=False, expected=False,
               allowNegativeSignals = False ):
        """ compute likelihood at maximum """
        if self.type == "simplified":
            nobs = [ x.dataInfo.observedN for x in self._datasets]
            if expected:
                # nobs = [ x.dataInfo.expectedBG for x in self._datasets]
                nobs = [ int(np.round(x.dataInfo.expectedBG)) for x in self._datasets]
            bg = [ x.dataInfo.expectedBG for x in self._datasets]
            cov = self.globalInfo.covariance
            if type(nsig) in [ list, tuple ]:
                nsig = np.array(nsig)
            computer = LikelihoodComputer(Data(nobs, bg, cov, None, nsig, deltas_rel=deltas_rel))
            mu_hat = computer.findMuHat ( nsig, allowNegativeSignals=allowNegativeSignals )
            musig = nsig * mu_hat
            return computer.likelihood ( musig, marginalize=marginalize, nll=nll )
        if self.type == "pyhf":
            ulcomputer = self.getPyhfComputer( nsig )
            combinations = ulcomputer.data.combinations
            return ulcomputer.lmax ( nll=False )
        return -1.

    def getPyhfComputer ( self, nsig ):
        """ create the pyhf ul computer object
        :returns: pyhf upper limit computer, and combinations of signal regions
        """
        # Getting the path to the json files
        jsonFiles = [js for js in self.globalInfo.jsonFiles]
        jsons = self.globalInfo.jsons.copy()
        datasets = [ds.getID() for ds in self._datasets]
        total = sum(nsig)
        if total == 0.: # all signals zero? can divide by anything!
            total = 1.
        nsig = [s/total for s in nsig] # Normalising signals to get an upper limit on the events count
        # Filtering the json files by looking at the available datasets
        for jsName in self.globalInfo.jsonFiles:
            if all([ds not in self.globalInfo.jsonFiles[jsName] for ds in datasets]):
                # No datasets found for this json combination
                jsIndex = jsonFiles.index(jsName)
                jsonFiles.pop(jsIndex)
                jsons.pop(jsIndex)
                continue
            if not all([ds in datasets for ds in self.globalInfo.jsonFiles[jsName]]):
                # Some SRs are missing for this json combination
                logger.error("Wrong json definition in globalInfo.jsonFiles for json : %s" % jsName)
        logger.debug("list of datasets: {}".format(datasets))
        logger.debug("jsonFiles after filtering: {}".format(jsonFiles))
        # Constructing the list of signals with subsignals matching each json
        nsignals = list()
        for jsName in jsonFiles:
            subSig = list()
            for srName in self.globalInfo.jsonFiles[jsName]:
                try:
                    index = datasets.index(srName)
                except ValueError:
                    logger.error("%s signal region provided in globalInfo is not in the list of datasets" % srName)
                sig = nsig[index]
                subSig.append(sig)
            nsignals.append(subSig)
        # Loading the jsonFiles, unless we already have them (because we pickled)
        from smodels.tools.pyhfInterface import PyhfData, PyhfUpperLimitComputer
        data = PyhfData(nsignals, jsons, jsonFiles )
        if data.errorFlag: return None
        ulcomputer = PyhfUpperLimitComputer(data)
        return ulcomputer

    def getCombinedUpperLimitFor(self, nsig, expected=False, deltas_rel=0.2):
        """
        Get combined upper limit. If covariances are given in globalInfo then simplified likelihood is used, else if json files are given pyhf cimbination is performed.

        :param nsig: list of signal events in each signal region/dataset. The list
                        should obey the ordering in globalInfo.datasetOrder.
        :param expected: return expected, not observed value
        :param deltas_rel: relative uncertainty in signal (float). Default value is 20%.

        :returns: upper limit on sigma*eff
        """

        if self.type == "simplified":
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

            if ret != None:
                #Convert limit on total number of signal events to a limit on sigma*eff
                ret = ret/self.getLumi()
            logger.debug("SL upper limit : {}".format(ret))
            return ret
        elif self.type == "pyhf":
            logger.debug("Using pyhf")
            if all([s == 0 for s in nsig]):
                logger.warning("All signals are empty")
                return None
            ulcomputer = self.getPyhfComputer( nsig )
            combinations = ulcomputer.data.combinations
            if ulcomputer.nWS == 1:
                ret = ulcomputer.ulSigma(expected=expected)
                ret = ret/self.getLumi()
                logger.debug("pyhf upper limit : {}".format(ret))
                return ret
            else:
                # Looking for the best combination
                logger.debug('self.bestCB : {}'.format(self.bestCB))
                if self.bestCB == None:
                    logger.debug("Performing best expected combination")
                    ulMin = float('+inf')
                    for i_ws in range(ulcomputer.nWS):
                        ul = ulcomputer.ulSigma(expected=True, workspace_index=i_ws)
                        if ul == None:
                            continue
                        if ul < ulMin:
                            ulMin = ul
                            i_best = i_ws
                    self.bestCB = combinations[i_best] # Keeping the index of the best combination for later
                    logger.debug('Best combination : %s' % self.bestCB)
                # Computing upper limit using best combination
                if expected:
                    try:
                        ret = ulMin/self.getLumi()
                    except NameError:
                        ret = ulcomputer.ulSigma(expected=True, workspace_index=combinations.index(self.bestCB))
                        ret = ret/self.getLumi()
                else:
                    ret = ulcomputer.ulSigma(expected=False, workspace_index=combinations.index(self.bestCB))
                    ret = ret/self.getLumi()
                logger.debug("pyhf upper limit : {}".format(ret))
                return ret
        else:
            logger.error ( "no covariance matrix or json file given in globalInfo.txt for %s" % self.globalInfo.id )
            raise SModelSError( "no covariance matrix or json file given in globalInfo.txt for %s" % self.globalInfo.id )

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



        if hasattr(self.globalInfo, "covariance" ):
            if len(self._datasets) == 1:
                if isinstance(nsig,list):
                    nsig = nsig[0]
                return self._datasets[0].likelihood(nsig,marginalize=marginalize)
            nobs = [ x.dataInfo.observedN for x in self._datasets]
            bg = [ x.dataInfo.expectedBG for x in self._datasets]
            cov = self.globalInfo.covariance
            computer = LikelihoodComputer(Data(nobs, bg, cov, None, nsig, deltas_rel=deltas_rel))
            return computer.likelihood(nsig, marginalize=marginalize )
        elif hasattr(self.globalInfo, "jsonFiles"):
            # Getting the path to the json files
            # Loading the jsonFiles
            ulcomputer = self.getPyhfComputer( nsig )
            combinations = ulcomputer.data.combinations
            if ulcomputer.nWS == 1:
                return ulcomputer.likelihood()
            else:
                # Looking for the best combination
                if self.bestCB == None:
                    ulMin = float('+inf')
                    for i_ws in range(ulcomputer.nWS):
                        logger.debug("Performing best expected combination")
                        ul = ulcomputer.ulSigma(expected=True, workspace_index=i_ws)
                        if  ul < ulMin:
                            ulMin = ul
                            i_best = i_ws
                    self.bestCB = combinations[i_best] # Keeping the index of the best combination for later
                    logger.debug('Best combination : %d' % self.bestCB)
                return ulcomputer.likelihood(workspace_index=combinations.index(self.bestCB))
        else:
            logger.error("Asked for combined likelihood, but no covariance or json file given." )
            return None

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

        if hasattr(self.globalInfo, "covariance" ):
            if len(self._datasets) == 1:
                if isinstance(nsig,list):
                    nsig = nsig[0]
                return self._datasets[0].chi2(nsig, marginalize=marginalize)
            nobs = [x.dataInfo.observedN for x in self._datasets ]
            bg = [x.dataInfo.expectedBG for x in self._datasets ]
            cov = self.globalInfo.covariance

            computer = LikelihoodComputer(Data(nobs, bg, cov, deltas_rel=deltas_rel))

            return computer.chi2(nsig, marginalize=marginalize)
        elif hasattr(self.globalInfo, "jsonFiles"):
            # Getting the path to the json files
            # Loading the jsonFiles
            ulcomputer = self.getPyhfComputer( nsig )
            combinations = ulcomputer.data.combinations
            if ulcomputer.nWS == 1:
                return ulcomputer.chi2()
            else:
                # Looking for the best combination
                if self.bestCB == None:
                    ulMin = float('+inf')
                    for i_ws in range(ulcomputer.nWS):
                        logger.debug("Performing best expected combination")
                        ul = ulcomputer.ulSigma(expected=True, workspace_index=i_ws)
                        if  ul < ulMin:
                            ulMin = ul
                            i_best = i_ws
                    self.bestCB = combinations[i_best] # Keeping the index of the best combination for later
                    logger.debug('Best combination : %d' % self.bestCB)
                return ulcomputer.chi2(workspace_index=combinations.index(self.bestCB))
        else:
            logger.error("Asked for combined likelihood, but no covariance error given." )
            return None
