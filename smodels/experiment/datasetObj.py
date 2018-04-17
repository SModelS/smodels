"""
.. module:: datasetObj
   :synopsis: Holds the classes and methods used to read and store the information in the
              data folders.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


import os,glob,sys
from smodels.experiment import txnameObj,infoObj
from smodels.tools.physicsUnits import fb
from smodels.tools.SimplifiedLikelihoods import LikelihoodComputer, Model, UpperLimitComputer
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from smodels.tools.smodelsLogging import logger
from smodels.theory.particleNames import elementsInStr
from smodels.theory.element import Element
import itertools

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

    def checkForRedundancy ( self ):
        """ In case of efficiency maps, check if any txnames have overlapping
            constraints. This would result in double counting, so we dont 
            allow it. """
        if self.dataInfo.dataType == "upperLimit": 
            return False
        logger.debug ( "checking for redundancy" )
        datasetElements = []
        for tx in self.txnameList:
            for el in elementsInStr(tx.constraint):
                datasetElements.append(Element(el))
        combos = itertools.combinations ( datasetElements, 2 )
        for x,y in combos:
            if x.particlesMatch ( y ):
                errmsg ="Constraints (%s) appearing in dataset %s, %s overlap "\
                        "(may result in double counting)." % \
                        (x,self.dataInfo.dataId,self.globalInfo.id )
                logger.error( errmsg )
                raise SModelSError ( errmsg )
#                return True
#        return False

    def __ne__ ( self, other ):
        return not self.__eq__ ( other )

    def __str__ ( self ):
        ret = "Dataset: %s" % ( ", ".join ( map ( str, self.txnameList ) ) )
        return ret

    def __eq__ ( self, other ):
        if type ( other ) != type ( self ):
            return False
        if self.dataInfo != other.dataInfo:
            return False
        if len(self.txnameList ) != len ( other.txnameList ):
            return False
        return True

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
        convenience function.
        same as self.getTxName(txname).getEfficiencyFor(m)
        """
        txname = self.getTxName(txname)
        if txname: return txname.getEfficiencyFor(mass)
        return None

    def getValuesFor(self,attribute=None):
        """
        Returns a list for the possible values appearing in the DataSet
        for the required attribute.


        :param attribute: name of a field in the database (string). If not defined
                          it will return a dictionary with all fields and
                          their respective values
        :return: list of values
        """


        fieldDict = list ( self.__dict__.items() )
        valuesDict = {}
        while fieldDict:
            for field,value in fieldDict[:]:
                if not '<smodels.experiment' in str(value):
                    if not field in valuesDict: valuesDict[field] = [value]
                    else: valuesDict[field].append(value)
                else:
                    if isinstance(value,list):
                        for entry in value: fieldDict += list ( entry.__dict__.items() )
                    else: fieldDict += list ( value.__dict__.items() )
                fieldDict.remove((field,value))

        #Try to keep only the set of unique values
        for key,val in valuesDict.items():
            try:
                valuesDict[key] = list(set(val))
            except TypeError:
                pass
        if not attribute: return valuesDict
        elif not attribute in valuesDict:
            logger.warning("Could not find field %s in database" % attribute)
            return False
        else:
            return valuesDict[attribute]

    def likelihood ( self, nsig, deltas=None):
        """
        Computes the likelihood to observe nobs events,
        given a predicted signal "nsig", assuming "deltas"
        error on the signal efficiency.
        The values observedN, expectedBG, and bgError are part of dataInfo.
        :param nsig: predicted signal (float)
        :param deltas: uncertainty on signal (float).  If None, 
        default value (20%) will be used.
        :returns: likelihood to observe nobs events (float)
        """

        m = Model ( self.dataInfo.observedN, self.dataInfo.expectedBG, self.dataInfo.bgError**2 )
        computer = LikelihoodComputer ( m )
        return computer.likelihood( nsig, deltas )

    def folderName ( self ):
        """
        Name of the folder in text database.
        """
        return os.path.basename ( self.path )

#    this feature is not yet ready
#    def isUncorrelatedWith ( self, other ):
#        can it be safely assumed that this dataset is approximately
#        uncorrelated with "other"?
#        "other" can be a dataset or an expResult, in which case it is
#        true only if we are uncorrelated with all datasets of "other".
#
#        Two datasets of the same exp Result are considered never to be
#        uncorrelated.
#
#        if other == self: return False
#        if type(other) == type(self): ## comparing with another dataset
#            if self.globalInfo.path == other.globalInfo.path:
#                return False ## same expResult? -> correlated!
#            if self.globalInfo.dirName ( 1 ) != other.globalInfo.dirName ( 1 ):
#                ## different folders? uncorrelated!
#                return True
#            ## different expResults
#            return None ## FIXME implement
                

    def chi2( self, nsig, deltas_rel=None ):
        """
        Computes the chi2 for a given number of observed events "nobs",
        given number of signal events "nsig", and error on signal "deltas".
        nobs, expectedBG and bgError are part of dataInfo.
        :param nsig: predicted signal (float)
        :param deltas_rel: relative uncertainty in signal (float). 
        If None, default value (20%) will be used.
        :return: chi2 (float)
        """
        m = Model ( self.dataInfo.observedN, self.dataInfo.expectedBG, 
                    self.dataInfo.bgError**2, deltas_rel = deltas_rel )
        computer = LikelihoodComputer ( m )
        ret = computer.chi2( nsig )
        return ret

    def getAttributes(self,showPrivate=False):
        """
        Checks for all the fields/attributes it contains as well as the
        attributes of its objects if they belong to smodels.experiment.

        :param showPrivate: if True, also returns the protected fields (_field)
        :return: list of field names (strings)
        """

        fields = self.getValuesFor().keys()
        fields = list(set(fields))

        if not showPrivate:
            for field in fields[:]:
                if "_" == field[0]: fields.remove(field)

        return fields


    def getSRUpperLimit(self,alpha = 0.05, expected = False, compute = False ):
        """
        Computes the 95% upper limit on the signal*efficiency for a given dataset (signal region).
        Only to be used for efficiency map type results.

        :param alpha: Can be used to change the C.L. value. The default value is 0.05 (= 95% C.L.)
        :param expected: Compute expected limit ( i.e. Nobserved = NexpectedBG )
        :param compute: If True, the upper limit will be computed
                        from expected and observed number of events. If False, the value listed
                        in the database will be used instead.

        :return: upper limit value
        """

        if not self.dataInfo.dataType == 'efficiencyMap':
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
        lumi = self.globalInfo.lumi
        if (lumi*fb).normalize()._unit:
            ID = self.globalInfo.id
            logger.error("Luminosity defined with wrong units for %s" %(ID) )
            return False

        m = Model ( Nobs,Nexp,bgError )
        computer = UpperLimitComputer ( lumi, cl=1.-alpha )
        maxSignalXsec = computer.ulSigma( m )


        return maxSignalXsec


