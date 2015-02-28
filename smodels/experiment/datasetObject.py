"""
.. module:: datasetObjects
   :synopsis: Holds the classes and methods used to read and store the information in the
              data folders.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


import logging,os,glob
from smodels.experiment import txnameObject,infoObject
from smodels.tools import statistics
from smodels.theory.auxiliaryFunctions import _memoize

FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.DEBUG)


class DataSet(object):
    """Holds the information to a data folder (TxName objects, dataInfo,...)
    """
    
    def __init__(self, path,infoObj):
        self.dataDir = path
        self.info = infoObj
        self.txnameList = []
        
        logger.debug('Creating object based on data folder : %s' %self.dataDir)
        
        #Get data folder info:
        if not os.path.isfile(os.path.join(path,"dataInfo.txt")):
            logger.error("dataInfo.txt file not found in " + path)
            raise TypeError
        self.dataInfo = infoObject.Info(os.path.join(path,"dataInfo.txt"))

        #Get list of TxName objects:
        for txtfile in glob.iglob(os.path.join(path,"*.txt")):
            try:
                txname = txnameObject.TxName(txtfile,self.info)
                self.txnameList.append(txname)
            except TypeError: continue
            

    def getValuesFor(self,attribute=None):
        """
        Returns a list for the possible values appearing in the DataSet
        for the required attribute.
        If there is a single value, returns the value itself.
        
        :param attribute: name of a field in the database (string). If not defined
                          it will return a dictionary with all fields and their respective
                          values
        :return: list of values or single value
        """
        
        
        fieldDict = self.__dict__.items()[:]
        valuesDict = {}
        while fieldDict:
            for field,value in fieldDict[:]:
                if not '<smodels.experiment' in str(value):
                    if not field in valuesDict: valuesDict[field] = [value]
                    else: valuesDict[field].append(value)              
                else:
                    if isinstance(value,list):
                        for entry in value: fieldDict += entry.__dict__.items()[:]
                    else: fieldDict += value.__dict__.items()[:]
                fieldDict.remove((field,value))                

        #Try to keep only the set of unique values
        for key,val in valuesDict.items():
            try: valuesDict[key] = list(set(val))
            except: pass
        if not attribute: return valuesDict
        elif not attribute in valuesDict:
            logger.warning("Could not find field %s in database" % attribute)
            return False
        else:
            if len(valuesDict[attribute]) == 1: return valuesDict[attribute][0]
            else:
                return valuesDict[attribute]
            

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
    
    @_memoize
    def getUpperLimit(self,alpha = 0.05, expected = False ):
        """
        Computes the 95% upper limit on the signal*efficiency for an efficiency-map
        type data set
        Only to be used for efficiency map type results.
        :param alpha: Can be used to change the C.L. value. The default value is 0.05 (= 95% C.L.)
        :param expected: Compute expected limit ( i.e. Nobserved = NexpectedBG )
        
        :return: upper limit value 
        """
        
        Nobs = self.dataInfo.observedN  #Number of observed events
        if expected: 
            Nobs = self.dataInfo.expectedBG 
        Nexp = self.dataInfo.expectedBG  #Number of expected BG events
        bgError = self.dataInfo.bgError # error on BG
        lumi = self.info.lumi
        maxSignalXsec = statistics.upperLimit (Nobs,Nexp,bgError,lumi,alpha)
                
        return maxSignalXsec
            
    
