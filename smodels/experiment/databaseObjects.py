"""
.. module:: experiment.databaseObjects
   :synopsis: Contains classes and methods to load the database and create the InfoFile and DataFile
              objects as well as the list of analyses.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import logging
import os
from smodels.experiment import infoObject
from smodels.experiment import txnameObject
from smodels.experiment import datasetObject
from smodels.theory.auxiliaryFunctions import _memoize
from smodels.experiment.exceptions import DatabaseNotFoundException
from smodels.tools.physicsUnits import fb

logger = logging.getLogger(__name__)


class ExpResult(object):
    """
    Object  containing the information and data corresponding to an
    experimental result (experimental conference note or publication).
    
    :ivar path: path to the experimental result folder (i.e. ATLAS-CONF-2013-047)
    :ivar globalInfo: Info object holding the data in <path>/globalInfo.txt
    :ivar datasets: List of DataSet objects corresponding to the dataset folders in <path>    
    """
        
    def __init__(self, path=None):
        if path and os.path.isdir(path):
            self.path = path
            if not os.path.isfile(os.path.join(path, "globalInfo.txt")):
                logger.error("globalInfo.txt file not found in " + path)
                raise TypeError
            self.globalInfo = infoObject.Info(os.path.join(path, "globalInfo.txt"))
            self.datasets = []
            for root, _, files in os.walk(path):
                if 'dataInfo.txt' in files:  # data folder found
                    # Build data set
                    try:
                        dataset = datasetObject.DataSet(root, self.globalInfo)
                        self.datasets.append(dataset)
                    except TypeError:
                        continue

    def __str__(self):
        label = self.globalInfo.getInfo('id') + ": "
        dataIDs = [dataset.dataInfo.dataId for dataset in self.datasets]
        if dataIDs:
            for dataid in dataIDs:
                if dataid:
                    label += dataid + ","
        label = label[:-1]
        label += ':'
        txnames = []
        for dataset in self.datasets:
            for txname in dataset.txnameList:
                tx = txname.txName
                if not tx in txnames:
                    txnames.append(tx)
        if isinstance(txnames, list):
            for txname in txnames:
                label += txname + ','
        else:
            label += txnames + ','
        return label[:-1]


    def getTxNames(self):
        """
        Returns a list of all TxName objects appearing in all datasets.        
        """
        txnames = []
        for dataset in self.datasets:
            txnames += dataset.txnameList
        return txnames

    
    def getUpperLimitFor(self, dataID=None, alpha=0.05, expected=False,
                          txname=None, mass=None, compute=False):
        """
        Computes the 95% upper limit (UL) on the signal cross-section according to the type of result.
        For an Efficiency Map type, returns  the UL for the signal*efficiency for the given dataSet ID (signal region).
        For  an Upper Limit type,  returns the UL for the signal*BR for for the given mass array and Txname.
        
        :param dataID: dataset ID (string) (only for efficiency-map type results)
        :param alpha: Can be used to change the C.L. value. The default value is 0.05 (= 95% C.L.)
                      (only for  efficiency-map results)
        :param expected: Compute expected limit, i.e. Nobserved = NexpectedBG
                         (only for efficiency-map results)
        :param txname: TxName object (only for UL-type results)
        :param mass: Mass array with units (only for UL-type results)
        :param compute: If True, the upper limit will be computed
                        from expected and observed number of events. If False, the value listed
                        in the database will be used instead.
        
        
        :return: upper limit (Unum object)
        
        """
        if self.datasets[0].dataInfo.dataType == 'efficiencyMap':
            if not dataID or not isinstance(dataID, str):
                logger.error("The data set ID must be defined when computing ULs for\
                            efficiency-map results.")
                return False
            
            useDataset = False
            for dataset in self.datasets:
                if dataset.dataInfo.dataId == dataID:
                    useDataset = dataset
                    break
            if useDataset is False:
                logger.error("The data set ID not found.")
                return False
                
            if compute:
                upperLimit = useDataset.getSRUpperLimit(alpha, expected, compute)
            else:
                upperLimit = useDataset.dataInfo.upperLimit
                if (upperLimit/fb).normalize()._unit:
                    logger.error("Upper limit defined with wrong units for %s and %s"
                                  %(dataset.globalInfo.id,dataset.dataInfo.dataId))
                    return False           
            
        elif self.datasets[0].dataInfo.dataType == 'upperLimit':
            if not txname or not mass:
                logger.error("A TxName and mass array must be defined when computing ULs for\
                            upper-limit results.")
                return False
            if not isinstance(txname, txnameObject.TxName) and not isinstance(txname, str):
                logger.error("txname must be a TxName object or a string")
                return False
            if not isinstance(mass, list):
                logger.error("mass must be a mass array")
                return False
            for tx in self.getTxNames():
                if tx == txname or tx.txName == txname:
                    upperLimit = tx.txnameData.getValueFor(mass)
        else:
            logger.warning("Unkown data type: %s. Data will be ignored.", 
                           self.datasets[0].dataInfo.dataType)

        
        return upperLimit
    

    def getValuesFor(self, attribute=None):
        """
        Returns a list for the possible values appearing in the ExpResult
        for the required attribute (sqrts,id,constraint,...).
        If there is a single value, returns the value itself.
        
        :param attribute: name of a field in the database (string). If not
                          defined it will return a dictionary with all fields and 
                          their respective values
        :return: list of values or value
        
        """
        fieldDict = self.__dict__.items()[:]
        valuesDict = {}
        while fieldDict:
            for field, value in fieldDict[:]:
                if not '<smodels.experiment' in str(value):
                    if not field in valuesDict:
                        valuesDict[field] = [value]
                    else: valuesDict[field].append(value)
                else:
                    if isinstance(value, list):
                        for entry in value:
                            fieldDict += entry.__dict__.items()[:]
                    else: fieldDict += value.__dict__.items()[:]
                fieldDict.remove((field, value))

        # Try to keep only the set of unique values
        for key, val in valuesDict.items():
            try:
                valuesDict[key] = list(set(val))
            except:
                pass
        if not attribute:
            return valuesDict
        elif not attribute in valuesDict:
            logger.warning("Could not find field %s in %s", attribute, self.path)
            return False
        else:
            return valuesDict[attribute]


    def getAttributes(self, showPrivate=False):
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
                if "_" == field[0]:
                    fields.remove(field)
        return fields
    

    def getTxnameWith(self, restrDict={}):
        """
        Returns a list of TxName objects satisfying the restrictions.
        The restrictions specified as a dictionary.
        
        :param restrDict: dictionary containing the fields and their allowed values.
                          E.g. {'txname' : 'T1', 'axes' : ....}
                          The dictionary values can be single entries or a list of values.
                          For the fields not listed, all values are assumed to be allowed.
        :return: list of TxName objects if more than one txname matches the selection
        criteria or a single TxName object, if only one matches the selection.
        
        """
        txnameList = []
        for tag, value in restrDict.items():
            for txname in self.getTxNames():
                txval = txname.getInfo(tag)
                if txval is False:
                    continue
                elif txval == value:
                    txnameList.append(txname)

        if len(txnameList) == 1:
            txnameList = txnameList[0]

        return txnameList


class Database(object):
    """
    Database object. Holds a list of ExpResult objects.
    
    :ivar base: path to the database (string)
    :ivar expResultList: list of ExpResult objects 
        
    """
    
    
    def __init__(self, base=None):
        self._base = self._validateBase(base)
        self._verbosity = 'error'
        self._databaseVersion = self._getDatabaseVersion
        self.expResultList = self._loadExpResults()


    @property
    def databaseVersion(self):
        """
        The version of the database, read from the 'version' file.
        
        """
        return self._databaseVersion


    @property
    def base(self):
        """
        This is the path to the base directory where to find the database.
        
        """
        return self._base


    def _validateBase(self, path):
        """
        Validates the base directory to locate the database. 
        Exits the script if something is wrong with the path.
    
        """
        logger.debug('Try to set the path for the database to: %s', path)
        path = os.path.realpath(path) + '/'
        if not os.path.exists(path):
            logger.error('%s is no valid path!', path)
            raise DatabaseNotFoundException("Database not found")
        return path


    def __str__(self):
        idList = "Database: " + self.databaseVersion + "\n---------- \n"
        for expRes in self.expResultList:
            idList += expRes.globalInfo.getInfo('id') + '\n'
        return idList


    @property
    def _getDatabaseVersion(self):
        """
        Retrieves the version of the database using the version file.
        
        """
        try:
            versionFile = open(self._base + '/version')
            content = versionFile.readlines()
            versionFile.close()
            logger.debug('Found version file %s with content %s', self._base + '/version', content)
            return content[0].strip()

        except IOError:
            logger.error('There is no version file %s', self._base + '/version')
            return 'unknown version'


    @property
    def verbosity(self):
        """
        Tells the level the logger is set to.
        
        """
        return self._verbosity


    @verbosity.setter
    def verbosity(self, level):
        """
        Set the logger to specified level.
        
        """
        level = self._validateLevel(level)
        self._verbosity = level
        self._setLogLevel(level)


    def _validateLevel(self, level):
        """
        Validates given level for Python's logger module.
        
        """
        if not level.lower() in ['debug', 'info', 'warning', 'error']:
            logger.error('No valid level for verbosity: %s! Browser will ' +
                         'use default setting!' % level)
            return 'error'
        return level.lower()


    def _setLogLevel(self, level='error'):
        if level == 'debug':
            logger.setLevel(level=logging.DEBUG)
        if level == 'info':
            logger.setLevel(level=logging.INFO)
        if level == 'warning':
            logger.setLevel(level=logging.WARNING)
        if level == 'error':
            pass

    
    def _loadExpResults(self):
        """
        Checks the database folder and generates a list of ExpResult objects for
        each (globalInfo.txt,sms.py) pair.
        
        :return: list of ExpResult objects 
  
        """
        resultsList = []
        for root, _, files in os.walk(self._base):
            if not 'globalInfo.txt' in files:
                logger.debug("Missing files in %s", root)
                continue
            else:
                expres = ExpResult(root)
                if expres:
                    resultsList.append(expres)

        if not resultsList:
            logger.warning("Zero results loaded.")

        return resultsList


    def getExpResults(self, analysisIDs=['all'], datasetIDs=['all'], txnames=['all'], dataTypes = ['all']):
        """
        Returns a list of ExpResult objects.
        
        Each object refers to an analysisID containing one (for UL) or more (for Efficiency maps)
        dataset (signal region) and each dataset containing one or more TxNames.
        If analysisIDs is defined, returns only the results matching one of the IDs in the list.
        If dataTypes is defined, returns only the results matching a dataType in the list.
        If datasetIDs is defined, returns only the results matching one of the IDs in the list.
        If txname is defined, returns only the results matching one of the Tx names in the list.
        
        :param analysisID: list of analysis ids ([CMS-SUS-13-006,...])
        :param dataType: dataType of the analysis (all, efficiencyMap or upperLimit)
        :param datasetIDs: list of dataset ids ([ANA-CUT0,...])
        :param txnames: list of txnames ([TChiWZ,...])
        :returns: list of ExpResult objects or the ExpResult object if the list contains
                   only one result
                   
        """
        
        expResultList = []
        for expResult in self.expResultList:
            ID = expResult.globalInfo.getInfo('id')
            # Skip analysis not containing any of the required ids:
            if analysisIDs != ['all']:
                if not ID in analysisIDs:
                    continue              
            newExpResult = ExpResult()
            newExpResult.path = expResult.path
            newExpResult.globalInfo = expResult.globalInfo
            newExpResult.datasets = []            
            for dataset in expResult.datasets:
                if dataTypes != ['all']:
                    if not dataset.dataInfo.dataType in dataTypes:
                        continue
                if datasetIDs != ['all']:
                    if not dataset.dataInfo.dataId in datasetIDs:
                        continue
                newDataSet = datasetObject.DataSet(dataset.path, dataset.globalInfo)
                newDataSet.dataInfo = dataset.dataInfo
                newDataSet.txnameList = []
                for txname in dataset.txnameList:
                    if txnames != ['all']:
                        if not txname.txName in txnames:
                            continue
                    newDataSet.txnameList.append(txname)
                # Skip data set not containing any of the required txnames:
                if not newDataSet.txnameList:
                    continue
                newExpResult.datasets.append(newDataSet)
            # Skip analysis not containing any of the required txnames:
            if not newExpResult.getTxNames():
                continue
            expResultList.append(newExpResult)
        return expResultList
