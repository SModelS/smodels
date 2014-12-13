"""
.. module:: databaseFactory
   :synopsis: Holds the classes and methods to load the database and create the Infotxt and Smspy
              objects.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import logging, os, sys
from smodels.experiment import infoObjects, dataObjects, analysisObjects

FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.INFO)


class ExpResult(object):
    """
    Simple object to hols a pair of (Infotxt,Smspy) objects corresponding to an
    experimental result.
    
    :ivar path: path to the result folder
    :ivar info: Infotxt object
    :ivar smspy: Smspy object
    """
    def __init__(self, path):
        self.path = path
        self.info = infoObjects.InfoFile(os.path.join(path,"info.txt"))
        self.smspy = dataObjects.DataFile(os.path.join(path,"sms.py"),self.info)
    

class DataBase(object):    
    """
    Database object. Holds a collection of Infotxt and Smspy objects containing
    all the metainfo from the info.txt files and the corresponding data from the sms.py
    files.
    
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
        """This is the path to the base directory where to find the database.
        
        """
        return self._base
        
    def _validateBase(self, path):
        """Validates the base directory to locate the database. 
        Exits the script if something is wrong with the path.
    
        """
        logger.debug('Try to set the path for the database to: %s' %path)
        path = os.path.realpath(path) + '/'
        if not os.path.exists(path):
            logger.error('%s is no valid path!' %path)
            sys.exit()        
        return path
    
    @property
    def _getDatabaseVersion(self):
        """Retrieves the version of the database using the version file.
        """
        try:
            versionFile = open(self._base + '/version')
            content = versionFile.readlines()
            versionFile.close()
            logger.debug('Found version file %s with content %s' \
            %(self._base + '/version', content))
            return content[0].strip()
            
        except IOError:
            logger.error('There is no version file %s' \
            %self._base + '/version')
            return 'unknown version'
        
    @property
    def verbosity(self):
        """Tells the level the logger is set to.
        
        """
        return self._verbosity
        
    @verbosity.setter
    def verbosity(self, level):
        """Set the logger to specified level.
        
        """
        level = self._validateLevel(level)
        self._verbosity = level
        self._setLogLevel(level)
        
    def _validateLevel(self, level):
        """Validates given level for pythons logger module.
        
        """
        if not level.lower() in ['debug', 'info', 'warning', 'error']:
            logger.error('No valid level for verbosity: %s! Browser will ' +
                         'use default setting!' %level)
            return 'error'
        return level.lower()
            
    def _setLogLevel(self, level = 'error'):
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
        each (info.txt,sms.py) pair.
        
        :return: list of ExpResult objects    
        """
        
        resultsList = []
        for root, dirs, files in os.walk(self._base):
            if not 'info.txt' in files or not 'sms.py' in files:
                logger.debug("Missing files in %s" % root)
                continue
            else:
                resultsList.append(ExpResult(root))
                
        return resultsList
    
    def _getExpResults(self,analysisIDs=[],txnames=[]):
        """
        Returns a list of ExpResult objects containing the experimental results information.
        If analysisIDs is defined, returns only the results matching one of the IDs in the list.
        If txname is defined, returns only the results which contains
        at least one of the Tx names in the list.
        :param analysisID: list of analysis ids ([CMS-SUS-13-006,...])
        :param txname: list of txnames ([TChiWZ,...])
        :returns: list of ExpResult objects    
        """
        
        
        resultsList = []
        for expResult in self.expResultList:
            ID = expResult.info.getInfo('id')            
            if analysisIDs and not ID in analysisIDs: continue
            commonNames = [txname for txname in txnames if txname in expResult.info.getTxNames()]
            if not commonNames: continue
            resultsList.append(expResult)
        
        return resultsList
        
    def getAnalyses(self,analysisIDs=[],txnames=[]):
        """
        Returns a list of ULanalysis or EManalysis objects.
        Each object refers to a single analysisID+Txname (for UL analyses)
        or analysisID+SignalRegion (for a EM analyses).
        If analysisIDs is defined, returns only the results matching one of the IDs in the list.
        If txname is defined, returns only the results matching one of the Tx names in the list.
        :param analysisID: list of analysis ids ([CMS-SUS-13-006,...])
        :param txname: list of txnames ([TChiWZ,...])
        :returns: list of ULanalysis or EManalysis objects    
        """
        
        
        analysisList = []
        for expResult in self.expResultList:
            info = expResult.info
            smspy = expResult.smspy
            ID = info.getInfo('id')
            analysisType = info.getInfo('analysisType')            
            if analysisIDs and not ID in analysisIDs: continue
            if analysisType == 'EfficiencyMap':
                data = smspy.getData()
                newAna = analysisObjects.EManalysis(info,data)
                logger.info('Only upper limits analyses are accepted. Skipping %s' % ID)
                continue
                analysisList.append(newAna)
            elif analysisType == 'UpperLimit':                
                for txnameInfo in info.txNameInfoList:
                    if txnames and not txnameInfo.name in txnames: continue
                    data = smspy.getData(txnameInfo.name)
                    newAna = analysisObjects.ULanalysis(info.globalInfo,data,txnameInfo)            
                    analysisList.append(newAna)
        
        return analysisList
        