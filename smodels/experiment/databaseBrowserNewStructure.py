"""
.. module:: databaseBrowserNewStructure
   :synopsis: Centralized facility to access smodels-database.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>

"""

import logging, os
#import setPath
import sys
from txObject import TxObject
from experimentIDObject import ExperimentIDObject
from analysisObject import AnalysisObject
from databaseBrowserException import DatabaseNotFoundException
from databaseBrowserException import InvalidRunRestrictionException
from databaseBrowserException import InvalidRunException
from databaseBrowserException import InvalidExperimentException
from databaseBrowserException import InvalidExperimentIDException
from databaseBrowserException import InvalidTxNameException
from databaseBrowserException import InvalidInfotxtFileException
from smodels.tools.physicsUnits import GeV
from infotxt import Infotxt


FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.ERROR)

class Browser(object):
    
    """Browses the database, exits if given path does not point to a valid 
    smodels-database. Browser can be restricted to specified run or experiment. 
    Verbosity can be set to specified level.
    
    """
    def __init__(self, base='/afs/hephy.at/user/w/walten/public/sms/'):
        self._allruns = ["8TeV", "7TeV"] # expand to 13TeV and 14TeV later
#        self._allexperiments = ["CMS", "ATLAS"]
#        self._artifacts = ['old', 'bad', 'missing', 'TODO', 'readme', 'SUCHI_RL_TEST']
        self._base = self._validateBase(base)
        self._experimentRestriction = None
        self._verbosity = 'error'
        self._databaseVersion = self._getDatabaseVersion
        self.database = self._getDatabase()
        self._runRestriction = None
        self._infos = {}
        self._experimentIDObjects = {}
        self._txObjects = {}
        self._analysisObjects = {}

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
            raise DatabaseNotFoundException("Invalid path")
        if not [run for run in os.listdir(path) if run in self._allruns]:
            logger.error('There is no valid database at %s' %path)
            raise DatabaseNotFoundException("Database not found")
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
        
    #@property
    #def experimentRestriction(self):
        #"""Tells if the browser is restricted to either CMS or ATLAS. 
        
        #"""
        #if self._experimentRestriction:
            #return self._experimentRestriction
        #return 'Browser will use CMS and ATLAS'
        
    #@experimentRestriction.setter
    #def experimentRestriction(self, detector):
        #"""Restricts the browser to either CMS or ATLAS.
        
        #"""
        #self._experimentRestriction = self._validateExperiment(detector)
        #for run in self.database:
            #self.database[run] = {detector: self.database[run][detector]}
        
    #@experimentRestriction.deleter
    #def experimentRestriction(self):
        #"""Removes the experimental restriction.
        
        #"""
        #self._experimentRestriction = None
        #self.database = self._getDatabase()
    
    #def _validateExperiment(self, detector):
        #"""Validates the given experiment. Raises an error if the given 
        #experiment is unknown.
        
        #"""
        #if not detector in ['CMS', 'ATLAS']:
            #logger.error('%s is no valid experiment!' %detector)
            #raise InvalidExperimentException("Invalid experiment")
        #logger.info('Focusing on experiment %s.' %detector)
        #return detector
        
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
        
    @property
    def _cache(self):
        """Retrieves all the dictionaries containing the experimental objects.
        Use the deleter to reset these.
        
        """
        expOb = {}
        expOb['experimentIDObjects'] = self._experimentIDObjects
        expOb['txObjects'] = self._txObjects
        expOb['analysisObjects'] = self._analysisObjects
        expOb['infos'] = self._infos

        return expOb
        
    @_cache.deleter
    def _cache(self):
        """Resets all the dictionaries containing the experimental objects in
        order to get all the logger messages from the building process again.
        
        """
        expOb['experimentIDObjects'] = {}
        expOb['txObjects'] = {}
        expOb['analysisObjects'] = {}
        expOb['infos'] = {}


    def _getDatabase(self):
        """Creates a nested dictionary containing all experiments, runs and all
        subdirectories resp. experimentIDs as keys and the info.txt objects and paths 
        to the sms.py as entries.
        :return: {run: {experiment: {experimentID: (infoObject, 'path/to/sms.py')}}}
    
        """
        data = {}
        for root, dirs, files in os.walk(self._base):
            if not 'info.txt' in files: continue
            if not 'sms.py' in files:
                pathToSms = None
            else:
                pathToSms = os.path.join(root, 'sms.py')
            pathToInfo = os.path.join(root, 'info.txt')
            logger.debug('Found info.txt in %s' %pathToInfo)
            try:
                info = Infotxt(pathToInfo)
            except InvalidInfotxtFileException: continue  
            run = info._run
            run = '%s*TeV' %float(run.split('*')[0])
            if not run in data:
                data[run] = {}
            experiment = info._experimentID.split('-')[0]
            if not experiment in data[run]:
                data[run][experiment] = {}    
            experimentID = info._experimentID
            data[run][experiment][experimentID] = [info, pathToSms]
        return data
        
    #@property
    #def runRestriction(self):
        #"""Tells if the browser is restricted to a specified run. 
        
        #"""
        #if not self._runRestriction:
            #return 'All runs are allowed.'
        #return self._runRestriction
        
    #@runRestriction.setter
    #def runRestriction(self, run):
        #"""Restricts the browser to one specified run. Don't use this lightly
        #it may cause serious changes in functionality.
        
        #"""
        #self._runRestriction = self._validateRun(run)
        #if not self._runRestriction:
            #logger.error('Failed to restrict browser to run: %s is not ' +
                         #'valid!' %run)
            #raise InvalidRunRestrictionException("Invalid run")
        #logger.warning('Browser restricted to run %s.' %run)
        #self.database = {key: self.database[key] for key in self.database if \
        #key == self._runRestriction}
    
    #@runRestriction.deleter
    #def runRestriction(self):
        #"""Deletes the restriction to one run.
        
        #"""
        #self._runRestriction = None
        #self.database = self._getDatabase()
        
    def _validateRun(self, run):
        """Validates the given run. Raises an error if the given run is unknown.
        
        """
        if not run: return None
        if not run in self.database:
            logger.error('%s is no valid run.' %run)
            raise InvalidRunException('Invalid run')

        return run
        
    def _validateExperimentID(self, experimentID):
        """Validates the given experimentID. Raises an error if the given experimentID 
        is unknown.
        
        """
        if not experimentID: return None
        for run in self.database:
            for exp in self.database[run]:
                for expID in self.database[run][exp]:
                    if expID == experimentID:
                        return experimentID
        logger.error('%s is no valid experimentID!' %experimentID)
        raise InvalidExperimentIDException
        
        
    def _validateTxName(self, txName):
        """Validates the given txName. Returns valid txNames and None if 
        the given txName is unknown.
        
        """
        
        if not txName in self.getTxNames():
            logger.error('%s is no valid txName!' %txName)
            raise InvalidTxNameException('Invalid txName')
            
        return txName
    
    def getRuns(self, experimentID = None, txName = None):
        """Retrieves all runs a given experimentID or txName is available for.
        Returns a list containing all runs or just a string when experimentID
        is given. 
        
        """
        
        if not experimentID and not txName:
            logger.info('No experimentID was given. Returnvalue will be list ' +
                        'containing all available runs!')
            return self.database.keys()
    
        #if self._runRestriction:
            #logger.warning('Cannot get all runs because browser is restricted ' +
                           #'to %s!' %self._runRestriction)
    
        #if self._experimentRestriction:
            #logger.warning('Browser is restricted to experiment %s!' \
            #%self._experimentRestriction)
        
        experimentID = self._validateExperimentID(experimentID)
        if experimentID and not txName:
            for run in self.database:
                for exp in self.database[run]:
                    if experimentID in self.database[run][exp]:
                        return run
        
        txName = self._validateTxName(txName)    
        if not experimentID and txName:
            runs = [key for key in self.database if \
            self.getTxNames(run = key) and txName in \
            self.getTxNames(run = key)]
            logger.info('No experimentID was given. There are %s runs ' %len(runs) +
                        'for given txName %s. ' %txName +
                        'Returnvalue will be list!' )
            return runs
        
        
    def getExperimentIDs(self, run = None, txName = None):
        """Retrieves all experimentIDs or all experimentIDs existing for given run or 
        run-txName-pair.
    
        """
    
        experimentIDs = []
    
        #if self._runRestriction:
            #logger.warning('Browser is restricted to run %s!' \
            #%self._runRestriction)
        
        if not run:
            runs = self.database.keys()
        
        if run and self._validateRun(run):
            runs = [run]
        
        for r in runs:
            for exp in self.database[r]:
                for expID in self.database[r][exp]:
                    if not txName:
                        experimentIDs.append(expID)
                    if txName and self._validateTxName(txName) and \
                        txName in self.database[r][exp][expID][0].txNames:
                            experimentIDs.append(expID)
                        
        return experimentIDs
    
    def getTxNames(self, run = None, experimentID = None):
        """Retrieves all txNames existing for given run or experimentID-run-pair.
    
        """
        
        topos = []
        
        if not run:
            runs = self.database.keys()
        
        if run and self._validateRun(run):
            runs = [run]
        
        for r in runs:
            for exp in self.database[r]:
                for expID in self.database[r][exp]:
                    if experimentID and self._validateExperimentID(experimentID):
                        if not expID == experimentID: continue
                    for t in self.database[r][exp][expID][0].txNames:
                        if not t in topos:
                            topos.append(t)
                        
        if not topos:
            logger.info('For runs %s and experimentIDs %s no txName could be found!' \
            %(runs, experimentIDs))
        
        return topos
    
    def getExperimentIDObject(self, experimentID):
        """This is the factory for the experimental experimentID object. 
        Returns None if it's not possible to build the experimental 
        experimentID object. 
        
        """
        if experimentID in self._experimentIDObjects:
            return self._experimentIDObjects[experimentID]
        
        run = self.getRuns(experimentID = experimentID)
        for exp in self.database[run]:
            if not experimentID in self.database[run][exp]:
                continue
        
            info = self.database[run][exp][experimentID][0]
            py = self.database[run][exp][experimentID][1]
        self._experimentIDObjects[experimentID] = ExperimentIDObject(experimentID, \
        info, run, None, py)
        return self._experimentIDObjects[experimentID]
        
    def getTxObject(self, txName):
        """This is the factory for the experimental txName object.
        
        """        
        if txName in self._txObjects:
            return self._txObjects[txName]
        txName = self._validateTxName(txName)
        txDict = self._txDict(txName)        
        self._txObjects[txName] = TxObject(txName, txDict)
        return self._txObjects[txName]
        
    def _txDict(self, txName):
        """Creates a nested dictionary that holds all the info.txt objects
        for each txName.
        :return: {'txName': {'experimentID': Infotxt(experimentID)}}
        
        """
        txDict = {}
        for r in self.getRuns(txName = txName):
            txDict[r] = {}
            for exp in self.database[r]:
                for expID in self.database[r][exp]:
                    if not expID in self.getExperimentIDs(run = r, txName = txName):
                        continue
                    txDict[r][expID] = self.database[r][exp][expID][0]
                
        return txDict
                
        
    def getAnalysisObject(self, experimentID, txName):
        """This is the factory for the experimental analysis object.
        It holds all plots for a pair of one experimentID and one txName
        :param experimentID: name of experimentID as string or experimentIDObject.name
        :param txName: name of txName as string or txObject.name
        
        """
        analysis = experimentID + '-' + txName
        logger.debug('Try to get analysis %s for %s-%s.' \
        %(analysis, experimentID, txName))
        if analysis in self._analysisObjects:
            logger.debug('Found analysis object for %s in dictionary.' \
            %analysis)
            return self._analysisObjects[analysis]
        experimentID = self._validateExperimentID(experimentID)
        expIDObj = self.getExperimentIDObject(experimentID)
        run = expIDObj.run
        exp = expIDObj.experiment
        py = self.database[run][exp][experimentID][1]
        if not txName in expIDObj.txNames:
            logger.error('There is no analysis object for ' +
                         'run-experimentID-txName: %s-%s-%s!' \
                         %(run, experimentID, txName))
            return None
        if not expIDObj.hasAxes:
            logger.error('There is no analysis object for ' +
                         'run-experimentID-txName: %s-%s-%s!' \
                         %(run, experimentID, txName))
            return None
        self._analysisObjects[analysis] = AnalysisObject(run, \
        expIDObj, self.getTxObject(txName), py)
        logger.debug('Built analysis object for %s-%s: %s' \
        %(experimentID, txName, self._analysisObjects[analysis]))
        return self._analysisObjects[analysis]
    
    def getPlotObject(self, experimentID, extendedTxName):
        """This is the factory for the experimental plot object.
        A plot referes to one slice in the 3D mass space for a pair of
        experimentID and txName.
        :param experimentID: name of experimentID as string or experimentIDObject.name
        :param extendedTxName: name of extended txName as string 
        (e.g. 'T6ttWWLSP050')
        
        """
        plot = experimentID + '-' + extendedTxName
        logger.debug('Try to get experimental plot %s for %s-%s.' \
        %(plot, experimentID, extendedTxName))
        expIDObj = self.getExperimentIDObject(experimentID)
        run = expIDObj.run
        if not expIDObj.hasAxes:
            logger.error('There is no experimental plot for ' +
            'run-experimentID-txName: %s-%s-%s!' %(run, experimentID, extendedTxName))
            return None
        for topo in expIDObj.extendedTxNames:
            if extendedTxName in expIDObj.extendedTxNames[topo]:
                return self.getAnalysisObject(experimentID, topo).plots[plot]
            else:
                continue
        logger.error('There is no experimental plot for ' +
                     'run-experimentID-txName: %s-%s-%s!' \
                     %(run, experimentID, extendedTxName))
        return None
    

            
