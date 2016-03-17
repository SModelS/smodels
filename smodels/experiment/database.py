#!/usr/bin/env python

"""
.. module:: experiment.database
   :synopsis: Contains Database class that represents the database of
              experimental results.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
import cPickle as pickle
import logging
import os
import time
from smodels.experiment import infoObject
from smodels.experiment import txnameObject
from smodels.experiment import datasetObject
from smodels.experiment.expResult import ExpResult
from smodels.theory.auxiliaryFunctions import _memoize
from smodels.experiment.exceptions import DatabaseNotFoundException
from smodels.tools.physicsUnits import fb

logger = logging.getLogger(__name__)

class Database(object):
    """
    Database object. Holds a list of ExpResult objects.
    
    :ivar base: path to the database (string)
    :ivar expResultList: list of ExpResult objects 
        
    """
    
    def __init__(self, base=None, force_load = None, verbosity='info' ):
        """
        :param force_load: force loading the text database ("txt"),
            or pickle database ("pcl"), dont force anything if None
        """
        self._base = self._validateBase(base)
        self._verbosity = verbosity
        self._databaseVersion = None
        self.expResultList = None
        self.txt_mtime = None, None
        self.pcl_mtime = None, None
        self.pcl_db = None
        self.pclfile = os.path.join ( self._base, "database.pcl" )
        self._setLogLevel ( self._verbosity )
        if force_load=="txt":
            self.loadTextDatabase()
        if force_load=="pcl":
            self.loadPickleFile()
        if force_load==None:
            self.loadDatabase()

    def __eq__ ( self, other ):
        """ compare two database 
        """
        if self.databaseVersion != other.databaseVersion:
            return False
        if len(self.expResultList ) != len (other.expResultList):
            return False
        for ( myres, otherres ) in zip ( self.expResultList, other.expResultList ):
            if myres != otherres:
                return False
        return True

    def loadDatabase ( self ):
        """ if no pickle file is available, then 
            load the database and create the pickle file.
            if pickle file is available, then check if
            it needs update, create new pickle file, in
            case it does need an update.
        """
        if not os.path.exists ( self.pclfile ):
            self.loadTextDatabase()
            self.createPickleFile()
        else:
            if self.needsUpdate():
                self.createPickleFile()
            else:
                self.loadPickleFile( lastm_only = False )

    def loadTextDatabase ( self ):
        """ simply loads the textdabase """
        self._databaseVersion = self._getDatabaseVersion
        self.expResultList = self._loadExpResults()

    def lastModifiedDir ( self, dirname, lastm ):
        """ return the last modified timestamp of dirname,
            working recursively 
        :param dirname: directory name that is checked
        :param lastm: the most recent timestamp so far
        :returns: the most recent timestamp, and the number of files
        """
        ret = lastm
        ctr=0
        for f in os.listdir ( dirname ):
            lf = os.path.join ( dirname, f )
            if os.path.isdir ( lf ):
                (ret,tctr) = self.lastModifiedDir ( lf, ret )
                ctr+=tctr+1
            else:
                ctr+=1
                tmp = os.stat ( lf ).st_mtime
                if tmp > ret:
                    ret = tmp
        return (ret,ctr)

    def lastModifiedAndFileCount( self ):
        if self.txt_mtime[0]:
            ## already evaluated
            return
        versionfile = os.path.join ( self._base, "version" ) 
        if not os.path.exists ( versionfile ):
            logger.error("%s does not exist." % versionfile )
            sys.exit()
        lastm = os.stat(versionfile).st_mtime
        count=1
        topdir = os.listdir ( self._base )
        for File in topdir:
            subdir = os.path.join ( self._base, File )
            if not os.path.isdir ( subdir ) or File in [ ".git" ]:
                continue
            (lastm,tcount) = self.lastModifiedDir ( subdir, lastm )
            count+=tcount+1
        self.txt_mtime = lastm, count

    def loadPickleFile ( self, lastm_only = False ):
        """ load a pickle file, returning
            last modified, file count, database.
        :param lastm_only: if true, the database itself is not read.
        :returns: database object, or None, lastm_only == True.
        """
        if lastm_only and self.pcl_mtime[0]:
            ## doesnt need to load database, and mtime is already
            ## loaded
            return None

        if self.pcl_db:
            return self.pcl_db

        with open ( self.pclfile, "r" ) as f:
            t0=time.time()
            self.pcl_mtime = pickle.load ( f )
            self._databaseVersion = pickle.load ( f )
            # self.pclfile = pickle.load ( f )
            if not lastm_only:
                logger.info ( "loading pickle file %s" % self.pclfile )
                self.expResultList = pickle.load ( f )
                t1=time.time()-t0
                logger.info ( "Loaded database from %s in %.1f secs." % \
                        ( self.pclfile, t1 ) )
        return self

    def checkPickleFile ( self ):
        nu=self.needsUpdate()
        logger.debug ( "Checking pickle file." )
        logger.debug ( "Pickle file dates to %s(%d)" % \
                      ( time.ctime(self.pcl_mtime[0]),self.pcl_mtime[1] ) )
        logger.debug ( "Database dates to %s(%d)" % \
                      ( time.ctime(self.txt_mtime[0]),self.txt_mtime[1] ) )
        if nu:
            logger.info ( "pickle file needs an update." )
        else:
            logger.info ( "pickle file does not need an update." )
        return nu

    def needsUpdate ( self ):
        """ does the pickle file need an update? """
        self.lastModifiedAndFileCount()
        self.loadPickleFile ( lastm_only = True )
        return ( self.txt_mtime[0] > self.pcl_mtime[0] or \
                 self.txt_mtime[1] != self.pcl_mtime[1] )


    def createPickleFile ( self, filename=None ):
        """ create a pcl file from the text database,
            potentially overwriting an old pcl file. """
        t0=time.time()
        logger.info ( "Creating pickle file (this may take a few minutes)" )
        logger.debug ( " * compute last modified timestamp." )
        self.lastModifiedAndFileCount()
        logger.debug (  " * compute timestamp: %s filecount: %d" % \
                ( time.ctime ( self.txt_mtime[0] ), self.txt_mtime[1] ) )
        pclfile = filename
        if pclfile == None:
            pclfile = self.pclfile
        logger.debug (  " * create %s" % self.pclfile )
        with open ( pclfile, "w" ) as f:
            pickle.dump ( self.txt_mtime, f )
            pickle.dump ( self._databaseVersion, f )
            # pickle.dump ( self.pclfile, f )
            logger.debug (  " * load text database" )
            self.loadTextDatabase() 
            logger.debug (  " * write %s" % self.pclfile )
            pickle.dump ( self.expResultList, f, protocol=2 )
            logger.info (  " * done writing %s in %.1f secs." % \
                    ( self.pclfile, time.time()-t0 ) )

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
        idList = "Database version: " + self.databaseVersion + "\n---------- \n"
        if self.expResultList == None:
            idList += "no experimental results available!"
        else:
            for expRes in self.expResultList:
                idList += expRes.globalInfo.getInfo('id') + '\n'
        return idList


    @property
    def _getDatabaseVersion(self):
        """
        Retrieves the version of the database using the version file.
        
        """
        try:
            vfile = os.path.join ( self._base, "version" )
            versionFile = open( vfile )
            content = versionFile.readlines()
            versionFile.close()
            line = content[0].strip()
            logger.debug("Found version file %s with content ``%s''" \
                   % ( vfile, line) )
            return line

        except IOError:
            logger.error('There is no version file %s', vfile )
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
        level = level.lower()
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
        
        :returns: list of ExpResult objects 
  
        """
        resultsList = []
        for root, _, files in os.walk(self._base):
            if "/.git/" in root:
                continue
            if root[-11:] == "/validation":
                continue
            if root[-5:] == "/orig":
                continue
            if not 'globalInfo.txt' in files:
                logger.debug("Missing globalInfo.txt in %s", root)
                continue
            else:
                expres = ExpResult(root)
                if expres:
                    resultsList.append(expres)

        if not resultsList:
            logger.warning("Zero results loaded.")

        return resultsList


    def getExpResults(self, analysisIDs=['all'], datasetIDs=['all'], txnames=['all'],
                      dataTypes = ['all']):
        """
        Returns a list of ExpResult objects.
        
        Each object refers to an analysisID containing one (for UL) or more
        (for Efficiency maps) dataset (signal region) and each dataset
        containing one or more TxNames.  If analysisIDs is defined, returns
        only the results matching one of the IDs in the list.  If dataTypes is
        defined, returns only the results matching a dataType in the list.  If
        datasetIDs is defined, returns only the results matching one of the IDs
        in the list.  If txname is defined, returns only the results matching
        one of the Tx names in the list.
        
        :param analysisID: list of analysis ids ([CMS-SUS-13-006,...])
        :param dataType: dataType of the analysis (all, efficiencyMap or upperLimit)
        :param datasetIDs: list of dataset ids ([ANA-CUT0,...])
        :param txnames: list of txnames ([TChiWZ,...])
        :returns: list of ExpResult objects or the ExpResult object if the list
                  contains only one result
                   
        """
        # print ".getExpResultList()"
        expResultList = []
        for expResult in self.expResultList:
          #  print "expResult=",expResult.__str__()[:20]
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
                # print "creating newDataSet %s %s " % ( dataset.path, dataset.globalInfo )
                # print "dataset=",type(dataset)
                newDataSet = datasetObject.DataSet(dataset.path, dataset.globalInfo,False)
                # print "Done creating new Dataset"
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
        # print "done .getExpResultList()"
        return expResultList

    def updatePickleFile ( self ):
        """ write a pickle file, but only if 
            necessary. """
        if self.needsUpdate():
            logger.debug ( "pickle file needs an update." )
            self.createPickleFile()
        else:
            logger.debug ( "pickle file does not need an update." )
if __name__ == "__main__":
    import argparse
    """ Run as a script, this checks and/or writes database.pcl files """
    argparser = argparse.ArgumentParser(description='simple script to check \
            and/or write database.pcl files')
    argparser.add_argument('-c', '--check', help='check pickle file',
                           action='store_true')
    argparser.add_argument('-t', '--time', help='time reading db',
                           action='store_true')
    argparser.add_argument('-r', '--read', help='read pickle file',
                           action='store_true')
    argparser.add_argument('-w', '--write', help='force writing pickle file',
                           action='store_true')
    argparser.add_argument('-u', '--update', help='update pickle file, if necessary',
                           action='store_true')
    argparser.add_argument('-d', '--debug', help='debug mode',
                           action='store_true')
    argparser.add_argument('-D', '--database', help='directory name of database', 
                            default="../../../smodels-database/" )
    args = argparser.parse_args()
    logger.setLevel(level=logging.INFO )
    if args.debug:
        logger.setLevel(level=logging.DEBUG )
    if args.write:
        db = Database ( args.database, force_load="txt" )
        if args.debug:
            db.verbosity = "debug"
        logger.debug ( "%s" % db )
        db.createPickleFile()
        sys.exit()
    db = Database ( args.database )
    if args.debug:
        db.verbosity = "debug"
    logger.debug ( "%s" % db )
    if args.update:
        db.updatePickleFile()
    if args.check:
        db.checkPickleFile()
    if args.time:
        t0=time.time()
        expResult = db.loadPickleFile ( lastm_only = False )
        t1=time.time()
        print "Time it took reading pickle file: %.1f s." % (t1-t0)
        txtdb = db.loadTextDatabase()
        t2=time.time()
        print "Time it took reading text   file: %.1f s." % (t2-t1)
    if args.read:
        db = db.loadPickleFile ( lastm_only = False )
        listOfExpRes = db.getExpResults() 
        for expResult in listOfExpRes:
            print expResult
        #txtdb=pickler.loadTextDatabase()
        #listOfExpRes = txtdb.getExpResults() 
        #for expResult in listOfExpRes:
        #    print expResult

