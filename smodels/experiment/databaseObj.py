#!/usr/bin/env python

"""
.. module:: databaseObj
   :synopsis: Contains Database class that represents the database of experimental results.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
import os
import time
from smodels.experiment import datasetObj
from smodels.experiment.expResultObj import ExpResult
from smodels.experiment.exceptions import DatabaseNotFoundException
from smodels.tools.physicsUnits import fb, TeV

try:
    import cPickle as serializer
except ImportError as e:
    import pickle as serializer

from smodels.tools.smodelsLogging import logger, setLogLevel

class Database(object):
    """
    Database object. Holds a list of ExpResult objects.
    
    :ivar base: path to the database (string)
    :ivar force_load: force loading the text database ("txt"),
        or binary database ("pcl"), dont force anything if None
    :ivar expResultList: list of ExpResult objects 
        
    """
    
    def __init__(self, base=None, force_load = None, verbosity=None,
                 pclfilename = "database.pcl" ):
        """
        :param force_load: force loading the text database ("txt"),
            or binary database ("pcl"), dont force anything if None
        :param pclfilename: filename of the binary (pickled) database file
        """
        self.force_load = force_load
        self.pclfilename = pclfilename
        self.hasFastLim = False # True if any ExpResult is from fastlim
        self._validateBase(base)
        self._verbosity = verbosity 
        self._databaseVersion = None
        self.expResultList = []
        self.txt_mtime = None, None
        self.pcl_mtime = None, None
        self.pcl_db = None
        self.sw_format_version = "115" ## what format does the software support?
        self.pcl_format_version = None ## what format is in the binary file?
        self.binfile = os.path.join ( self._base, self.pclfilename )
        setLogLevel ( self._verbosity )
        if self.force_load=="txt":
            self.loadTextDatabase()
            self.printFastlimBanner()
            return
        if self.force_load=="pcl":
            self.loadBinaryFile()
            self.printFastlimBanner()
            return
        if self.force_load in [ None, "none", "None" ]:
            self.loadDatabase()
            self.printFastlimBanner()
            return
        logger.error ( "when initialising database: force_load=%s is not " \
                       "recognized. Valid values are: pcl, txt, None." % force_load )
        sys.exit()

    def printFastlimBanner ( self ):
        """ check if fastlim appears in data.
            If yes, print a statement to stdout. """
        if not self.hasFastLim: return
        logger.info ( "FastLim v1.1 efficiencies loaded. Please cite: arXiv:1402.0492, EPJC74 (2014) 11" )

    def __eq__ ( self, other ):
        """ compare two databases
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
        """ if no binary file is available, then 
            load the database and create the binary file.
            if binary file is available, then check if
            it needs update, create new binary file, in
            case it does need an update.
        """
        if not os.path.exists ( self.binfile ):
#             self.loadTextDatabase()
            self.createBinaryFile()
        else:
            if self.needsUpdate():
                self.createBinaryFile()
            else:
                self.loadBinaryFile( lastm_only = False )

    def loadTextDatabase ( self ):
        """ simply loads the textdabase """
        if self._databaseVersion and len(self.expResultList)>0:
            logger.debug ( "Asked to load database, but has already been loaded. Ignore." )
            return
        logger.info ( "Parsing text database at %s" % self._base )
        self._databaseVersion = self._getDatabaseVersion
        self.expResultList = self._loadExpResults()

    def lastModifiedDir ( self, dirname, lastm ):
        """
        Return the last modified timestamp of dirname, working recursively
         
        :param dirname: directory name that is checked
        :param lastm: the most recent timestamp so far
        :returns: the most recent timestamp, and the number of files
        """
        ret = lastm
        ctr=0
        for f in os.listdir ( dirname ):
            if f in [ "orig", "sms.root", "validation" ]:
                continue
            if f[-1:]=="~":
                continue
            if f[0]==".":
                continue
            if f[-3:]==".py":
                continue
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

    def loadBinaryFile ( self, lastm_only = False ):
        """
        Load a binary database, returning last modified, file count, database.
        
        :param lastm_only: if true, the database itself is not read.
        :returns: database object, or None, if lastm_only == True.
        """
        if lastm_only and self.pcl_mtime[0]:
            ## doesnt need to load database, and mtime is already
            ## loaded
            return None

        if self.pcl_db:
            return self.pcl_db

        if not os.path.exists ( self.binfile ):
            return None

        try:
            with open ( self.binfile, "rb" ) as f:
                t0=time.time()
                self.pcl_python = serializer.load ( f )
                self.pcl_format_version = serializer.load ( f )
                self.pcl_mtime = serializer.load ( f )
                self._databaseVersion = serializer.load ( f )
                if not lastm_only:
                    if self.pcl_python != sys.version:
                        logger.warning ( "binary file was written with a different "
                                         "python version. Regenerating." )
                        self.createBinaryFile()
                        return self
                    if self.pcl_format_version != self.sw_format_version:
                        logger.warning ( "binary file format (%s) and format " 
                                         "supported by software (%s) disagree." % 
                           ( self.pcl_format_version, self.sw_format_version ) )
                        logger.warning ( "will recreate binary." )
                        self.createBinaryFile()
                        return self

                    logger.info ( "loading binary db file %s format version %s" % 
                            ( self.binfile, self.pcl_format_version ) )
                    self.hasFastLim = serializer.load ( f )
                    self.expResultList = serializer.load ( f )
                    t1=time.time()-t0
                    logger.info ( "Loaded database from %s in %.1f secs." % \
                            ( self.binfile, t1 ) )
        except EOFError as e:
            os.unlink ( self.binfile )
            if lastm_only:
                self.pcl_format_version = -1
                self.pcl_mtime = 0
                return self
            logger.error ( "%s is not a binary database file! recreate it!" % self.binfile )
            self.createBinaryFile()
        return self

    def checkBinaryFile ( self ):
        nu=self.needsUpdate()
        logger.debug ( "Checking binary db file." )
        logger.debug ( "Binary file dates to %s(%d)" % \
                      ( time.ctime(self.pcl_mtime[0]),self.pcl_mtime[1] ) )
        logger.debug ( "Database dates to %s(%d)" % \
                      ( time.ctime(self.txt_mtime[0]),self.txt_mtime[1] ) )
        if nu:
            logger.info ( "Binary db file needs an update." )
        else:
            logger.info ( "Binary db file does not need an update." )
        return nu

    def needsUpdate ( self ):
        """ does the binary db file need an update? """
        try:
            # logger.debug ( "needsUpdate?" )
            self.lastModifiedAndFileCount()
            self.loadBinaryFile ( lastm_only = True )
            return ( self.txt_mtime[0] > self.pcl_mtime[0] or \
                     self.txt_mtime[1] != self.pcl_mtime[1]  or \
                     self.sw_format_version != self.pcl_format_version
                )
        except (IOError,DatabaseNotFoundException,TypeError,ValueError):
            # if we encounter a problem, we rebuild the database.
            return True


    def createBinaryFile ( self, filename=None ):
        """ create a pcl file from the text database,
            potentially overwriting an old pcl file. """
        t0=time.time()
        logger.info ( "Creating binary database " )
        logger.info ( "(this may take a few minutes, but it's done only once!)" )
        logger.debug ( " * compute last modified timestamp." )
        self.lastModifiedAndFileCount()
        logger.debug (  " * compute timestamp: %s filecount: %d" % \
                ( time.ctime ( self.txt_mtime[0] ), self.txt_mtime[1] ) )
        binfile = filename
        if binfile == None:
            binfile = self.binfile
        logger.debug (  " * create %s" % self.binfile )
        with open ( binfile, "wb" ) as f:
            logger.debug (  " * load text database" )
            self.loadTextDatabase() 
            logger.debug (  " * write %s version %s" % ( self.binfile,
                       self.sw_format_version ) )
            ptcl = serializer.HIGHEST_PROTOCOL
            self.pcl_python = sys.version
            serializer.dump ( self.pcl_python, f, protocol=ptcl )
            serializer.dump ( self.sw_format_version, f, protocol=ptcl )
            serializer.dump ( self.txt_mtime, f, protocol=ptcl )
            serializer.dump ( self._databaseVersion, f, protocol=ptcl )
            serializer.dump ( self.hasFastLim, f, protocol=ptcl )
            serializer.dump ( self.expResultList, f, protocol=ptcl )
            logger.info (  " * done writing %s in %.1f secs." % \
                    ( binfile, time.time()-t0 ) )

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
        Raises an exception if something is wrong with the path.
    
        """
        logger.debug('Try to set the path for the database to: %s', path)
        tmp = os.path.realpath(path)
        if os.path.isfile ( tmp ):
            self._base = os.path.dirname ( tmp )
            self.force_load = "pcl" 
            self.pclfilename = os.path.basename ( tmp )
            return

        if tmp[-4:]==".pcl":
            if not os.path.exists ( tmp ):
                logger.error ( "File not found: %s" % tmp )
                sys.exit()
            logger.error ( "Supplied a pcl filename, but %s is not a file." % tmp )
            sys.exit()

        path = tmp + '/'
        if not os.path.exists(path):
            logger.error('%s is no valid path!' % path)
            raise DatabaseNotFoundException("Database not found")
        self._base = path

    def __str__(self):
        idList = "Database version: " + self.databaseVersion
        idList += "\n"
        idList += "-" * len(idList) + "\n"
        if self.expResultList == None:
            idList += "no experimental results available! "
            return idList
        idList += "%d experimental results: " % \
                   len ( self.expResultList ) 
        atlas,cms = [],[]
        datasets = 0
        txnames = 0
        s = { 8:0, 13:0  }
        for expRes in self.expResultList:
            Id = expRes.globalInfo.getInfo('id')
            sqrts = expRes.globalInfo.getInfo('sqrts').asNumber ( TeV )
            if not sqrts in s.keys():
                s[sqrts] = 0
            s[sqrts]+=1
            datasets += len ( expRes.datasets )
            for ds in expRes.datasets:
                txnames += len ( ds.txnameList )
            if "ATLAS" in Id:
                atlas.append ( expRes )
            if "CMS" in Id:
                cms.append ( expRes )
        idList += "%d CMS, %d ATLAS, " % ( len(cms), len(atlas) )
        for sqrts in s.keys():
            idList += "%d @ %d TeV, " % ( s[sqrts], sqrts )
            # idList += expRes.globalInfo.getInfo('id') + ', '
        idList = idList[:-2] + '\n'
        idList += "%d datasets, %d txnames.\n" % ( datasets, txnames )
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

    def _loadExpResults(self):
        """
        Checks the database folder and generates a list of ExpResult objects for
        each (globalInfo.txt,sms.py) pair.
       
        :returns: list of ExpResult objects 
  
        """
        folders=[]
        for root, _, files in os.walk(self._base):
            folders.append ( (root, files) )
        folders.sort()

        roots = []
        for root,files in folders:
            if "/.git/" in root:
                continue
            if root[-11:] == "/validation":
                continue
            if root[-5:] == "/orig":
                continue
            if not 'globalInfo.txt' in files:
          #      logger.debug("Missing globalInfo.txt in %s", root)
                continue
            else:
                roots.append ( root )

        resultsList = []
        for root in roots:
            expres = ExpResult(root)
            if expres:
                resultsList.append(expres)
                contact = expres.globalInfo.getInfo("contact")
                if contact and "fastlim" in contact.lower():
                    self.hasFastLim = True

        if not resultsList:
            logger.warning("Zero results loaded.")

        return resultsList


    def getExpResults(self, analysisIDs=['all'], datasetIDs=['all'], txnames=['all'],
                    dataTypes = ['all'], useSuperseded=False, useNonValidated=False):
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
        :param useSuperseded: If False, the supersededBy results will not be included
        :param useNonValidated: If False, the results with validated = False 
                                will not be included
        :returns: list of ExpResult objects or the ExpResult object if the list
                  contains only one result
                   
        """
        expResultList = []
        for expResult in self.expResultList:
            superseded = None
            if hasattr(expResult.globalInfo,'supersededBy'):
                superseded = expResult.globalInfo.supersededBy.replace(" ","")
            if superseded and (not useSuperseded):
                continue
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
                newDataSet = datasetObj.DataSet(dataset.path, dataset.globalInfo,False)
                newDataSet.dataInfo = dataset.dataInfo
                newDataSet.txnameList = []
                for txname in dataset.txnameList:
                    if txname.validated is False and (not useNonValidated):
                        continue
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

    def updateBinaryFile ( self ):
        """ write a binar db file, but only if 
            necessary. """
        if self.needsUpdate():
            logger.debug ( "Binary db file needs an update." )
            self.createBinaryFile()
        else:
            logger.debug ( "Binary db file does not need an update." )

class ExpResultList(object):
    """
    Holds a list of ExpResult objects for printout.
    
    :ivar expResultList: list of ExpResult objects 
        
    """
    def __init__(self, expResList):
        self.expResultList = expResList

if __name__ == "__main__":
    import argparse
    """ Run as a script, this checks and/or writes database.pcl files """
    argparser = argparse.ArgumentParser(description='simple script to check \
            and/or write database.pcl files')
    argparser.add_argument('-c', '--check', help='check binary db file',
                           action='store_true')
    argparser.add_argument('-t', '--time', help='time reading db',
                           action='store_true')
    argparser.add_argument('-r', '--read', help='read binary db file',
                           action='store_true')
    argparser.add_argument('-w', '--write', help='force writing binary db file',
                           action='store_true')
    argparser.add_argument('-u', '--update', help='update binary db file, if necessary',
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
        db.createBinaryFile()
        sys.exit()
    db = Database ( args.database )
    if args.debug:
        db.verbosity = "debug"
    logger.debug ( "%s" % db )
    if args.update:
        db.updateBinaryFile()
    if args.check:
        db.checkBinaryFile()
    if args.time:
        t0=time.time()
        expResult = db.loadBinaryFile ( lastm_only = False )
        t1=time.time()
        print ( "Time it took reading binary db file: %.1f s." % (t1-t0) )
        txtdb = db.loadTextDatabase()
        t2=time.time()
        print ( "Time it took reading text   file: %.1f s." % (t2-t1) )
    if args.read:
        db = db.loadBinaryFile ( lastm_only = False )
        listOfExpRes = db.getExpResults() 
        for expResult in listOfExpRes:
            print (expResult)

