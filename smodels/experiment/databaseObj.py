#!/usr/bin/env python

"""
.. module:: databaseObj
   :synopsis: Contains Database class that represents the database of experimental results.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
import os
## sweet spot for numpy multi-threading is 2? More threads
## make some weaker machines freeze when building the pickle file.
## Anyhow, we parallelize at a higher level.
os.environ["OMP_NUM_THREADS"]="2"
import sys
import time
from smodels.experiment import datasetObj
from smodels.experiment.metaObj import Meta
from smodels.experiment.expResultObj import ExpResult
from smodels.experiment.exceptions import DatabaseNotFoundException
from smodels.tools.physicsUnits import TeV
from smodels.tools.smodelsLogging import logger
import logging

try:
    import cPickle as serializer
except ImportError as e:
    import pickle as serializer

class Database(object):
    """
    Database object. Holds a list of ExpResult objects.

    :ivar base: path to the database (string)
    :ivar force_load: force loading the text database ("txt"),
        or binary database ("pcl"), dont force anything if None
    :ivar expResultList: list of ExpResult objects

    """

    def __init__( self, base=None, force_load = None, discard_zeroes = True,
                  progressbar = False, subpickle = False ):
        """
        :param base: path to the database, or pickle file (string)
        :param force_load: force loading the text database ("txt"),
            or binary database ("pcl"), dont force anything if None
        :param discard_zeroes: discard txnames with only zeroes as entries.
        :param progressbar: show a progressbar when building pickle file
                            (needs the python-progressbar module)
        :param subpickle: produce small pickle files per exp result.
            Should only be used when working on the database.
        """
        self.source=None
        self.force_load = force_load
        self.subpickle = subpickle
        base, pclfile = self.checkPathName(base, discard_zeroes )
        self.pcl_meta = Meta( pclfile )
        self.expResultList = []
        self.txt_meta = self.pcl_meta
        if not self.force_load == "pcl":
            self.txt_meta = Meta ( base, discard_zeroes = discard_zeroes )
        self.progressbar = None
        if progressbar:
            try:
                import progressbar as P
                self.progressbar = P.ProgressBar( widgets= 
                        [ "Building Database ", P.Percentage(), 
                          P.Bar( marker=P.RotatingMarker() ), P.ETA() ] )
            except ImportError as e:
                logger.warning ( "progressbar requested, but python-progressbar is not installed." )

        if self.force_load=="txt":
            self.loadTextDatabase()
            self.txt_meta.printFastlimBanner()
            return
        if self.force_load=="pcl":
            self.loadBinaryFile()
            self.pcl_meta.printFastlimBanner()
            return
        if self.force_load in [ None, "none", "None" ]:
            self.loadDatabase()
            self.txt_meta.printFastlimBanner()
            return
        logger.error ( "when initialising database: force_load=%s is not " \
                       "recognized. Valid values are: pcl, txt, None." % force_load )
        sys.exit()

    def __eq__ ( self, other ):
        """ compare two databases """
        if type ( self ) != type ( other ):
            return False
        if not self.txt_meta.sameAs ( other.txt_meta ):
            return False
        if len( self.expResultList ) != len (other.expResultList):
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
        if not os.path.exists ( self.pcl_meta.pathname ):
            logger.info ( "Creating binary database " )
            logger.info ( "(this may take a few minutes, but it's done only once!)" )
            self.loadTextDatabase()
            self.createBinaryFile()
        else:
            if self.needsUpdate():
                self.createBinaryFile()
            else:
                self.loadBinaryFile( lastm_only = False )

    def loadTextDatabase ( self ):
        """ simply loads the textdabase """
        if self.txt_meta.databaseVersion and len(self.expResultList)>0:
            logger.debug ( "Asked to load database, but has already been loaded. Ignore." )
            return
        logger.info ( "Parsing text database at %s" % self.txt_meta.pathname )
        self.expResultList = self._loadExpResults()

    def loadBinaryFile ( self, lastm_only = False ):
        """
        Load a binary database, returning last modified, file count, database.

        :param lastm_only: if true, the database itself is not read.
        :returns: database object, or None, if lastm_only == True.
        """
        if lastm_only and self.pcl_meta.mtime:
            ## doesnt need to load database, and mtime is already
            ## loaded
            return None

        if not os.path.exists ( self.pcl_meta.pathname ):
            return None

        try:
            with open ( self.pcl_meta.pathname, "rb" ) as f:
                t0=time.time()
                pclfilename = self.pcl_meta.pathname
                self.pcl_meta = serializer.load ( f )
                self.pcl_meta.pathname = pclfilename
                if self.force_load == "pcl":
                    self.txt_meta = self.pcl_meta
                if not lastm_only:
                    if not self.force_load == "pcl" and self.pcl_meta.needsUpdate ( self.txt_meta ):
                        logger.warning ( "Something changed in the environment."
                                         "Regenerating." )
                        self.createBinaryFile()
                        return self
                    logger.info ( "loading binary db file %s format version %s" %
                            ( self.pcl_meta.pathname, self.pcl_meta.format_version ) )
                    if sys.version[0]=="2":
                        self.expResultList = serializer.load ( f )
                    else:
                        self.expResultList = serializer.load ( f, encoding="latin1" )
                    t1=time.time()-t0
                    logger.info ( "Loaded database from %s in %.1f secs." % \
                            ( self.pcl_meta.pathname, t1 ) )
        except (EOFError,ValueError) as e:
            os.unlink ( self.pcl_meta.pathname )
            if lastm_only:
                self.pcl_meta.format_version = -1
                self.pcl_meta.mtime = 0
                return self
            logger.error ( "%s is not a binary database file, recreate it." % \
                            ( self.pcl_meta.pathname ) )
            self.createBinaryFile()
        # self.txt_meta = self.pcl_meta
        return self

    def checkBinaryFile ( self ):
        nu=self.needsUpdate()
        logger.debug ( "Checking binary db file." )
        logger.debug ( "Binary file dates to %s(%d)" % \
                      ( time.ctime(self.pcl_meta.mtime),self.pcl_meta.filecount ) )
        logger.debug ( "Database dates to %s(%d)" % \
                      ( time.ctime(self.txt_meta.mtime),self.txt_meta.filecount ) )
        if nu:
            logger.info ( "Binary db file needs an update." )
        else:
            logger.info ( "Binary db file does not need an update." )
        return nu

    def needsUpdate ( self ):
        """ does the binary db file need an update? """
        try:
            self.loadBinaryFile ( lastm_only = True )
            # logger.error ( "needs update?" )
            return ( self.pcl_meta.needsUpdate ( self.txt_meta ) )
        except (IOError,DatabaseNotFoundException,TypeError,ValueError):
            # if we encounter a problem, we rebuild the database.
            return True


    def createBinaryFile ( self, filename=None ):
        """ create a pcl file from the text database,
            potentially overwriting an old pcl file. """
        if self.txt_meta == None:
            logger.error ( "Trying to create database pickle, but no txt_meta defined." )
            sys.exit()
        logger.debug ( "database timestamp: %s, filecount: %d" % \
                     ( time.ctime ( self.txt_meta.mtime ), self.txt_meta.filecount ) )
        binfile = filename
        if binfile == None:
            binfile = self.pcl_meta.pathname
        logger.debug (  " * create %s" % binfile )
        with open ( binfile, "wb" ) as f:
            logger.debug (  " * load text database" )
            self.loadTextDatabase()
            logger.debug (  " * write %s db version %s, format version %s, %s" % \
                    ( binfile, self.txt_meta.databaseVersion,
                      self.txt_meta.format_version, self.txt_meta.cTime() ) )
            ptcl = serializer.HIGHEST_PROTOCOL
            serializer.dump ( self.txt_meta, f, protocol=ptcl )
            serializer.dump ( self.expResultList, f, protocol=ptcl )
            logger.info (  "%s created." % ( binfile ) )

    @property
    def databaseVersion(self):
        """
        The version of the database, read from the 'version' file.

        """
        return self.txt_meta.databaseVersion


    @property
    def base(self):
        """
        This is the path to the base directory.
        """
        return self.txt_meta.pathname

    def fetchFromScratch ( self, path, store, discard_zeroes ):
        """ fetch database from scratch, together with
            description.
            :param store: filename to store json file.
        """
        def sizeof_fmt(num, suffix='B'):
            for unit in [ '','K','M','G','T','P' ]:
                if abs(num) < 1024.:
                    return "%3.1f%s%s" % (num, unit, suffix)
                num /= 1024.0
            return "%.1f%s%s" % (num, 'Yi', suffix)

        import requests
        r = requests.get( path )
        if r.status_code != 200:
            logger.error ( "Error %d: could not fetch %s from server." % \
                           ( r.status_code, path ) )
            sys.exit()
        ## its new so store the description
        with open( store, "w" ) as f:
            f.write ( r.text )
        if not "url" in r.json().keys():
            logger.error ( "cannot parse json file %s." % path )
            sys.exit()
        size = r.json()["size"]
        logger.info ( "need to fetch %s. size is %s." % \
                      ( r.json()["url"], sizeof_fmt ( size ) ) )
        t0=time.time()
        r2=requests.get ( r.json()["url"], stream=True )
        filename= "./" + r2.url.split("/")[-1]
        with open ( filename, "wb" ) as dump:
            print ( "         " + " "*51 + "<", end="\r" )
            print ( "loading >", end="" )
            for x in r2.iter_content(chunk_size=int ( size / 50 ) ):
                dump.write ( x )
                dump.flush ()
                print ( ".", end="" )
                sys.stdout.flush()
            print()
            dump.close()
        logger.info ( "fetched %s in %d secs." % ( r2.url, time.time()-t0 ) ) 
        logger.debug ( "store as %s" % filename )
        #with open( filename, "wb" ) as f:
        #    f.write ( r2.content )
        #    f.close()
        self.force_load = "pcl"
        return ( "./", "%s" % filename )


    def fetchFromServer ( self, path, discard_zeroes ):
        import requests, time, json
        logger.debug ( "need fetch from server: %s" % path )
        store = "." + path.replace ( ":","_" ).replace( "/", "_" ).replace(".","_" )
        if not os.path.isfile ( store ):
            ## completely new! fetch the description and the db!
            return self.fetchFromScratch ( path, store, discard_zeroes )
        with open(store,"r") as f:
            jsn = json.load(f)
        filename= "./" + jsn["url"].split("/")[-1]
        class _: ## pseudo class for pseudo requests
            def __init__ ( self ): self.status_code = -1
        r=_()
        try:
            r = requests.get( path )
        except Exception:
            pass
        if r.status_code != 200:
            logger.warning ( "Error %d: could not fetch %s from server." % \
                           ( r.status_code, path ) )
            if not os.path.isfile ( filename ):
                logger.error ( "Cant find a local copy of the pickle file. Exit." )
                sys.exit()
            logger.warning ( "I do however have a local copy of the file. I work with that." )
            self.force_load = "pcl"
            # next step: check the timestamps
            return ( "./", filename )

        if r.json()["lastchanged"] > jsn["lastchanged"]:
            ## has changed! redownload everything!
            return self.fetchFromScratch ( path, store, discard_zeroes )
        
        if not os.path.isfile ( filename ):
            return self.fetchFromScratch ( path, store, discard_zeroes )
        self.force_load = "pcl"
        # next step: check the timestamps
        return ( "./", filename )

    def checkPathName( self, path, discard_zeroes ):
        """
        checks the path name,
        returns the base directory and the pickle file name.
        If path starts with http or ftp, fetch the description file
        and the database.
        :returns: directory name, (full) pickle file name
        """
        logger.debug('Try to set the path for the database to: %s', path)
        if path.startswith( ( "http://", "https://", "ftp://" ) ):
            return self.fetchFromServer ( path, discard_zeroes )
        if path.startswith( ( "file://" ) ):
            path=path[7:]
            
        tmp = os.path.realpath(path)
        if os.path.isfile ( tmp ):
            base = os.path.dirname ( tmp )
            return ( base, tmp )

        if tmp[-4:]==".pcl":
            if not os.path.exists ( tmp ):
                if self.force_load == "pcl":
                    logger.error ( "File not found: %s" % tmp )
                    sys.exit()
                logger.info ( "File not found: %s. Will generate." % tmp )
                base = os.path.dirname ( tmp )
                return ( base, tmp )
            logger.error ( "Supplied a pcl filename, but %s is not a file." % tmp )
            sys.exit()

        path = tmp + '/'
        if not os.path.exists(path):
            logger.error('%s is no valid path!' % path)
            raise DatabaseNotFoundException("Database not found")
        m=Meta ( path, discard_zeroes = discard_zeroes )
        return ( path, path + m.getPickleFileName() )

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

    def _loadExpResults(self):
        """
        Checks the database folder and generates a list of ExpResult objects for
        each (globalInfo.txt,sms.py) pair.

        :returns: list of ExpResult objects

        """
        folders=[]
        for root, _, files in os.walk(self.txt_meta.pathname):
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
                continue
            else:
                roots.append ( root )

        if self.progressbar:
            self.progressbar.maxval = len ( roots )
            self.progressbar.start()
        resultsList = []
        for ctr,root in enumerate(roots):
            if self.progressbar:
                self.progressbar.update(ctr)
            expres = self.createExpResult ( root )
            if expres:
                resultsList.append(expres)

        if not resultsList:
            logger.warning("Zero results loaded.")
        if self.progressbar:
            self.progressbar.finish()

        return resultsList

    def createExpResult ( self, root ):
        """ create, from pickle file or text files """
        txtmeta = Meta ( root, discard_zeroes = self.txt_meta.discard_zeroes, 
                         hasFastLim=None, databaseVersion = self.databaseVersion )
        pclfile = "%s/.%s" % ( root, txtmeta.getPickleFileName() )
        logger.debug ( "Creating %s, pcl=%s" % (root,pclfile ) )
        expres = None
        try:
            # logger.info ( "%s exists? %d" % ( pclfile,os.path.exists ( pclfile ) ) )
            if not self.force_load=="txt" and os.path.exists ( pclfile ):
                # logger.info ( "%s exists" % ( pclfile ) )
                with open(pclfile,"rb" ) as f:
                    logger.debug ( "Loading: %s" % pclfile )
                    ## read meta from pickle
                    pclmeta = serializer.load ( f )
                    if not pclmeta.needsUpdate ( txtmeta ):
                        logger.debug ( "we can use expres from pickle file %s" % pclfile )
                        expres = serializer.load ( f )
                    else:
                        logger.debug ( "we cannot use expres from pickle file %s" % pclfile )
                        logger.debug ( "txt meta %s" % txtmeta )
                        logger.debug ( "pcl meta %s" % pclmeta )
                        logger.debug ( "pcl meta needs update %s" % pclmeta.needsUpdate ( txtmeta ) )
        except IOError as e:
            logger.error ( "exception %s" % e )
        if not expres: ## create from text file
            expres = ExpResult(root, discard_zeroes = self.txt_meta.discard_zeroes )
            if self.subpickle and expres: expres.writePickle( self.databaseVersion )
        if expres:
            contact = expres.globalInfo.getInfo("contact")
            if contact and "fastlim" in contact.lower():
                self.txt_meta.hasFastLim = True
        return expres

    def getExpResults(self, analysisIDs=['all'], datasetIDs=['all'], txnames=['all'],
                    dataTypes = ['all'], useSuperseded=False, useNonValidated=False,
                    onlyWithExpected = False ):
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
        :param onlyWithExpected: Return only those results that have expected values
                 also. Note that this is trivially fulfilled for all efficiency maps.
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
                newDataSet = datasetObj.DataSet( dataset.path, dataset.globalInfo,
                       False, discard_zeroes=self.txt_meta.discard_zeroes )
                newDataSet.dataInfo = dataset.dataInfo
                newDataSet.txnameList = []
                for txname in dataset.txnameList:
                    if type(txname.validated) == str:
                        txname.validated = txname.validated.lower()
                    # print ( "txname",txname.validated,type(txname.validated) )
                    if (txname.validated not in [True, False, "true", "false", "n/a", "tbd", None, "none"]):
                        logger.error("value of validated field '%s' in %s unknown." % (txname.validated, expResult))
                    if txname.validated in [None, "none"]: ## FIXME after 1.1.1 this becomes a warning msg?
                        logger.debug("validated is None in %s/%s/%s. Please set to True, False, N/A, or tbd." % \
                            ( expResult.globalInfo.id, dataset.dataInfo.dataId, txname ) )
                    if txname.validated not in [ None, True, "true", "n/a", "tbd" ] and (not useNonValidated ):
#                    if txname.validated is False and (not useNonValidated):
                        continue
                    if txnames != ['all']:
                        if not txname.txName in txnames:
                            continue
                    if onlyWithExpected and dataset.dataInfo.dataType == \
                        "upperLimit" and not txname.txnameDataExp:
                        continue
                    newDataSet.txnameList.append(txname)
                # Skip data set not containing any of the required txnames:
                if not newDataSet.txnameList or newDataSet.txnameList == []:
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
    from smodels.tools.smodelsLogging import setLogLevel
    """ Run as a script, this checks and/or writes dbX.pcl files """
    argparser = argparse.ArgumentParser(description='simple script to check \
            and/or write dbX.pcl files')
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
        setLogLevel(level=logging.DEBUG )
    if args.write:
        db = Database ( args.database, force_load="txt" )
        db.createBinaryFile()
        sys.exit()
    db = Database ( args.database )
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

