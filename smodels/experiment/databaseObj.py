#!/usr/bin/env python3

"""
.. module:: databaseObj
   :synopsis: Contains Database class that represents the database of experimental results.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Matthias Wolf <matthias.wolf@wot.at>

"""

from __future__ import print_function
import os
import hashlib
import pathlib
## sweet spot for numpy multi-threading is 2? More threads
## make some weaker machines freeze when building the pickle file.
## Anyhow, we parallelize at a higher level.
os.environ["OMP_NUM_THREADS"]="2"
import sys
import time
import copy
from smodels.experiment import datasetObj
from smodels.installation import cacheDirectory
from smodels.experiment.metaObj import Meta
from smodels.experiment.expResultObj import ExpResult
from smodels.experiment.exceptions import DatabaseNotFoundException
from smodels.tools.physicsUnits import TeV
from smodels.tools.stringTools import cleanWalk
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from smodels.tools.smodelsLogging import logger
import logging

try:
    import cPickle as serializer
except ImportError as e:
    import pickle as serializer

def _getSHA1 ( filename ):
    return hashlib.sha1( pathlib.Path(filename).read_bytes() ).hexdigest()

class Database(object):
    """ Database object. Holds a list of SubDatabases.
        Delegates all calls to SubDatabases.
    """
    def __init__(self, base=None, force_load = None, discard_zeroes = True,
                  progressbar = False, subpickle = True, combinationsmatrix = None ):
        """
        :param base: path to the database, or pickle file (string), or http
                     address. If None, "official", or "official_fastlim",
                     use the official database for your code version
                     (including fastlim results, if specified).
                     If "latest", or "latest_fastlim", check for the latest database.
                     Multiple databases may be specified using `+' as a delimiter.
        :param force_load: force loading the text database ("txt"),
                           or binary database ("pcl"), dont force anything if None
        :param discard_zeroes: discard txnames with only zeroes as entries.
        :param progressbar: show a progressbar when building pickle file
                            (needs the python-progressbar module)
        :param subpickle: produce small pickle files per exp result.
                          Should only be used when working on the database.
        :param combinationsmatrix: an optional dictionary that contains info
                     about combinable analyses, e.g. { "anaid1": ( "anaid2", "anaid3" ) }
                     optionally specifying signal regions, e.g. { "anaid1:SR1":
                     ( "anaid2:SR2", "anaid3" ) }
        """
        self.subs = []
        if "_fastlim" in base: ## for backwards compatibility
            base = base.replace("_fastlim","+fastlim")
        sstrings = base.split ( "+" )
        for ss in sstrings:
            self.subs.append ( SubDatabase ( ss, force_load, discard_zeroes,
                                             progressbar, subpickle, combinationsmatrix ) )

    @property
    def expResultList(self):
        """
        The combined list, compiled from the individual lists

        """
        if len(self.subs)==0:
            return []

        lists = [ x.expResultList for x in self.subs ]
        return self.mergeLists ( lists )

    def mergeLists ( self, lists ):
        """ small function, merges lists of ERs """
        D = {}
        for tmp in lists:
            for t in tmp:
                if len(t.datasets)== 0: # skip empty entries
                    logger.warning ( f"Analysis {t.globalInfo.id} has no datasets. Will remove it." )
                    continue
                anaid = t.globalInfo.id + t.datasets[0].getType()
                if not anaid in D:
                    D[anaid]=t
                else: ## FIXME merge expResults
                    D[anaid]=self.mergeERs ( D[anaid], t )
        return list ( D.values() )

    def mergeERs ( self, o1, r2 ):
        """ merge the content of exp res r1 and r2 """
        r1 = copy.deepcopy ( o1 )
        dids = [ x.getID() for x in o1.datasets ]
        for ds in r2.datasets:
            if not ds.getID() in dids: ## completely new dataset
                r1.datasets.append ( ds )
            else: ## just overwrite the old txnames
                idx = dids.index ( ds.getID() ) ## ds index
                r2txs = ds.txnameList
                r1txnames = [ x.txName for x in  r1.datasets[idx].txnameList ]
                for txn in r2txs:
                     if txn.txName in r1txnames:
                        tidx = r1txnames.index ( txn.txName ) ## overwrite
                        r1.datasets[idx].txnameList[tidx]=txn
                     else:
                        # a new txname
                        r1.datasets[idx].txnameList.append ( txn )
        return r1

    def createBinaryFile(self, filename=None ):
        """ create a pcl file from all the subs """
        ## make sure we have a model to pickle with the database!
        logger.debug(  " * create %s" % filename )
        if filename == None:
            filename = self.pcl_meta.pathname
        with open( filename, "wb" ) as f:
            logger.debug(  " * load text database" )
            logger.debug(  " * write %s db version %s" % \
                    ( filename, self.databaseVersion ) )
            ptcl = min ( 4, serializer.HIGHEST_PROTOCOL )
            ## 4 is default protocol in python3.8, and highest protocol in 3.7
            serializer.dump(self.txt_meta, f, protocol=ptcl)
            serializer.dump(self.expResultList, f, protocol=ptcl)
            serializer.dump(self.databaseParticles, f, protocol=ptcl )
            logger.info(  "%s created." % ( filename ) )

    def __str__(self):
        # r = [ str(x) for x in self.subs ]
        # return "+".join(r)
        idList = "Database version: " + self.databaseVersion
        idList += "\n"
        idList += "-" * len(idList) + "\n"
        if self.expResultList == None:
            idList += "no experimental results available! "
            return idList
        idList += "%d experimental results: " % \
                   len( self.expResultList )
        atlas,cms = [],[]
        datasets = 0
        txnames = 0
        s = { 8:0, 13:0  }
        for expRes in self.expResultList:
            Id = expRes.globalInfo.getInfo('id')
            sqrts = expRes.globalInfo.getInfo('sqrts').asNumber( TeV )
            if not sqrts in s.keys():
                s[sqrts] = 0
            s[sqrts]+=1
            datasets += len( expRes.datasets )
            for ds in expRes.datasets:
                txnames += len( ds.txnameList )
            if "ATLAS" in Id:
                atlas.append( expRes )
            if "CMS" in Id:
                cms.append( expRes )
        idList += "%d CMS, %d ATLAS, " % ( len(cms), len(atlas) )
        for sqrts in s.keys():
            idList += "%d @ %d TeV, " % ( s[sqrts], sqrts )
            # idList += expRes.globalInfo.getInfo('id') + ', '
        idList = idList[:-2] + '\n'
        idList += "%d datasets, %d txnames.\n" % ( datasets, txnames )
        return idList

    def __eq__( self, other ):
        if type(other) != type(self):
            return False
        for x,y in zip ( self.subs, other.subs ):
            if x != y:
                return False
        return True

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

        :param analysisIDs: list of analysis ids ([CMS-SUS-13-006,...]). Can
                            be wildcarded with usual shell wildcards: * ? [<letters>]
                            Furthermore, the centre-of-mass energy can be chosen
                            as suffix, e.g. ":13*TeV". Note that the asterisk
                            in the suffix is not a wildcard.
        :param datasetIDs: list of dataset ids ([ANA-CUT0,...]). Can be wildcarded
                            with usual shell wildcards: * ? [<letters>]
        :param txnames: list of txnames ([TChiWZ,...]). Can be wildcarded with
                            usual shell wildcards: * ? [<letters>]
        :param dataTypes: dataType of the analysis (all, efficiencyMap or upperLimit)
                            Can be wildcarded with usual shell wildcards: * ? [<letters>]
        :param useSuperseded: If False, the supersededBy results will not be included
                              (deprecated)
        :param useNonValidated: If False, the results with validated = False
                                will not be included
        :param onlyWithExpected: Return only those results that have expected values
                 also. Note that this is trivially fulfilled for all efficiency maps.
        :returns: list of ExpResult objects or the ExpResult object if the list
                  contains only one result

        """
        if useSuperseded:
            hasSuperseded = False
            for s in self.subs:
                if "superseded" in s.url:
                    hasSuperseded = True
                    break
            ss = ""
            if hasSuperseded:
                ss = " - which you seem to have already done"
            logger.warning ( "the useSuperseded flag is deprecated from smodels v2.1 onwards. if you wish to use superseded results, please simply add them to your database path%s, e.g. 'official+superseded'." % ss )
        ret = []
        for sub in self.subs:
            tmp = sub.getExpResults( analysisIDs, datasetIDs, txnames, dataTypes,
                    True, useNonValidated, onlyWithExpected )
            ret.append ( tmp )
        return self.mergeLists ( ret )

    @property
    def databaseParticles(self):
        """
        Database particles, a list, one entry per sub
        """
        r = [ x.databaseParticles for x in self.subs ]
        return r[0] ## FIXME do sth smarter?

    @property
    def databaseVersion(self):
        """
        The version of the database, concatenation of the individual versions

        """
        r = [ x.databaseVersion for x in self.subs ]
        for i,ri in enumerate(r): # avoid repetitions
            for j,rj in enumerate(r[i+1:]):
                if ri in rj:
                    r[i+j+1]=rj.replace(ri,"")

        return "+".join ( r )

    @property
    def txt_meta(self):
        """
        The meta info of the text version, a merger of the original ones

        """
        r = [ x.txt_meta for x in self.subs ]
        ret = r[0]
        return ret

    @property
    def pcl_meta(self):
        """
        The meta info of the text version, a merger of the original ones

        """
        ret = None
        r = []
        for x in self.subs:
            if hasattr ( x, "pcl_meta" ):
                r.append ( x.pcl_meta )
        ret = r[0]
        return ret

    def createLinksToCombinationsMatrix ( self ):
        """ in all globalInfo objects, create a shallow link to the
            combinations matrix """
        for x in self.subs:
            if not hasattr ( x, "combinationsmatrix" ) or x.combinationsmatrix == None:
                x.combinationsmatrix = self.combinationsmatrix
            x.createLinksToCombinationsMatrix()

    def clearLinksToCombinationsMatrix ( self ):
        """ clear all shallow links to the combinations matrix """
        self.combinationsmatrix = None
        for x in self.subs:
            x.combinationsmatrix = None
            x.clearLinksToCombinationsMatrix()

class SubDatabase(object):
    """
    SubDatabase object. Holds a list of ExpResult objects.
    """

    def __init__(self, base=None, force_load = None, discard_zeroes = True,
                  progressbar = False, subpickle = True, combinationsmatrix=None ):
        """
        :param base: path to the database, or pickle file (string), or http
                     address. If None, "official", or "official_fastlim",
                     use the official database for your code version
                     (including fastlim results, if specified).
                     If "latest", or "latest_fastlim", check for the latest database.
                     Multiple databases may be named, use "+" as delimiter.
                     Order matters: Results with same name will overwritten
                     according to sequence
        :param force_load: force loading the text database ("txt"),
                           or binary database ("pcl"), dont force anything if None
        :param discard_zeroes: discard txnames with only zeroes as entries.
        :param progressbar: show a progressbar when building pickle file
                            (needs the python-progressbar module)
        :param subpickle: produce small pickle files per exp result.
                          Should only be used when working on the database.
        :param combinationsmatrix: an optional dictionary that contains info
                     about combinable analyses, e.g. { "anaid1": ( "anaid2", "anaid3" ) }
                     optionally specifying signal regions, e.g. { "anaid1:SR1":
                     ( "anaid2:SR2", "anaid3" ) }
        """

        self.url = base
        self.combinationsmatrix = combinationsmatrix
        self.source=""
        if force_load == None and base.endswith(".pcl"):
            force_load = "pcl"
        self.force_load = force_load
        self.subpickle = subpickle
        obase = base ## keep old name for more checks for 'latest'
        from smodels.installation import __dblabels__
        if base in __dblabels__:
            from smodels.installation import databasePath
            base = databasePath( base )
        base, pclfile = self.checkPathName(base, discard_zeroes )
        self.pcl_meta = Meta( pclfile )
        self.expResultList = []
        self.txt_meta = self.pcl_meta
        if not self.force_load == "pcl":
            self.txt_meta = Meta( base, discard_zeroes = discard_zeroes )
        self.progressbar = None
        if progressbar:
            try:
                import progressbar as P
                self.progressbar = P.ProgressBar( widgets=
                        [ "Building Database ", P.Percentage(),
                          P.Bar( marker=P.RotatingMarker() ), P.ETA() ] )
            except ImportError as e:
                logger.warning( "progressbar requested, but python-progressbar is not installed." )

        if self.force_load=="txt":
            self._setParticles()
            self.loadTextDatabase()
            self.txt_meta.printFastlimBanner()
            return
        if self.force_load=="pcl":
            self.loadBinaryFile()
            self._setParticles()
            self.pcl_meta.printFastlimBanner()
            if "latest" in obase:
                from smodels import installation
                codeVersion = installation.version()
                pclVersion = self.pcl_meta.databaseVersion
                if codeVersion[0]!=pclVersion[0]:
                    logger.error ( "major versions of code and database differ! code=%s, database=%s" % ( codeVersion[0], pclVersion[0] ) )
            return
        if self.force_load in [ None, "none", "None" ]:
            self.loadDatabase()
            self._setParticles()
            self.txt_meta.printFastlimBanner()
            return
        logger.error( "when initialising database: force_load=%s is not " \
                       "recognized. Valid values are: pcl, txt, None." % force_load )
        raise SModelSError()

    def __eq__( self, other ):
        """ compare two databases """
        if type( self ) != type( other ):
            return False
        if not self.txt_meta.sameAs( other.txt_meta ):
            return False
        if len( self.expResultList ) != len(other.expResultList):
            return False
        for ( myres, otherres ) in zip( self.expResultList, other.expResultList ):
            if myres != otherres:
                return False
        return True

    def loadDatabase( self ):
        """ if no binary file is available, then
            load the database and create the binary file.
            if binary file is available, then check if
            it needs update, create new binary file, in
            case it does need an update.
        """
        if not os.path.exists( self.pcl_meta.pathname ):
            logger.info( "Creating binary database " )
            logger.info( "(this may take a few minutes, but it's done only once!)" )
            self.loadTextDatabase()
            self.createBinaryFile()
        else:
            if self.needsUpdate():
                self.createBinaryFile()
            else:
                self.loadBinaryFile( lastm_only = False )

    def loadTextDatabase( self ):
        """ simply loads the textdabase """
        if self.txt_meta.databaseVersion and len(self.expResultList)>0:
            logger.debug( "Asked to load database, but has already been loaded. Ignore." )
            return
        logger.info( "Parsing text database at %s" % self.txt_meta.pathname )
        self.expResultList = self._loadExpResults()
        self.createLinksToModel()
        self.createLinksToCombinationsMatrix()

    def createLinksToModel( self ):
        """ in all globalInfo objects, create links to self.databaseParticles """
        if not hasattr ( self, "databaseParticles" ):
            return
        if type(self.databaseParticles) == type(None):
            return
        for ctr,er in enumerate(self.expResultList):
            if not hasattr ( er.globalInfo, "_databaseParticles" ):
                er.globalInfo._databaseParticles = self.databaseParticles
            elif type(er.globalInfo._databaseParticles) == type(None):
                er.globalInfo._databaseParticles = self.databaseParticles

    def createLinksToCombinationsMatrix( self ):
        """ in all globalInfo objects, create links to self.combinationsmatrix  """
        if not hasattr ( self, "combinationsmatrix" ):
            return
        if type(self.combinationsmatrix) == type(None):
            return
        for ctr,er in enumerate(self.expResultList):
            if not hasattr ( er.globalInfo, "_combinationsmatrix" ):
                er.globalInfo._combinationsmatrix = self.combinationsmatrix
            elif type(er.globalInfo._combinationsmatrix) == type(None):
                er.globalInfo._combinationsmatrix = self.combinationsmatrix

    def clearLinksToCombinationsMatrix( self ):
        for ctr,er in enumerate(self.expResultList):
            if hasattr ( er.globalInfo, "_combinationsmatrix" ):
                del er.globalInfo._combinationsmatrix


    def removeLinksToModel ( self ):
        """ remove the links of globalInfo._databaseParticles to the model.
            Currently not used. """
        for ctr,er in enumerate(self.expResultList):
            if hasattr ( er.globalInfo, "_databaseParticles" ):
                del er.globalInfo._databaseParticles

    def loadBinaryFile( self, lastm_only = False ):
        """
        Load a binary database, returning last modified, file count, database.

        :param lastm_only: if true, the database itself is not read.
        :returns: database object, or None, if lastm_only == True.
        """
        if lastm_only and self.pcl_meta.mtime:
            ## doesnt need to load database, and mtime is already
            ## loaded
            return None

        if not os.path.exists( self.pcl_meta.pathname ):
            return None

        try:
            with open( self.pcl_meta.pathname, "rb" ) as f:
                t0=time.time()
                pclfilename = self.pcl_meta.pathname
                self.pcl_meta = serializer.load( f )
                self.pcl_meta.pathname = pclfilename
                if self.force_load == "pcl":
                    self.txt_meta = self.pcl_meta
                if not lastm_only:
                    if not self.force_load == "pcl" and self.pcl_meta.needsUpdate( self.txt_meta ):
                        logger.warning( "Something changed in the environment."
                                         "Regenerating." )
                        self.createBinaryFile()
                        return self
                    logger.info( "loading binary db file %s format version %s" %
                           ( self.pcl_meta.pathname, self.pcl_meta.format_version ) )
                    if sys.version[0]=="2":
                        self.expResultList = serializer.load( f )
                    else:
                        self.expResultList = serializer.load( f, encoding="latin1" )
                    t1=time.time()-t0
                    logger.info( "Loaded database from %s in %.1f secs." % \
                            ( self.pcl_meta.pathname, t1 ) )
                    self.databaseParticles = None
                    try:
                        self.databaseParticles = serializer.load ( f )
                    except EOFError as e:
                        pass ## a model does not *have* to be defined
                    self.createLinksToModel()
                    self.createLinksToCombinationsMatrix()
        except(EOFError,ValueError) as e:
            os.unlink( self.pcl_meta.pathname )
            if lastm_only:
                self.pcl_meta.format_version = -1
                self.pcl_meta.mtime = 0
                return self
            logger.error( "%s is not readable (%s)." % \
                            ( self.pcl_meta.pathname, str(e) ) )
            if self.source in [ "http", "ftp", "pcl" ]:
                logger.error( "source cannot be rebuilt. supply a different path to the database in your ini file." )
                raise SModelSError()
            self.createBinaryFile()
        # self.txt_meta = self.pcl_meta
        return self

    def checkBinaryFile( self ):
        nu=self.needsUpdate()
        logger.debug( "Checking binary db file." )
        logger.debug( "Binary file dates to %s(%d)" % \
                     ( time.ctime(self.pcl_meta.mtime),self.pcl_meta.filecount ) )
        logger.debug( "Database dates to %s(%d)" % \
                      ( time.ctime(self.txt_meta.mtime),self.txt_meta.filecount ) )
        if nu:
            logger.info( "Binary db file needs an update." )
        else:
            logger.info( "Binary db file does not need an update." )
        return nu

    def needsUpdate( self ):
        """ does the binary db file need an update? """
        try:
            self.loadBinaryFile( lastm_only = True )
            # logger.error( "needs update?" )
            return( self.pcl_meta.needsUpdate( self.txt_meta ) )
        except(IOError,DatabaseNotFoundException,TypeError,ValueError):
            # if we encounter a problem, we rebuild the database.
            return True

    def createBinaryFile(self, filename=None):
        """ create a pcl file from the text database,
            potentially overwriting an old pcl file. """
        ## make sure we have a model to pickle with the database!
        if self.txt_meta == None:
            logger.error("Trying to create database pickle, but no txt_meta defined." )
            raise SModelSError()
        logger.debug( "database timestamp: %s, filecount: %s" % \
                     ( time.ctime( self.txt_meta.mtime ), self.txt_meta.filecount ) )
        binfile = filename
        if binfile == None:
            binfile = self.pcl_meta.pathname
        if not hasattr(self,'databaseParticles') or \
            type(self.databaseParticles) == type(None):
           self._setParticles(self._getParticles())
        logger.debug(  " * create %s" % binfile )
        with open( binfile, "wb" ) as f:
            logger.debug(  " * load text database" )
            self.loadTextDatabase()
            logger.debug(  " * write %s db version %s, format version %s, %s" % \
                    ( binfile, self.txt_meta.databaseVersion,
                      self.txt_meta.format_version, self.txt_meta.cTime() ) )
            # ptcl = serializer.HIGHEST_PROTOCOL
            ptcl = min ( 4, serializer.HIGHEST_PROTOCOL ) ## 4 is default protocol in python3.8, and highest protocol in 3.7
            serializer.dump(self.txt_meta, f, protocol=ptcl)
            serializer.dump(self.expResultList, f, protocol=ptcl)
            serializer.dump(self.databaseParticles, f, protocol=ptcl )
            logger.info(  "%s created." % ( binfile ) )

    @property
    def databaseVersion(self):
        """
        The version of the database, read from the 'version' file.

        """
        return self.txt_meta.databaseVersion

    @databaseVersion.setter
    def databaseVersion(self, x ):
        self.txt_meta.databaseVersion = x
        self.pcl_meta.databaseVersion = x

    def inNotebook(self):
        """
        Are we running within a notebook? Has an effect on the
        progressbar we wish to use.
        """

        try:
            cfg = get_ipython().config
            if 'IPKernelApp' in cfg.keys():
                return True
            else:
                return False
        except NameError:
            return False

    @property
    def base(self):
        """
        This is the path to the base directory.
        """
        return self.txt_meta.pathname

    def fetchFromScratch( self, path, store, discard_zeroes ):
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
        try:
            r = requests.get( path, timeout=5 )
        except requests.exceptions.RequestException as e:
            logger.error( "Exception when trying to fetch database: %s" % e )
            logger.error( "Consider supplying a different database path in the ini file (possibly a local one)" )
            raise SModelSError()
        if r.status_code != 200:
            line = "Error %d: could not fetch '%s' from server: '%s'" % \
                           ( r.status_code, path, r.reason )
            logger.error( line )
            raise SModelSError( line )
        ## its new so store the description
        with open( store, "w" ) as f:
            f.write( r.text )
        if not "url" in r.json().keys():
            logger.error( "cannot parse json file %s." % path )
            raise SModelSError()
        size = r.json()["size"]
        cDir,defused = cacheDirectory ( create=True, reportIfDefault=True )
        t0=time.time()
        filename =  os.path.join ( cDir, r.json()["url"].split("/")[-1] )
        if os.path.exists ( filename ):
            # if file exists and checksums match, we dont download
            if "sha1" in r.json():
                sha = _getSHA1 ( filename )
                if sha == r.json()["sha1"]:
                    ## seems it hasnt changed
                    self.force_load = "pcl"
                    return ( "./", "%s" % filename )
        r2=requests.get ( r.json()["url"], stream=True, timeout=(500,500) )
        # filename= os.path.join ( cDir, r2.url.split("/")[-1] )
        msg = "caching the downloaded database in %s." % cDir
        if defused:
            msg += " If you want the pickled database file to be cached in a different location, set the environment variable SMODELS_CACHEDIR, e.g. to '/tmp'."
        logger.warning ( msg )
        logger.info ( "need to fetch %s and store in %s. size is %s." % \
                      ( r.json()["url"], filename, sizeof_fmt ( size ) ) )
        with open( filename, "wb" ) as dump:
            import fcntl
            fcntl.lockf ( dump, fcntl.LOCK_EX )
            if not self.inNotebook(): ## \r doesnt work in notebook
                print( "         " + " "*51 + "<", end="\r" )
            print( "loading >", end="" )
            for x in r2.iter_content(chunk_size=int( size / 50 ) ):
                dump.write( x )
                dump.flush()
                print( ".", end="" )
                sys.stdout.flush()
            if self.inNotebook():
                print( "done." )
            else:
                print( "" )
            fcntl.lockf ( dump, fcntl.LOCK_UN )
            dump.close()
        logger.info( "fetched %s in %d secs." % ( r2.url, time.time()-t0 ) )
        logger.debug( "store as %s" % filename )
        self.force_load = "pcl"
        return ( "./", "%s" % filename )

    def fetchFromServer( self, path, discard_zeroes ):
        import requests, time, json
        self.source = "http"
        if "ftp://" in path:
            self.source = "ftp"
        cDir = cacheDirectory ( create=True )
        store = os.path.join ( cDir, path.replace ( ":","_" ).replace( "/", "_" ).replace(".","_" ) )
        logger.debug ( "need to fetch from server: %s and store to %s" % ( path, store ) )
        if not os.path.isfile( store ):
            ## completely new! fetch the description and the db!
            return self.fetchFromScratch( path, store, discard_zeroes )
        with open(store,"r") as f:
            jsn = json.load(f)
        filename= os.path.join ( cDir, jsn["url"].split("/")[-1] )
        class _: ## pseudo class for pseudo requests
            def __init__( self ): self.status_code = -1
        r=_()
        try:
            r = requests.get( path, timeout=2 )
        except requests.exceptions.RequestException as e:
            pass
        if r.status_code != 200:
            logger.warning( "Error %d: could not fetch %s from server." % \
                           ( r.status_code, path ) )
            if not os.path.isfile( filename ):
                logger.error( "Cant find a local copy of the pickle file. Exit." )
                sys.exit()
            logger.warning ( "I do however have a local copy of the file at %s. I work with that." % filename )
            self.force_load = "pcl"
            return ( cDir, filename )
            #return ( cDir, os.path.basename ( filename ) )

        if not os.path.exists ( filename ):
            return self.fetchFromScratch ( path, store, discard_zeroes )
        stats = os.stat ( filename )
        if abs ( stats.st_size - jsn["size"]) > 4096:
            ## size doesnt match (4096 is to allow for slightly different file
            ## sizes reported by the OS). redownload!
            return self.fetchFromScratch ( path, store, discard_zeroes )
        """
        # dont do this b/c its slowish
        if "sha1" in r.json():
            t0 = time.time()
            sha = _getSHA1 ( filename )
            print ( "it took", time.time()-t0 )
            if sha != r.json()["sha1"]:
                return self.fetchFromScratch ( path, store, discard_zeroes )
        """
        if r.json()["lastchanged"] > jsn["lastchanged"]:
            ## has changed! redownload everything!
            return self.fetchFromScratch( path, store, discard_zeroes )

        if not os.path.isfile( filename ):
            return self.fetchFromScratch( path, store, discard_zeroes )
        self.force_load = "pcl"
        return ( "./", filename )

    def checkPathName( self, path, discard_zeroes ):
        """
        checks the path name,
        returns the base directory and the pickle file name.
        If path starts with http or ftp, fetch the description file
        and the database.
        returns the base directory and the pickle file name
        """
        logger.debug('Try to set the path for the database to: %s', path)
        if path.startswith( ( "http://", "https://", "ftp://" ) ):
            return self.fetchFromServer( path, discard_zeroes )
        if path.startswith( ( "file://" ) ):
            path=path[7:]

        tmp = os.path.realpath(path)
        if os.path.isfile( tmp ):
            base = os.path.dirname( tmp )
            return ( base, tmp )

        if tmp[-4:]==".pcl":
            self.source="pcl"
            if not os.path.exists( tmp ):
                if self.force_load == "pcl":
                    logger.error( "File not found: %s" % tmp )
                    raise SModelSError()
                logger.info( "File not found: %s. Will generate." % tmp )
                base = os.path.dirname( tmp )
                return ( base, tmp )
            logger.error( "Supplied a pcl filename, but %s is not a file." % tmp )
            raise SModelSError()

        path = tmp + '/'
        if not os.path.exists(path):
            logger.error('%s is no valid path!' % path)
            raise DatabaseNotFoundException("Database not found")
        m=Meta( path, discard_zeroes = discard_zeroes )
        self.source="txt"
        return ( path, path + m.getPickleFileName() )

    def __str__(self):
        idList = "Database version: " + self.databaseVersion
        idList += "\n"
        idList += "-" * len(idList) + "\n"
        if self.expResultList == None:
            idList += "no experimental results available! "
            return idList
        idList += "%d experimental results: " % \
                   len( self.expResultList )
        atlas,cms = [],[]
        datasets = 0
        txnames = 0
        s = { 8:0, 13:0  }
        for expRes in self.expResultList:
            Id = expRes.globalInfo.getInfo('id')
            sqrts = expRes.globalInfo.getInfo('sqrts').asNumber( TeV )
            if not sqrts in s.keys():
                s[sqrts] = 0
            s[sqrts]+=1
            datasets += len( expRes.datasets )
            for ds in expRes.datasets:
                txnames += len( ds.txnameList )
            if "ATLAS" in Id:
                atlas.append( expRes )
            if "CMS" in Id:
                cms.append( expRes )
        idList += "%d CMS, %d ATLAS, " % ( len(cms), len(atlas) )
        for sqrts in s.keys():
            idList += "%d @ %d TeV, " % ( s[sqrts], sqrts )
            # idList += expRes.globalInfo.getInfo('id') + ', '
        idList = idList[:-2] + '\n'
        idList += "%d datasets, %d txnames.\n" % ( datasets, txnames )
        return idList

    def _setParticles(self,databaseParticles=None):
        """
        Set the databaseParticles attribute.

        If databaseParticles is None and the self.databaseParticles is None,
        try to use the particles stored in the first ExpResult
        in the database (ExptResult.globalInfo._databaseParticles).
        If not found, fallback to the final states defined in defaultFinalStates.py.
        :param databaseParticles: Model object containing the final state particles
                                  used in the database.
        """
        #If not yet defined, set the attribute to None:
        if not hasattr(self,'databaseParticles'):
            self.databaseParticles = None
        #If input is given, use it to set the databaseParticles attribute:
        if databaseParticles:
            logger.debug("Setting database particles from %s" %str(databaseParticles))
            self.databaseParticles = databaseParticles

        #If still None, fallback to default:
        if self.databaseParticles is None:
            logging.debug("databaseParticles not found. Using default state.")
            from smodels.experiment.defaultFinalStates import finalStates
            self.databaseParticles = finalStates

    def _getParticles(self, particlesFile='databaseParticles.py'):
        """
        Load the particle objects used in the database.

        The particles are searched for in the database folder.
        If not found, the default particles will be loaded.
        """
        fulldir = os.path.join(self.txt_meta.pathname,particlesFile)
        if os.path.isfile(fulldir):
            from importlib import import_module
            sys.path.append(self.txt_meta.pathname)
            pFile = os.path.splitext(particlesFile)[0]
            logger.debug("Loading database particles from: %s" %fulldir)
            modelFile = import_module(pFile, package='smodels')
            if not hasattr(modelFile,'finalStates'):
                logger.error("Model definition (finalStates) not found in" % fulldir)
            else:
                #set model name to file location:
                modelFile.finalStates.label = os.path.basename(fulldir)
                return modelFile.finalStates

        return None

    def _loadExpResults(self):
        """
        Checks the database folder and generates a list of ExpResult objects for
        each (globalInfo.txt,sms.py) pair.

        :returns: list of ExpResult objects
        """
        #Try to load particles from databaseParticles.py
        self._setParticles(self._getParticles())
        folders=[]
        #for root, _, files in os.walk(self.txt_meta.pathname):
        # for root, _, files in cleanWalk(self._base):
        for root, _, files in cleanWalk(self.txt_meta.pathname):
            folders.append( (root, files) )
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
                roots.append( root )

        if self.progressbar:
            self.progressbar.maxval = len( roots )
            self.progressbar.start()
        resultsList = []
        for ctr,root in enumerate(roots):
            if self.progressbar:
                self.progressbar.update(ctr)
            expres = self.createExpResult( root )
            if expres:
                resultsList.append(expres)

        if not resultsList:
            logger.warning("Zero results loaded.")
        if self.progressbar:
            self.progressbar.finish()

        return resultsList

    def createExpResult( self, root ):
        """ create, from pickle file or text files """
        txtmeta = Meta( root, discard_zeroes = self.txt_meta.discard_zeroes,
                         hasFastLim=None, databaseVersion = self.databaseVersion )
        pclfile = "%s/.%s" % ( root, txtmeta.getPickleFileName() )
        logger.debug( "Creating %s, pcl=%s" % (root,pclfile ) )
        expres = None
        try:
            # logger.info( "%s exists? %d" % ( pclfile,os.path.exists( pclfile ) ) )
            if not self.force_load=="txt" and os.path.exists( pclfile ):
                # logger.info( "%s exists" % ( pclfile ) )
                with open(pclfile,"rb" ) as f:
                    logger.debug( "Loading: %s" % pclfile )
                    ## read meta from pickle
                    pclmeta = serializer.load( f )
                    if not pclmeta.needsUpdate( txtmeta ):
                        logger.debug( "we can use expres from pickle file %s" % pclfile )
                        expres = serializer.load( f )
                    else:
                        logger.debug( "we cannot use expres from pickle file %s" % pclfile )
                        logger.debug( "txt meta %s" % txtmeta )
                        logger.debug( "pcl meta %s" % pclmeta )
                        logger.debug( "pcl meta needs update %s" % pclmeta.needsUpdate( txtmeta ) )
        except IOError as e:
            logger.error( "exception %s" % e )
        if not expres: ## create from text file
            expres = ExpResult(root, discard_zeroes = self.txt_meta.discard_zeroes,
                databaseParticles = self.databaseParticles)
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

        :param analysisIDs: list of analysis ids ([CMS-SUS-13-006,...]). Can
                            be wildcarded with usual shell wildcards: * ? [<letters>]
                            Furthermore, the centre-of-mass energy can be chosen
                            as suffix, e.g. ":13*TeV". Note that the asterisk
                            in the suffix is not a wildcard.
        :param datasetIDs: list of dataset ids ([ANA-CUT0,...]). Can be wildcarded
                            with usual shell wildcards: * ? [<letters>]
        :param txnames: list of txnames ([TChiWZ,...]). Can be wildcarded with
                            usual shell wildcards: * ? [<letters>]
        :param dataTypes: dataType of the analysis (all, efficiencyMap or upperLimit)
                            Can be wildcarded with usual shell wildcards: * ? [<letters>]
        :param useSuperseded: If False, the supersededBy results will not be included
                              (deprecated)
        :param useNonValidated: If False, the results with validated = False
                                will not be included
        :param onlyWithExpected: Return only those results that have expected values
                 also. Note that this is trivially fulfilled for all efficiency maps.
        :returns: list of ExpResult objects or the ExpResult object if the list
                  contains only one result

        """
        if type(analysisIDs)==str: analysisIDs=[analysisIDs]
        if type(datasetIDs)==str: datasetIDs=[datasetIDs]
        if type(txnames)==str: txnames=[txnames]
        if type(dataTypes)==str: dataTypes=[dataTypes]

        import fnmatch
        expResultList = []
        for expResult in self.expResultList:
            superseded = None
            if hasattr(expResult.globalInfo,'supersededBy'):
                superseded = expResult.globalInfo.supersededBy.replace(" ","")
            if superseded and (not useSuperseded):
                continue

            analysisID = expResult.globalInfo.getInfo('id')
            sqrts = expResult.globalInfo.getInfo('sqrts')

            # Skip analysis not containing any of the required ids:
            if analysisIDs != ['all']:
                hits=False
                for patternString in analysisIDs:
                    # Extract centre-of-mass energy
                    # Assuming 0 or 1 colons.
                    pattern = patternString.split(':')
                    hits = fnmatch.filter( [ analysisID ], pattern[0] )
                    if len( pattern ) > 1:
                        # Parse suffix
                        # Accepted Strings: ":13", ":13*TeV", ":13TeV", ":13 TeV"
                        # Everything else will yield an error at the unum-conversion (eval())
                        if pattern[1].endswith('TeV'):
                            pattern[1] = pattern[1][:-3]
                        if pattern[1][-1] in [" ", "*"]:
                            pattern[1] = pattern[1][:-1]
                        pattern[1] += "*TeV"
                        if sqrts != eval(pattern[1]):
                            hits = False
                    if hits:
                        break
                        # continue
                if not hits:
                    continue

            newExpResult = ExpResult()
            newExpResult.path = expResult.path
            newExpResult.globalInfo = expResult.globalInfo
            newExpResult.datasets = []

            for dataset in expResult.datasets:
                if dataTypes != ['all']:
                    hits=False
                    for pattern in dataTypes:
                        hits = fnmatch.filter( [ dataset.dataInfo.dataType ], pattern )
                        if hits:
                            break
                            #continue
                    if not hits:
                        continue

                if hasattr(dataset.dataInfo, 'dataId') and datasetIDs != ['all']:
                    hits=False
                    if datasetIDs == None:
                        datasetIDs = [ None ]
                    for pattern in datasetIDs:
                        hits = fnmatch.filter( [ str(dataset.dataInfo.dataId) ], str(pattern) )
                        if hits:
                            break
                            # continue
                    if not hits:
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
                        #Replaced by wildcard-evaluation below (2018-04-06 mat)
                        hits=False
                        for pattern in txnames:
                            hits = fnmatch.filter( [ txname.txName ], pattern )
                            if hits: # one match is enough
                                break
                        if not hits:
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

    def updateBinaryFile( self ):
        """ write a binar db file, but only if
            necessary. """
        if self.needsUpdate():
            logger.debug( "Binary db file needs an update." )
            self.createBinaryFile()
        else:
            logger.debug( "Binary db file does not need an update." )

class ExpResultList(object):
    """
    Holds a list of ExpResult objects for printout.
    """

    def __init__(self, expResList):
        """
        :param expResultList: list of ExpResult objects
        """

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
        db = Database( args.database, force_load="txt" )
        db.createBinaryFile()
        sys.exit()
    db = Database( args.database )
    if args.update:
        db.updateBinaryFile()
    if args.check:
        db.checkBinaryFile()
    if args.time:
        t0=time.time()
        expResult = db.loadBinaryFile( lastm_only = False )
        t1=time.time()
        print( "Time it took reading binary db file: %.1f s." % (t1-t0) )
        txtdb = db.loadTextDatabase()
        t2=time.time()
        print( "Time it took reading text   file: %.1f s." % (t2-t1) )
    if args.read:
        db = db.loadBinaryFile( lastm_only = False )
        listOfExpRes = db.getExpResults()
        for expResult in listOfExpRes:
            print(expResult)
