#!/usr/bin/env python

"""
.. module:: tools.dbPickler
   :synopsis: Contains the facility to pickle the database (i.e.
              store it in an efficient binary format), and
              determine if the database needs to be re-parsed
              by looking at the last-changed timestamps of all
              files

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import time
import sys
import os
# import pickle
import cPickle as pickle
from smodels.experiment.exceptions import DatabaseNotFoundException
from smodels.experiment.databaseObjects import Database
import logging
logger = logging.getLogger(__name__)

class DbPickler(object):
    """
        a class that encapsulates creating and reading database pickle files.
    """

    def __init__ ( self, db_dir ):
        """
        :param database_dir: Directory containing the database 
        """
        self.db_dir = db_dir
        self.txt_db = None
        self.pclfile = os.path.join ( self.db_dir, "database.pcl" )
        self.pcl_mtime = None, None
        self.pcl_db = None
        self.txt_mtime = None, None

    def __str__ ( self ):
        ret = "Database at: %s\n" % ( self.db_dir )
        return  ret
         
    def loadTextDatabase ( self ):
        """ loads the database from the text files """
        if self.txt_db: 
            # already loaded!
            return
        try:
            self.txt_db = Database( self.db_dir )
        except DatabaseNotFoundException:
            logger.error("Database not found in %s" % os.path.realpath(self.db_dir))
            sys.exit()
        return self.txt_db

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
        versionfile = os.path.join ( self.db_dir, "version" ) 
        if not os.path.exists ( versionfile ):
            logger.error("%s does not exist." % versionfile )
            sys.exit()
        lastm = os.stat(versionfile).st_mtime
        count=1
        topdir = os.listdir ( self.db_dir )
        for File in topdir:
            subdir = os.path.join ( self.db_dir, File )
            if not os.path.isdir ( subdir ) or File in [ ".git" ]:
                continue
            print subdir
            (lastm,tcount) = self.lastModifiedDir ( subdir, lastm )
            count+=tcount+1
        self.txt_mtime = lastm, count

    def createPickleFile ( self ):
        """ create a pcl file from the text database,
            potentially overwriting an old pcl file. """
        logger.debug ( "Creating pickle file:" )
        logger.debug ( " * compute last modified timestamp." )
        self.lastModifiedAndFileCount()
        logger.debug (  " * compute timestamp: %s filecount: %d" % \
                ( time.ctime ( self.txt_mtime[0] ), self.txt_mtime[1] ) )
        logger.debug (  " * create %s" % self.pclfile )
        with open ( self.pclfile, "w" ) as f:
            pickle.dump ( self.txt_mtime, f )
            logger.debug (  " * load text database" )
            self.loadTextDatabase() 
            logger.debug (  " * write %s" % self.pclfile )
            pickle.dump ( self.txt_db, f, protocol=2 )
            logger.debug (  " * done writing %s" % self.pclfile )


    def needsUpdate ( self ):
        """ does the pickle file need an update? """
        self.lastModifiedAndFileCount()
        self.loadPickleFile ( lastm_only = True )
        return ( self.txt_mtime[0] > self.pcl_mtime[0] or \
                 self.txt_mtime[1] != self.pcl_mtime[1] )

    def updatePickleFile ( self ):
        """ write a pickle file, but only if 
            necessary. """
        if self.needsUpdate():
            logger.debug ( "pickle file needs an update." )
            self.createPickleFile()
        else:
            logger.debug ( "pickle file does not need an update." )

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
            self.pcl_mtime = pickle.load ( f )
            if not lastm_only:
                self.pcl_db = pickle.load ( f )
        return self.pcl_db

    def checkPickleFile ( self ):
        nu=self.needsUpdate()
        logger.debug ( "Checking pickle file." )
        logger.debug ( "Pickle file dates to %s(%d)" % \
                      ( time.ctime(self.pcl_mtime[0]),self.pcl_mtime[1] ) )
        logger.debug ( "Database dates to %s(%d)" % \
                      ( time.ctime(self.txt_mtime[0]),self.txt_mtime[1] ) )
        if nu:
            logger.debug ( "pickle file needs an update." )
        else:
            logger.debug ( "pickle file does not need an update." )
        return nu

            
        
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
    argparser.add_argument('-d', '--database', help='directory name of database', 
                            default="../../../smodels-database/" )
    args = argparser.parse_args()
    logger.setLevel(level=logging.DEBUG )
    pickler = DbPickler ( args.database )
    logger.debug ( "%s" % pickler )
    if args.write:
        pickler.createPickleFile()
    if args.update:
        pickler.updatePickleFile()
    if args.check:
        pickler.checkPickleFile()
    if args.time:
        t0=time.time()
        db = pickler.loadPickleFile ( lastm_only = False )
        t1=time.time()
        print "Time it took reading pickle file: %.1f s." % (t1-t0)
        txtdb = pickler.loadTextDatabase()
        t2=time.time()
        print "Time it took reading text   file: %.1f s." % (t2-t1)
    if args.read:
        db = pickler.loadPickleFile ( lastm_only = False )
        listOfExpRes = db.getExpResults() 
        for expResult in listOfExpRes:
            print expResult
        #txtdb=pickler.loadTextDatabase()
        #listOfExpRes = txtdb.getExpResults() 
        #for expResult in listOfExpRes:
        #    print expResult
