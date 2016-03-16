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
import pickle
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
        # logger.info ( "Last modified time stamp" )
        # logger.info ( "db at %s " % self.db_dir )
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
        return lastm,count

    def createPickleFile ( self ):
        """ create a pcl file from the text database,
            potentially overwriting an old pcl file. """
        logger.info ( "Creating pickle file:" )
        logger.info ( " * compute last modified timestamp." )
        lastm,count=self.lastModifiedAndFileCount()
        logger.info (  " * compute timestamp: %s filecount: %d" % \
                ( time.ctime ( lastm ), count ) )
        logger.info (  " * create %s" % self.pclfile )
        with open ( self.pclfile, "w" ) as f:
            pickle.dump ( lastm, f )
            pickle.dump ( count, f )
            logger.info (  " * load text database" )
            self.loadTextDatabase() 
            pickle.dump ( self.txt_db, f )

    def loadPickleFile ( self, lastm_only = False ):
        """ load a pickle file, returning
            last modified, file count, database.
        :param lastm_only: if true, the database itself is not read.
        """
        with open ( self.pclfile, "r" ) as f:
            lastm = pickle.load ( f )
            count = pickle.load ( f )
            db=None
            if not lastm_only:
                db = pickle.load ( f )
        return lastm,count,db
        
if __name__ == "__main__":
    import argparse
    """ Run as a script, this checks and/or writes database.pcl files """
    argparser = argparse.ArgumentParser(description='simple script to check \
            and/or write database.pcl files')
    argparser.add_argument('-c', '--check', help='check pickle file',
                           action='store_true')
    argparser.add_argument('-w', '--write', help='force writing pickle file',
                           action='store_true')
    argparser.add_argument('-u', '--update', help='update pickle file, if necessary',
                           action='store_true')
    argparser.add_argument('-d', '--database', help='directory name of database', default="../../../smodels-database/" )
    args = argparser.parse_args()
    logger.setLevel(level=logging.DEBUG )
    pickler = DbPickler ( args.database )
    print pickler
    if args.write:
        pickler.createPickleFile()
    
