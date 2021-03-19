#!/usr/bin/env python3

"""
.. module:: metaObj
   :synopsis: Contains the Meta object used to track the provenance of
              the pickle file.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os
import sys
import time
from smodels.tools.smodelsLogging import logger

class Meta(object):
    current_version = 214 ## the current format version

    """ The Meta object holds all meta information regarding the
        database, like number of analyses, last time of modification, ...
        This info is needed to understand if we have to re-pickle. """

    def __init__(self, pathname, discard_zeroes=None, mtime=None, filecount=None,
                   hasFastLim=None, databaseVersion=None, format_version=current_version,
                   python=sys.version):
        """
        :param pathname: filename of pickle file, or dirname of text files
        :param discard_zeroes: do we discard zeroes?
        :param mtime: last modification time stamps
        :param filecount: number of files
        :param hasFastLim: fastlim in the database?
        :param databaseVersion: version of database
        :param format_version: format version of pickle file
        :param python: python version
        """
        self.pathname = pathname
        self.discard_zeroes = discard_zeroes
        self.mtime = mtime
        self.filecount = filecount
        self.hasFastLim = hasFastLim
        self.format_version = format_version
        self.python = python
        self.databaseVersion = databaseVersion
        self.versionFromFile()
        self.determineLastModified()

    def getPickleFileName ( self ):
        """ get canonical pickle file name """
        hfl=""
        if self.hasFastLim:
            hfl="1"
        if self.hasFastLim==False:
            hfl="0"
        return "db%s%d.pcl" % ( self.python[0], self.discard_zeroes )
        # return "db%s%d%s.pcl" % ( self.python[0], self.discard_zeroes, hfl )

    def versionFromFile ( self ):
        """
        Retrieves the version of the database using the version file.
        """
        # logger.error ( "versionFromFile, %s, %s, %s" % ( self.pathname, self.databaseVersion, self.isPickle() ) )
        if self.databaseVersion or self.isPickle():
            return self.databaseVersion
        try:
            vfile = os.path.join ( self.pathname, "version" )
            versionFile = open( vfile )
            content = versionFile.readlines()
            versionFile.close()
            line = content[0].strip()
            #logger.debug("Found version file %s with content ``%s''" \
            #       % ( vfile, line) )
            self.databaseVersion = line
            # return line
        except IOError:
            pass
            # logger.error('There is no version file %s', vfile )
            # return 'unknown version'

    def __str__ ( self ):
        ret  = "Meta: path =%s\n" % self.pathname
        ret += "      mtime=%s" % time.ctime ( self.mtime )
        ret += ", filecount=%d" % self.filecount
        ret += ", discard_0=%d" % self.discard_zeroes
        ret += ", fl=%s" % self.hasFastLim
        ret += ", format_version=%d" % self.format_version
        ret += ", dbVersion=%s" % self.databaseVersion
        return ret

    def isPickle ( self ):
        """ is this meta info from a pickle file? """
        if not os.path.isdir ( self.pathname ):
            return True

    def cTime ( self ):
        return time.ctime ( self.mtime )

    def determineLastModified ( self, force=False ):
        """ compute the last modified timestamp, plus count
            number of files. Only if text db """
        if self.isPickle():
            return self.mtime, self.filecount
        if force==False and self.mtime:
            return self.mtime, self.filecount
        self.mtime = 0
        self.filecount = 0
        versionfile = os.path.join ( self.pathname, "version" )
        if os.path.exists ( versionfile ):
            self.mtime = os.stat(versionfile).st_mtime
            self.filecount=1
        topdir = os.listdir ( self.pathname )
        self.lastModifiedSubDir ( self.pathname )
        # (self.mtime,self.filecount) = self.lastModifiedSubDir ( self.pathname, lastm )
        # return self.mtime, self.filecount

    def lastModifiedSubDir ( self, subdir ):
        """
        Return the last modified timestamp of subdir (working recursively)
        plus the number of files.

        :param subdir: directory name that is checked
        :param lastm: the most recent timestamp so far, plus number of files
        :returns: the most recent timestamp, and the number of files
        """
        #ret = lastm
        #ctr=0
        for f in os.listdir ( subdir ):
            if f in [ "orig", "sms.root", "validation", ".git" ]:
                continue
            if f[-1:]=="~":
                continue
            if f[0]==".":
                continue
            if f[-3:]==".py":
                continue
            if f[-4:]==".pcl":
                continue
            lf = os.path.join ( subdir, f )
            if os.path.isdir ( lf ):
                self.lastModifiedSubDir ( lf )
                # (ret,tctr) = self.lastModifiedSubDir ( lf, ret )
                self.filecount+=1
                # ctr+=tctr+1
            else:
                self.filecount+=1
                #ctr+=1
                tmp = os.stat ( lf ).st_mtime
                if tmp > self.mtime:
                    self.mtime = tmp
                #    ret = tmp
                
        #return (ret,ctr)

    def sameAs ( self, other ):
        """ check if it is the same database version """
        if type(self) != type(other):
            return False
        return (self.databaseVersion == other.databaseVersion)

    def printFastlimBanner ( self ):
        """ check if fastlim appears in data.
            If yes, print a statement to stdout. """
        if not self.hasFastLim: return
        # print ( "FastLim v1.1 efficiencies loaded. Please cite: arXiv:1402.0492, EPJC74 (2014) 11" )
        logger.info ( "FastLim v1.1 efficiencies loaded. Please cite: arXiv:1402.0492, EPJC74 (2014) 11" )

    def __eq__ ( self, other ):
        if other == None: return False
        if self.pathname != other.pathname:
            return False
        if self.discard_zeroes != other.discard_zeroes:
            return False
        if self.mtime != other.mtime:
            return False
        if self.filecount != other.filecount:
            return False
        if self.hasFastLim != other.hasFastLim:
            return False
        if self.format_version != other.format_version:
            return False
        if self.python != other.python:
            return False
        if self.databaseVersion != other.databaseVersion:
            return False
        return True

    def needsUpdate ( self, current ):
        """ do we need an update, with respect to <current>.
            so <current> is the text database, <self> the pcl.
        """
        if self.mtime < current.mtime: ## someone tinkered
            return True
        if self.filecount != current.filecount:
            return True ## number of files changed
        #if self.databaseVersion != current.databaseVersion:
        #    return True ## database version changed
        if self.discard_zeroes != current.discard_zeroes:
            return True ## flag changed
        if self.format_version != current.format_version:
            return True ## pickle file format version changed
        if self.python != current.python:
            return True ## different python
        return False
