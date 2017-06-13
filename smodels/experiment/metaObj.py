#!/usr/bin/env python

"""
.. module:: metaObj
   :synopsis: Contains the Meta object used to track the provenance of
              the pickle file.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os
import sys
from smodels.tools.smodelsLogging import logger

class Meta(object):
    current_version = 201 ## the current format version

    """ The Meta object holds all meta information regarding the
        database, like number of analyses, last time of modification, ...
        This info is needed to understand if we have to re-pickle. """

    def __init__ ( self, pathname, discard_zeroes=None, mtime=None, filecount=None,
                   hasFastLim=None, databaseVersion=None, format_version=current_version,
                   python=sys.version ):
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

    def versionFromFile ( self ):
        """
        Retrieves the version of the database using the version file.
        """
        if self.databaseVersion or self.isPickle():
            return self.databaseVersion
        try:
            vfile = os.path.join ( self.pathname, "version" )
            versionFile = open( vfile )
            content = versionFile.readlines()
            versionFile.close()
            line = content[0].strip()
            logger.debug("Found version file %s with content ``%s''" \
                   % ( vfile, line) )
            return line
        except IOError:
            # logger.error('There is no version file %s', vfile )
            return 'unknown version'

    def __str__ ( self ):
        ret  = "Meta: path =%s\n" % self.pathname
        ret += "      mtime=%.1f" % self.mtime
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

    def determineLastModified ( self ):
        """ compute the last modified timestamp, plus count
            number of files. Only if text db """
        if self.isPickle() or self.mtime:
            return
        lastm = 0
        count = 0
        versionfile = os.path.join ( self.pathname, "version" )
        if not os.path.exists ( versionfile ):
            logger.debug("%s does not exist." % versionfile )
            # sys.exit()
        else:
            lastm = os.stat(versionfile).st_mtime
            count=1
        topdir = os.listdir ( self.pathname )
        for File in topdir:
            subdir = os.path.join ( self.pathname, File )
            if not os.path.isdir ( subdir ) or File in [ ".git" ]:
                continue
            (lastm,tcount) = self.lastModifiedSubDir ( subdir, lastm )
            count+=tcount+1
        self.mtime=lastm
        self.filecount = count

    def lastModifiedSubDir ( self, subdir, lastm ):
        """
        Return the last modified timestamp of subdir (working recursively)
        plus the number of files.

        :param subdir: directory name that is checked
        :param lastm: the most recent timestamp so far, plus number of files
        :returns: the most recent timestamp, and the number of files
        """
        ret = lastm
        ctr=0
        for f in os.listdir ( subdir ):
            if f in [ "orig", "sms.root", "validation", ".git" ]:
                continue
            if f[-1:]=="~":
                continue
            if f[0]==".":
                continue
            if f[-3:]==".py":
                continue
            lf = os.path.join ( subdir, f )
            if os.path.isdir ( lf ):
                (ret,tctr) = self.lastModifiedSubDir ( lf, ret )
                ctr+=tctr+1
            else:
                ctr+=1
                tmp = os.stat ( lf ).st_mtime
                if tmp > ret:
                    ret = tmp
        return (ret,ctr)

    def sameAs ( self, other ):
        """ check if it is the same database version """
        if type(self) != type(other):
            return False
        return (self.databaseVersion == other.databaseVersion)

    def printFastlimBanner ( self ):
        """ check if fastlim appears in data.
            If yes, print a statement to stdout. """
        if not self.hasFastLim: return
        logger.info ( "FastLim v1.1 efficiencies loaded. Please cite: arXiv:1402.0492, EPJC74 (2014) 11" )

    def __eq__ ( self, other ):
        return self == other

    def needsUpdate ( self, current ):
        """ do we need an update, with respect to <current>.
            so <current> is the text database, <self> the pcl.
        """
        if self.mtime < current.mtime: ## someone tinkered
            return True
        if self.filecount != current.filecount:
            return True ## number of files changed
        if self.databaseVersion != current.databaseVersion:
            return True ## number of files changed
        if self.discard_zeroes != current.discard_zeroes:
            return True ## flag changed
        if self.format_version != current.format_version:
            return True ## pickle file format version changed
        if self.python != current.python:
            return True ## different python
        return False
