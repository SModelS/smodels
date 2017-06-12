#!/usr/bin/env python

"""
.. module:: metaObj
   :synopsis: Contains the Meta object used to track the provenance of 
              the pickle file.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys

class Meta(object):
    current_version = 200 ## the current format version

    """ The Meta object holds all meta information regarding the
        database, like number of analyses, last time of modification, ...
        Needed to understand if we have to re-pickle. """

    def __init__ ( self, pathname=None, mtime=None, filecount=None, 
                   hasFastLim=None, discard_zeroes=None,
                   format_version=current_version, 
                   python=sys.version, databaseVersion=None ):
        """
        :param pathname: filename of pickle file, or dirname of text files
        :param mtime: last modification time stamps
        :param filecount: number of files
        :param hasFastLim: fastlim in the database?
        :param discard_zeroes: do we discard zeroes?
        :param format_version: format version of pickle file
        :param python: python version
        :param databaseVersion: version of database
        """
        self.pathname = pathname
        self.mtime = mtime
        self.filecount = filecount
        self.hasFastLim = hasFastLim
        self.discard_zeroes = discard_zeroes
        self.format_version = format_version
        self.python = python
        self.databaseVersion = databaseVersion

    def sameAs ( self, other ):
        """ check if it is the same database version """
        return (self.databaseVersion == other.databaseVersion)

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
