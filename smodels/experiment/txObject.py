#!/usr/bin/env python

"""
.. module:: experimentalTopology
   :synopsis: Holds the txName object retrieved from smodels-database 
   in order to produce summaryplots.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Michael Traub <michael.traub@gmx.at>

"""    

import logging
#import setPath
#import sys
from txNames import TxDecay

FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.ERROR)

class TxObject(TxDecay):
    """Contains all txName-specific information (e.g. experimentIDs and runs that 
    contain this txName, category, particles resp. production mode, ...)
    
    """
    def __init__ (self, txName, txDict): #give the info.txt-objects!
        self._name = txName
        self._txDict = txDict
        self._runs = [key for key in self._txDict]
        self._experimentIDs = self._expIDs
        self._verbosity = 'error'
        TxDecay.__init__(self,self._name)
     
    @property
    def verbosity(self):
        """Tells the level the logger is set to.
        
        """
        return self._verbosity
        
    @verbosity.setter
    def verbosity(self, level):
        """Set the logger to specified level.
        
        """
        self._verbosity = level
        self._setLogLevel(level)

    def __str__(self):
        ret = "%s [%s]: " % ( self.name, self.category )
        ret += ", ".join ( self.constraints ).replace("'","") 
        return ret
    
    @property    
    def name(self):
        return self._name
    
    @property
    def experimentIDs(self):
        return self._experimentIDs
    
    @property
    def runs(self):
        return self._runs
        
    @property
    def category(self):
        return self._category
    
    @property
    def constraints(self):
        return self._constraints
    
    @property
    def massParametrizations(self):
        return self._massParametrizations
    
    @property
    def extensions(self):
        return self._extensions
    
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
    def _expIDs(self):
        """Extracts all the experimentIDs given as inner keys of nested txDict.
        
        """
        expIDs = []
        for r in self._runs:
            for a in self._txDict[r]:
                expIDs.append(a)
        return expIDs
    
    def _getInfoProperty(self, info, requested):
        """Retrieves the requested property of the given info.txt object.
        
        """
        return getattr(info, requested)
    
    @property    
    def _category(self):
        """Takes the category for this txName from every info.txt, 
        compares them and returns the string if they are all the same. 
        Raises an ERROR and returns None if they are not!
        
        """
        
        cats = []
        for run in self._runs:
            for ana in self._expIDs:
                try:
                    category = self._getInfoProperty(self._txDict[run][ana], 'category')[self.name]
                    if cats.count(category) == 0:
                        cats.append(category)
                    if cats and cats.count(category) == 0:
                        logger.error('There are different categories for txName %s! \
                        Please check the database entry %s-%s!' %(self._name, run, ana))
                except KeyError:
                    logger.warning('The category for %s is missing! Please \
                    check the database entry %s-%s!' %(self._name, run, ana))
        logger.debug('List of categories: %s.' %cats)
        if len(cats) == 0:
            logger.error('Could not get any category information for %s.' % \
                    self._name )
            return None
        if len(cats) == 1:
            return cats[0]
        
        logger.error('Unable to get consistent category for txName %s: %s' % \
                      (self._name,cats) )
        return None
    
    @property    
    def _constraints(self):
        """Takes the constraints for this txName from every info.txt, 
        returns a list containing all available constraints.
        
        """
        
        const = []
        for run in self._runs:
            for ana in self._expIDs:
                try:
                    c = self._getInfoProperty(self._txDict[run][ana], 'constraints')[self.name]
                    if not c in const:
                        const.append(c)
                except KeyError:
                    logger.warning('The constraint for %s is missing! \
                    Please check the database entry %s-%s!' %(self._name, run, ana))
        logger.debug('List of constraints: %s.' %const)
        return const
        
    
    @property
    def _massParametrizations(self):
        """Retrieves all conditions for the third mass, available for this txName.
        
        """
        massConds = {}
        for run in self._runs:
            for ana in self._expIDs:
                try:
                    axes = self._getInfoProperty(self._txDict[run][ana], 'axes')[self.name]
                    axes = [ax for ax in axes if 'mz' in ax and ax['mz']]
                    for ax in axes:
                        if not ax['extension']:
                            ax['extension'] = self.name
                        else:
                            ax['extension'] = self.name + ax['extension']
                        if not ax['extension'] in massConds:
                            massConds[ax['extension']] = ax['mz']
                except KeyError:
                    logger.warning('The axes for %s are missing! \
                    Please check the database entry %s-%s!' %(self._name, run, ana))
        return massConds
    
    @property
    def _extensions(self):
        """Retrieves all extensions, available for this txName.
        
        """
        ext = []
        for run in self._runs:
            for ana in self._expIDs:
                try:
                    extensions = self._getInfoProperty(self._txDict[run][ana], 'extensions')[self.name]
                    for ex in extensions:
                        if not ex in ext:
                            ext.append(ex)
                except KeyError:
                    logger.warning('The extensions for %s are missing! \
                    Please check the database entry %s-%s!' %(self._name, run, ana))
        return ext
        
  
    #def getPrettyName       # particles resp. productionmode
    #def treatMasssplitting
    #def setAnalyses
    #def refreshAnalyses
