#!/usr/bin/env python

"""
.. module:: ExperimentIDObject
   :synopsis: Holds the experimentID object retrieved from smodels-database.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Michael Traub <michael.traub@gmx.at>

"""    

import logging
#import setPath
#import sys

FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.ERROR)
    
class ExperimentIDObject(object):
    
    """Contains all experimentID-specific information 
    (e.g. PAS, lumi, publication-url, ...)
    
    """
        
    def __init__(self, experimentID, infotxt, run, smsroot, smspy):
        self._name = experimentID
        self._infotxt = infotxt
        self._run = run
        self._smsroot = smsroot
        self._smspy = smspy
        self._verbosity = 'error'
        self._dummies = ['not yet assigned', 'not yet assigned', '', 'not yet assignet']

    def __str__(self ):
        ret = "%s [%s, %s] <<%s>>" % \
             (self.name, self.experiment, self.run, ", ".join(self.txNames))
        return ret
     
    @property
    def verbosity(self):
        """Tells the level the logger is set to.
        
        """
        return self._verbosity
        
    @verbosity.setter
    def verbosity(self, level):
        """Set the logger to specified level.
        
        """
        level = self._validateLevel(level)
        self._verbosity = level
        self._setLogLevel(level)
    
    @property
    def lumi(self):
        return self._lumi

    @property
    def publishedData(self):
        return eval(str(self._parseMetaInfo('publisheddata')))
    
    @property
    def extendedTxNames(self):
        return self._extendedTxNames
    
    @property
    def pas(self):
        return self._parseMetaInfo('pas')
     
    @property
    def ID(self):
        return self._parseMetaInfo('id') 
     
    @property    
    def url(self):
        if self._parseMetaInfo('url'):
            return self._parseMetaInfo('url')
        return None
    
    @property    
    def hasUrl(self):
        if self._parseMetaInfo('url') and not \
        self._parseMetaInfo('url') in self._dummies:
            return True
        return False
    
    @property    
    def experiment(self):
        if 'ATLAS' in self.ID:
            return 'ATLAS'
        if 'CMS' in self.ID:
            return 'CMS'
        else:
            return 'unknown experiment'
        
    @property    
    def comment(self):
        return self._parseMetaInfo('comment')
    
    @property
    def prettyName(self):
        return self._parseMetaInfo('prettyname')
        
    @property    
    def hasConstraints(self):
        """Checks if there are any constraints for this experimentID.
        
        """
        
        if self._parseInfo('constraint') and not \
        self._parseInfo('constraint') in self._dummies:
            return True
        return False
        
    @property
    def constraints(self):
        """Retrieves all the constraints stored in the info.txt file.
        
        """
        if self._parseInfo('constraint'):
            return self._parseInfo('constraint')
        return None
     
    @property
    def hasConditions(self):
        """Checks is there are any conditions for this experimentID.
        
        """
        if self._parseInfo('condition'): return True
        return False
        
        
    @property
    def conditions(self):
        """Retrieves all the conditions stored in the info.txt file.
        
        """
        return self._parseInfo('condition')
    
    @property
    def hasFuzzyConditions(self):
        """Checks is there are any conditions for this experimentID.
        
        """
        if self._parseInfo('fuzzycondition'): return True
        return False
        
        
    @property
    def fuzzyConditions(self):
        """Retrieves all the conditions stored in the info.txt file.
        
        """
        return self._parseInfo('fuzzycondition')
    
    @property    
    def private(self):
        """States if the experimentID is private (True) or public (False).

        """
        t = self._parseMetaInfo('private')
        if not t: 
            return False
        t = t.lower()
        return t in [ "1", "yes", "true" ]
    
    @property
    def public(self):
        """States if the experimentID is public (True) i.e. NOT private 
        or private (False) i.e. private.
        
        """
        if self.private:
            return False
        return True
    
    @property    
    def hasArxiv(self):
        if self._parseMetaInfo('arxiv') and not \
        self._parseMetaInfo('arxiv') in self._dummies : return True
        return False
        
    @property        
    def arxiv(self):
        if self.hasArxiv:
            return self._parseMetaInfo('arxiv')
        return None
    
    @property    
    def hasPublication(self):
        if self._parseMetaInfo('publication') and not \
        self._parseMetaInfo('publication') in self._dummies: return True
        return False

    @property    
    def publication(self):
        if self.hasPublication:
            self._parseMetaInfo('publication')
        return None
    
    @property
    def hasAxes(self):
        if self._parseMetaInfo('axes'): return True
        return False
        
    @property
    def axes(self):
        """Retrieves the axes information stored in the axes-labeled line of 
        info.txt.
        :return: {'txName': [{'mx': 'x-axis', 'my': 'y-axis', 
        'mz': ('mass condition', value)}]}
       
    
        """
        if self.hasAxes:
            return self._axes
        return None
    
    @property    
    def isChecked(self):
        if self._parseMetaInfo('checked') and not \
        self._parseMetaInfo('checked') in self._dummies: return True
        return False
        
    @property        
    def checked(self):
        return self._parseMetaInfo('checked')
        
    @property        
    def supersedes(self):
        return self._parseMetaInfo('supersedes')
        
    @property        
    def superseded(self):
        return self._parseMetaInfo('superseded_by')
    
    @property
    def isPublished(self):
        if self._parseMetaInfo('arxiv') or self._parseMetaInfo('publication'):
            return True
        return False
    
    @property    
    def name(self):
        return self._name
    
    @property    
    def sqrts(self):
        return self._sqrts
    
    @property    
    def run(self):
        return self._run
        
    @property
    def txNames(self):
        """Retrieves all the txNames this experimentID has results for as strings.
        
        """
        return self._getInfoProperty('txNames')
        
    @property    
    def extendedTxNames(self):
        """Retrieves all the txNames with their particular extensions 
        (referring to possible mass conditions) this experimentID has results 
        for as strings.
        
        """
        return self._extendedTxNames
    
    @property    
    def massParametrizations(self):
        """Retrieves a dictionary giving the correlation between all the 
        txNames, their particular extensions and their mass parametrization.
        :return: {txName: [{extended txName: (parametrization, value)}]}
        
        """
        return self._massParametrizations
    
    
    @property
    def exclusions(self):
        """Retrieves all the exclusion values stored in the info.txt file.
        
        """
        return self._getInfoProperty('exclusions')
    
    @property    
    def hasROOT(self):
        if self._smsroot: return True
        return False
    
    @property    
    def hasPY(self):
        if self._smspy: return True
        return False
    
    @property    
    def ROOT(self):
        return self._smsroot
        
    @property    
    def PY(self):
        return self._smspy
    
    def _validateLevel(self, level):
        """Validates given level for pythons logger module.
        
        """
        if not level.lower() in ['debug', 'info', 'warning', 'error']:
            logger.error('No valid level for verbosity: %s! Browser will \
            use default setting!' %level)
            return 'error'
        return level.lower()
            
    def _setLogLevel(self, level = 'error'):
        if level == 'debug':
            logger.setLevel(level=logging.DEBUG)
        if level == 'info':
            logger.setLevel(level=logging.INFO)
        if level == 'warning':
            logger.setLevel(level=logging.WARNING)
        if level == 'error':
            pass
    
    def _getInfoProperty(self, requested):
        """Retrieves the requested property of the given info.txt object.
        
        """
        return getattr(self._infotxt, requested)
    
    def _parseMetaInfo(self, requested):
        metaInf = self._getInfoProperty('metaInfo')
        if not requested in metaInf:
            logger.warning('Requested keyword %s could not be found for %s!' \
            %(requested, self._name))
            return None
        return metaInf[requested]
        
    def _parseInfo(self, requested):
        inf = self._getInfoProperty('info')
        content = [line for line in inf if requested == line.split(':')[0].strip()]
        if not content:
            logger.warning('Requested lines %s could not be found for %s!' \
            %(requested, self._name))
            return None
        content = [line.split(':')[1].strip() for line in content]
        return content
    
    @property
    def _axes(self):
        """
        Creates the axes dictionary from the infotxt object's axes method.
        
        """
        axesDict = {}
        axes = self._getInfoProperty('axes')
        for tx in axes:
            axesDict[tx] = []
            for a in axes[tx]:
                axe = {key: a[key] for key in a if key != 'extension'} 
                axesDict[tx].append(axe)
        return axesDict        
    
    @property
    def _extendedTxNames(self):
        """Creates all the extended txNames by concatenating the name of
        the txName and possible extensions and feeds them into a dictionary.
        
        """
        extTopos = {}
        extensions = self._getInfoProperty('extensions')
        if not extensions:
            logger.warning('Could not get extended txNames for %s' %self.name)
            return None
        for tx in self.txNames:
            extTopos[tx] = []
            if not extensions[tx]:
                extTopos[tx].append(tx)
                continue
            for ext in extensions[tx]:
                extTopos[tx].append(tx + ext)
        return extTopos        
    
    @property
    def _massParametrizations(self):
        paramDict = {}
        axes = self._getInfoProperty('axes')
        for tx in axes:
            for a in axes[tx]:
                if a['extension']:
                    ext = a['extension']
                else:
                    ext = ''
                if 'mz' in a:
                    paramDict[tx + ext] = a['mz']
                else:
                    paramDict[tx + ext] = 'more than one'
        for tx in self.txNames:
            if not tx in paramDict:
                logger.error('Mass parametrization for txName %s is missing!' \
                %tx)
        return paramDict
        
    @property
    def _sqrts(self):
        s = self._parseMetaInfo('sqrts')
        try:
            return float(s)
        except ValueError:
            try:
                return float(s.split('*')[0])
            except TypeError:
                if '8' in s: return 8.0
                if '7' in s: return 7.0
                if not s: return None
                
    @property
    def _lumi(self):
        l = self._parseMetaInfo('lumi')
        try:
            return float(l)
        except ValueError:
            try:
                return float(l.split('/')[0])
            except TypeError:
                return None
        
    
        
    #def getRestOfInfo => contact, arxiv, publisheddata ### check something missing?
