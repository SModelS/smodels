 #!/usr/bin/env python

"""
.. module:: analysisObject
   :synopsis: Holds the analysis object retrieved from smodels-database.
   Can be used for smodels only.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Michael Traub <michael.traub@gmx.at>

"""    

import logging, os
#import setPath
#import sys
from smodels.tools.physicsUnits import GeV
from unum import Unum
from databaseBrowserException import MassParametrizationException

FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.ERROR)


# ### FIX ME: When there is no intermediate particle in a given txName and one provides "condition" and "value" currently it throws an error in line 782: 
    #if self._txObject.name + ax['extension'] != self._txName:
    #TypeError: cannot concatenate 'str' and 'NoneType' objects
# this should be turned into a more human readable exception (e.g. MassParametrizationError), a warning should be given and the "extensionless" txName should be used to get e.g. upper limits or exclusion lines!


class AnalysisObject (object):
    """Contains all plot-specific information and objects (e.g. 
    exclusion lines, histograms, ...). Encapsules the plot objects to 
    handle different mass assumptions for given txName and analysis.
    
    """
    
    def __init__ (self, run, experimentIDObject, txObject, smspy):
        """Sets all private variables, especially self._plots 
        as list containing all available plots as objects.
    
        """
        self._txObject = txObject
        self._experimentIDObject = experimentIDObject
        self._smspy = smspy
        self._txName = txObject.name
        self._experimentID = experimentIDObject.name
        self._run = run
        logger.info('Creating analysis object for %s-%s-%s!' \
        %(self._run, self._experimentID, self._txName))
        self._extendedTxNames = self._experimentIDObject.extendedTxNames[self._txName]
        self._plots = self._getPlots
        self._verbosity = 'error'
    
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
    
    def __str__(self):
        ret = "ExperimentID: %s \nTxName: %s" %(str(self._experimentIDObject).split('<')[0], self._txObject)
        return ret
        
    @property
    def name(self):
        return self._experimentID + '-' + self._txName
        
        
    @property    
    def experimentIDObject(self):
        """Returns the analysis object linked to this analysis object.
        
        """
        return self._experimentIDObject
    
    @property    
    def txObject(self):
        """Returns the txName object linked to this analysis object.
        
        """
        return self._txObject
        
    @property    
    def isChecked(self):
        """Is this analysis object checked?
        
        """
        if self.checked: return True
        return False
    
    @property    
    def checked(self):
        """Retrieves checked_by entry from info.txt.
        
        """
        return self._checked
    
    @property    
    def txNameSet(self):
        """Returns all the extended txNames linked to this analysis.
        
        """
        return self._extendedTxNames
        
    @property
    def plots(self):
        """Returns a dictionary containing all available plot objects.
        
        """
        return self._plotDict    
    
    @property
    def plotNames(self):
        """Returns a list containing all available plot object names.
        
        """
        return self._plotNames    
    
    
    @property    
    def condition(self):
        """Retrieves the condition for this analysis object.
        
        """
        return self._condition
    
    @property    
    def fuzzyCondition(self):
        """Retrieves the fuzzy condition for this analysis object.
        
        """
        return self._fuzzyCondition
     
    @property    
    def constraint(self):
        """Retrieves constraint for this analysis object.
        
        """
        return self._constraint
    
    @property
    def axes(self):
        """Retrieves the axes for this analysis object.
        
        """
        return self._axes
        
    @property
    def members(self):
        """Retrieves the members of this analysis object.
        :return: {'extended txName': ('condition', value)}
        """
        return self._members
        
    def hasUpperLimitDicts(self, expected = False):
        """Checks which observed/expected upper limit dictionaries there are  
        for this analysis.
        
        """
        if self.upperLimitDicts(expected):
            return [key for key in self.upperLimitDicts(expected)]
        return None
        
    def upperLimitDicts(self, expected = False):
        """Retrieves all the observed/expected upper limit dictionaries 
        available for this analysis.
        # ### FIX ME: yields list -> for every extTxName => compare to exclusions to fix!
        """ 
        ulDicts = {}
        for plot in self._plots:
            ulDicts[plot.name] = plot.upperLimitDict(expected)
        return ulDicts
        
    def upperLimitDict(self, expected = False, condition = None, value = None):
        """Retrieves one observed/expected upper limit dictionary (out of all 
        upper limit dictionaries available for this txName). 
        If condition and value are None, the default mass assumptions will be used.
        Condition and value as a tuple specify the plot (out of this set) to be taken, e.g. ('fixedLSP', 50), ('massSplitting', 0.25), ...
        :param condition: condition for the third mass 
        :param value: value of the condition as either float or integer
        """
        
        extTxName = self._getExtendedTxName(condition = condition, value = value)
        plotName = self.name.replace(self._txName, extTxName)
        if not plotName in self.upperLimitDicts(expected = expected):
            if expected:
                logger.error('No expected upper limit dictionary could be found' + \
                             'for %s.' %extTxName)
                return None
            logger.error('No upper limit dictionary could be found for %s.' \
            %extTxName)
            return None
        return self.upperLimitDicts(expected = expected)[plotName]
    
        
    def _validateLevel(self, level):
        """Validates given level for pythons logger module.
        
        """
        if not level.lower() in ['debug', 'info', 'warning', 'error']:
            logger.error('No valid level for verbosity: %s.' + \
                         'Browser will use default setting.' %level)
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
        
    @property    
    def _checked(self):
        """Retrieves checked_by entry from info.txt.
        
        """
        infoLine = self._experimentIDObject.checked
        logger.debug('Got infoLine from ExperimentID-object: %s.' %infoLine)
        if not infoLine: return None
        # ### FIX ME: the IF below will be obsolet when 
        #the checked flag is fixed in every info.txt
        if len(infoLine.split()) == 1: 
            logger.warning('There is no information about single txNames.')
            return infoLine[0]
        
        infoLine = infoLine.split(',')
        logger.debug('First preprocessed infoLine: %s.' %infoLine)
        infoLine = [c.split(':')[1].strip() for c in infoLine if \
        c.split(':')[0].strip() == self._txName]
        logger.debug('Second preprocessed infoLine: %s.' %infoLine)
        if not infoLine:
            logger.warning('This Result is not checked.')
            return None
        logger.debug('Return value of infoLine: %s.' %infoLine[0])
        return infoLine[0].strip()
    
    @property
    def _condition(self):
        """Retrieves the condition for this plot.
        
        """
        cond = []
        if not self._experimentIDObject.hasConditions:
            logger.warning('No conditions available for analysis %s.' \
            %self._experimentID)
            return None
        cond = [c.split('->')[1].strip() for c in self._experimentIDObject.conditions \
        if c.split('->')[0].strip() == self._txName]
        if not cond:
            logger.warning('No condition available for plot %s.' \
            %self.name)
            return cond
        return cond[0]
    
    @property
    def _fuzzyCondition(self):
        """Retrieves the fuzzycondition for this plot.
        
        """
        fuzcond = []
        if not self._experimentIDObject.hasFuzzyConditions:
            logger.warning('No fuzzy conditions available for analysis %s.' \
            %self._experimentID)
            return None
        fuzcond = [c.split('->')[1].strip() for c in self._experimentIDObject.conditions \
        if c.split('->')[0].strip() == self._txName]
        if not fuzcond:
            logger.warning('No fuzzy condition available for plot %s.' \
            %self.name)
            return fuzcond
        return fuzcond[0]
    
    @property
    def _constraint(self):
        """Retrieves the constraint for this plot.
        
        """
        cons = []
        if not self._experimentIDObject.hasConstraints:
            logger.warning('No constraints available for analysis %s.' \
            %self._experimentID)
            return None
        cons = [c.split('->')[1].strip() for c in self._experimentIDObject.constraints \
        if c.split('->')[0].strip() == self._txName]
        if not cons:
            logger.warning('No constraints available for plot %s.' \
            %self.name)
            return cons
        return cons[0]
    
    

    @property
    def _getPlots(self):
        """Retrieves a list of all extended plots we have for this 
        analysis txName pair.
        
        """
        res = [Plot(extop, self._experimentIDObject, self._txObject, \
        self._smspy) for extop in self._extendedTxNames]
        return res
        
    @property
    def _plotDict(self):
        return {p.name: p for p in self._plots}
    
    @property
    def _plotNames(self):
        return [p.name for p in self._plots]
    
    
    def _getExtendedTxName(self, condition = None, value = None):
        """Creates the name of the extended txName (e.g. 'T6ttWWLSP050')
        :param condition: condition for the third mass as string (e.g. 'massSplitting')
        :param value: value for the condition as float or as integer if has dimension of
        a mass (e.g. 0.25, 2.0, 150)
        :return: 'extended txName'
        
        """
        #print '*****************', condition, value
        if not isinstance(value, Unum):
            if not condition or not value:
                return self._getDefaultExtendedTxName
       
        if isinstance(value, int):
            #value = addunit(value, 'GeV')
            #print '*****************', value
            value = value * GeV
                
        for plot in self._plots:
            if plot.axes['mz'] == (condition, value):
                return plot._txName
            else: continue 
        logger.warning('Unknown condition for third mass %s = %s.' + \
                       'Returning default.' %(condition, value))
        return self._getDefaultExtendedTxName
    
    @property
    def _getDefaultExtendedTxName(self):
        """Retrieves the default txName settings for this analysis object.
        :return: 'extended txName'
    
        """
        first = self._experimentIDObject._infotxt.axes[self._txName][0]
        if first['extension']:
            return self._txName + first['extension']
        else:
            return self._txName
        
    
    @property
    def _axes(self):
        """Retrieves the axes information for this plot.
        :return: {extended txName: 
        {'mx': mass on x-axis, 'my': mass on y-axes, 
        'mz': condition for intermediate mass}}
        
        """
        if not self._experimentIDObject.hasAxes:
            logger.warning('No axes information available for experimentID %s.' \
            %self._experimentID)
            return None
        try:
            return self._experimentIDObject.axes[self._txName]
        except KeyError:
            logger.warning('No axes information available for analysis' + \
                           'object %s.' %self.name)
            return None
    
    @property
    def _members(self):
        """Retrieves (condition, value) tuples for all the plots in this set.
        # ### FIX ME: if there is no information about mz this gives {'Tx': None}
        this is not very nice?
        """
        axes = self._experimentIDObject._infotxt.axes[self._txName]
        mems ={}
        for ax in axes:
            if ax['extension']:
                mems[self._txName + ax['extension']] = ax['mz']
            else:
                mems[self._txName] = ax['mz']
        return mems  


class Plot(object):
    """Contains all specific informations linked to one plot,
    where a plot denotes one mass plane for a pair of experimentID and 
    txName with a specific parametrization for the third mass (e.g. x-value = 050, 
    mass of LSP = 50 GeV, ...).
    
    """
    
    def __init__(self, txName, experimentIDObject, txObject, smspy):
        """Sets all private variables and retrieves the upper limit 
        dictionary.
        
        """
        self._txName = txName
        self._experimentIDObject = experimentIDObject
        self._txObject = txObject
        self._experimentID = experimentIDObject.name
        self._run = experimentIDObject.run
        self._smspy = smspy
    
    def __str__(self):
        ret = "%s" %self.name
        return ret
        
    @property
    def name(self):
        """Returns the name of this experimental plot as concatenated string.
        
        """
        return self._experimentID + '-' + self._txName
        
    @property
    def experimentIDObject(self):
        """Retrieves the experimental analysis object.
        """
        return self._experimentIDObject
    
    @property
    def txObject(self):
        """Retrieves the experimental txName object.
        """
        return self._txObject
    
    @property
    def siblings(self):
        """Retrieves the names of all the related plots.
        
        """
        return self._siblings
    
    @property
    def axes(self):
        """Retrieves the x- and y- axis of the upper limit histogram
        end the additional condition for the third mass, if there is any.
        """
        return self._axes

    
    def upperLimitDict(self, expected = False):
        """Retrieves the observed/expected cross section upper limit dictionary for this 
        plot from the sms.py file located in the database.
        
        """
        return self._upperLimitDict(expected)
    
    @property
    def _siblings(self):
        sibs = []
        for t in self._experimentIDObject.extendedTxNames[self._txObject.name]:
            sibs.append(self._experimentID + '-' + t)
        return sibs
            
    def _upperLimitDict(self, expected):
        """Retrieves the observed/expected cross section upper limit dictionary for this 
        plot from the sms.py file located in the database.
        
        """
        localSms = {}
        fakeDicts = [None, [None], [None, None], [None, None, None]]
        if self._smspy:
            execfile(self._smspy, localSms)
        else:
            return None
        if not expected:
            if 'Dict' in localSms and self._txName in localSms['Dict']:
                if localSms['Dict'][self._txName] in fakeDicts:
                    logger.warning('No useful upper limit dictionary was found' + \
                                   'for plot %s.' %self.name)
                    return None
                return localSms['Dict'][self._txName]
            logger.warning('No upper limit dictionary was found for plot %s.' %self.name)
            return None
        if expected:
            if 'ExpectedDict' in localSms and self._txName in \
            localSms['ExpectedDict']:
                if localSms['ExpectedDict'][self._txName] in fakeDicts:
                    logger.warning('No useful expected upper limit dictionary' + \
                                   'was found for plot %s.' %self.name)
                    return None
                return localSms['ExpectedDict'][self._txName]    
            logger.warning('No expected upper limit dictionary was found' + \
                           'for plot %s.' %self.name)  
            return None
            
    @property
    def _axes(self):
        if self._txObject.name == self._txName:
            ax = self._experimentIDObject._infotxt.axes[self._txObject.name][0]
            a = {key: ax[key] for key in ax if key != 'extension'}
            return a
        for ax in self._experimentIDObject._infotxt.axes[self._txObject.name]:
            #### FIXME what about extensions = None for tx without intermediate?
            if self._txObject.name + ax['extension'] != self._txName:
                continue
            a = {key: ax[key] for key in ax if key != 'extension'}
            return a
