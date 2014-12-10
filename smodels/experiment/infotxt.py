"""
.. module:: infotxt
   :synopsis: Holds the info.txt object

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>


"""

import logging
from databaseBrowserException import InvalidInfotxtFileException
from smodels.tools.physicsUnits import GeV

FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.ERROR)


class Infotxt(object):
    """Holds all the lines, stored in the info.txt file. 
    Provides the required information about txNames, results and all the 
    meta-information needed for the experimental objects.

    """
    
    def __init__(self, path):
        
        self._path = path
        logger.debug('Creating object based on info.txt: %s' %self._path)
        self._exceptions = ['constraint', 'condition', 'fuzzycondition', \
        'unconstraint', 'exclusions', 'expectedexclusions', 'exclusionsp1', \
        'expectedexclusionsp1','exclusionsm1', 'expectedexclusionsm1', \
        'category', 'branchcondition']
        try:
            self._experimentID = self.metaInfo['id']
        except KeyError:
            try:
                self._experimentID = self.metaInfo['pas']
            except KeyError:
                raise InvalidInfotxtFileException('no id or pas found')
        self._run = self.metaInfo['sqrts']
        self._dummies = ['not yet assigned', 'not jet assigned', '', 'not yet assignet']
    
    @property
    def _readInfo(self):
        """Reads the whole info.txt file, returns a tuple containing a 
        dictionary holding the meta-information (e.g. PAS, lumi, comment, ...) 
        and a list holding all the lines with keywords that show up several 
        times (e.g. 'constraint', 'condition', 'exclusions', ...).
    
        """
        info = []
        infoFile = open(self._path)
        content = infoFile.readlines()
        infoFile.close()
        
        #logger.debug('Found info.txt for %s.' %self._path)
        metaInfo = {line.split(':', 1)[0].strip(): line.split(':', 1)[1].strip() \
        for line in content if not line.split(':')[0].strip() in self._exceptions}
        #logger.debug('Meta info is %s' %metaInfo)
        for key in self._exceptions:
            for line in content:
                if key in line:
                    info.append(line.strip())
        #logger.debug('Info is %s' %info)    
        return [metaInfo, info]
    
    @property
    def metaInfo(self):
        """Returns the meta info dictionary (contains the axes line too).
        
        """
        return self._readInfo[0]
        
    @property
    def info(self):
        """Returns the list of lines connected with the txNames.
        
        """
        return self._readInfo[1]
    
    def _txNameInfo(self, requested):
        """Creates a dictionary for txName related information.
        :return: {'txName': 'line of info.txt'}
        
        """
        dic = {}
        logger.debug('Look for requested keyword %s.' %requested)
        content = self.info
        content = [string.strip() for string in content if requested in string \
        and not 'assigned' in string]
        content = [string.split(':')[1] for string in content] 
        for c in content:
            dic[c.split('->')[0].strip()] = c.split('->')[1].strip()
        return dic
        
    @property
    def category(self):
        cat = self._txNameInfo('category')
        return cat
    
    @property
    def constraints(self):
        const = self._txNameInfo('constraint')
        return const
    
    # takes every line holding constraint and unconstraint
    #@property
    #def txNames(self):
        #txNames = []
        #content = [c for c in self.info if not 'assigned' in c]
        #content = [string.strip() for string in content if 'constraint' \
        #in string or 'unconstraint' in string]
        
        #for c in content:
            #if not c.split(' ')[1] in txNames:
                #txNames.append(c.split(' ')[1])                
        #return txNames

    # gives only properly assigned constraint lines and avoids T7 and T8 txNames
    @property
    def txNames(self):
        txNames = []
        for key in self._txNameInfo('constraint'):
            if '8' in key or '7' in key: continue
            content = self._txNameInfo('constraint')[key]
            if content in self._dummies: continue
            if not key in txNames:
                txNames.append(key)
        for key in self._txNameInfo('unconstraint'):
            if '8' in key or '7' in key: continue
            content = self._txNameInfo('unconstraint')[key]
            if not key in txNames:
                txNames.append(key)
        return txNames        
                
    @property    
    def exclusions(self):
        """Retrieves all the exclusions for every extended txName stored in 
        the info.txt and returns them as simple list.
        
        """
        infList = self.info
        exList = []
        keys = ['exclusions', 'expectedexclusions', 'exclusionsp1', \
        'expectedexclusionsp1','exclusionsm1', 'expectedexclusionsm1']
        infList = [l for l in infList for k in keys if k in l]
        for l in infList:
            if exList.count(l) == 0:
                exList.append(l)
        logger.debug('List of exclusions for %s-%s: %s.' \
        %(self._run, self._experimentID, exList))
        return exList
        
    @property    
    def _preprocessAxesLine(self):
        """Handles the information stored in the axes-labeled line of info.txt, 
        therefor this line has to be preprocessed.
    
        """
        if not 'axes' in self.metaInfo:
            logger.warning('There is no axes information for %s!' %self._experimentID)
            return None
        infoLine = self.metaInfo['axes'].split(',')
        infoLine = [ax.strip() for ax in infoLine]
        logger.debug('axes-information: %s' %infoLine)
        return infoLine
        
    def _axesDict(self, axesLines):
        """Splits the axes line and retrieves all the txNames to form
        a dictionary.
        :return: {'txName': ['axes entry', 'axes entry', ...]}
    
        """
        axDic = {}
        for axesLine in axesLines:
            txName = axesLine.split(' ')[0].replace(':', '').strip()
            axDic[txName] = axesLine.replace(txName + ':', '').split('-')
            axDic[txName] = [c.strip() for c in axDic[txName]]
            logger.debug('For %s there are %s masses.' \
            %(txName, len(axDic[txName])))
            logger.debug('For %s the axes dictionary is: %s.' %(txName, axDic[txName]))
        return axDic
        
    
    def _massDict(self, axesEntry):
        """Retrieves the axes information for given entry as dictionary.
        :param axesEntry: one axes entry for one txName
        :return:  {'mx': mass on x-axis, 'my': mass on y-axes, 
        'mz': condition for intermediate mass} or 
        {'mx': mass on x-axis, 'my': mass on y-axes, 
        'm3': condition for first intermediate mass,
        'm4': condition for second intermediate mass}
        
        """
        
        axDict = {}
        
        logger.debug('Axes entry: %s.' %axesEntry.split())
        axDict['mx'] = axesEntry.split()[0].strip()
        axDict['my'] = axesEntry.split()[1].strip()
        try:
            axesEntry.split()[3]
            logger.info('There are more then three masses!'
                        'Keys will be mx, my, m3 and m4!')
            axDict['m3'] = axesEntry.split()[2].strip()
            axDict['m4'] = axesEntry.split()[3].strip()
            axDict['extension'] = '%s%s' %(axDict['m3'], axDict['m4'])
        except IndexError:
            try:
                axDict['mz'] = axesEntry.split()[2].strip()
                axDict['extension'] = axesEntry.split()[2].strip()
            except IndexError:
                logger.debug('No intermediate mass mz.')
                axDict['mz'] = None
                axDict['extension'] = None
        if axDict['extension'] and 'D' in axDict['extension']:
            axDict['extension'] = 'D' + axDict['extension'].split('=')[-1].strip()
        return axDict
        
    def _massCondition(self, mz):
        """Takes the axes entry for the third mass and splits it into
        condition for this mass and its value with units added if not unitless:
        -) fixed LSP, Chargino or other mass in GeV
        -) mass splitting (as well 'xvalue'): M2=x*M1+(1-x)*M0
        -) difference between masses (e.g. M1-M0) in GeV
        -) ratio between masses (e.g. M2/M0) in percent
        :param mz: third item of an axes-entry for one txName
        :return: ('massCondition', value)
        
        """
        
        try:
            value = float(mz) / (10. ** (len(mz)-1))
            condition = 'massSplitting'
        except TypeError:
            logger.debug('Got no mz!')
            return (None, None)
        except ValueError:
            if 'D' in mz:
                #value = addunit(int(mz.split('=')[-1].strip()), 'GeV')
                value = int(mz.split('=')[-1].strip()) * GeV
                if mz.split('=')[0].strip() == 'D':
                    logger.error('There is something wrong with the "D-entry"!\n \
                    Check database for %s!' %self.name)
                    condition = 'unknownDifference' 
                else:
                    condition = mz.split('=')[0].strip().replace('D(', '')
                    condition = condition.replace(')', '')
                    condition = condition.split('/')
                    condition = '%s-%s' %(condition[0], condition[1])
            elif 'LSP' in mz:
                #value = addunit(int(mz.replace('LSP', '')), 'GeV')
                value = int(mz.replace('LSP', '')) * GeV
                condition = 'fixedLSP'
            elif 'M' in mz:
                #value = addunit(int(mz[2:]), 'GeV')
                value = int(mz[2:]) * GeV
                condition = 'fixed%s' %mz[:2]
            elif 'C' in mz:
                #value = addunit(int(mz.replace('C', '')), 'GeV')
                value = int(mz.replace('C', '')) * GeV
                condition = 'fixedM2'
            elif 'x' in mz:
                value = (float(mz.replace('x', ''))/100.)
                condition = 'M2/M0'
            else:
                logger.error('Unknown third mass entry %s!' %mz)
                value = None
                condition = None
        return (condition, value)
                        
    @property
    def axes(self):
        """Runs all the preprocessing methods and retrieves the wrought axes 
        information.
        :return: {'txName': [{'mx': 'mass on x-axis', 'my': 'mass on y-axis',
        'mz': ('condition for third mass', int(value for this condition)), 
        'extension': 'extension glued to txName name'}]}
  
        """
        
        axLines = self._preprocessAxesLine
        if not axLines:
            return None
        axDict = self._axesDict(axLines)
        # Just to cross check these two fields of the info.txt,
        # print some warnings.
        for t in self.txNames:
            if not t in axDict:
                logger.warning('There is no axes entry for %s-%s! Check database!' \
                %(self._experimentID, t))
        for t in axDict:
            if not t in self.txNames:
                logger.debug('There is an axes entry for %s-%s' %(self._experimentID, t) +
                             ', but this is no known txName! Check database!' \
                             )
        axDict = {t: axDict[t] for t in axDict if t in self.txNames}
        for t in axDict:
            entries = []
            for entry in axDict[t]:
                entries.append(self._massDict(entry))

            for entry in entries:
                for key in ['mz', 'm3', 'm4']:
                    if key in entry:
                        entry[key] = self._massCondition(entry[key])
                    
            axDict[t] = entries    
            logger.debug('Axes information for %s is: %s' \
            %(t, axDict[t]))
        return axDict
    
    @property
    def extensions(self):
        """Retrieves a dictionary with the txNames as keys and the 
        extensions for these txNames to find exclusion lines etc.
        :return: {'txName': ['extension']}
        
        """
        if not self.axes:
            return None
        extDic = {}
        for t in self.txNames:
            extDic[t] = []
            try:
                axes = self.axes[t]
            except KeyError:
                continue
            for ax in axes:
                if ax['extension']:
                    extDic[t].append(ax['extension'])
                        
        return extDic