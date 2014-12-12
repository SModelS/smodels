"""
.. module:: infotxt
   :synopsis: Holds the classes and methods used to generate InfoTxt objects.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>


"""

import logging
from databaseBrowserException import InvalidInfoFieldException,InvalidFieldValueException
from smodels.tools.physicsUnits import GeV, fb, TeV, pb

FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.ERROR)


class TxNameInfo(object):
    
    def __init__(self, txname=None):
        
        self.name = txname    
        
    def addInfo(self,tag,value):
        """
        Adds the info field labeled by tag with value value to the object.
        :param tag: information label (string)
        :param value: value for the field in string format 
        """

        setattr(self,tag,value)        
    

class GlobalInfo(object):
    """
    Holds the global information containing in a info.txt file
    (luminosity, sqrts, experimentID,...)
    :ivar infofile: path to the info.txt file
    :ivar numericalAttr: list of properties which should be evaluated
                         when added to the object
    """
    
    def __init__(self, path):        
        self.infofile = path
        
        self._numericalAttr = ['lum','sqrts']
        
    def addInfo(self,tag,value):
        """
        Adds the info field labeled by tag with value value to the object.
        If the attribute contains
        :param tag: information label (string)
        :param value: value for the field in string format 
        """
                
        if tag in self._numericalAttr:
            try: setattr(self,tag,eval(value))
            except NameError:
                logger.error("The value for % should be numerical" % tag)
                raise InvalidFieldValueException("Non numerical value for numerical field.")
        else: setattr(self,tag,value)
    
class Infotxt(object):
    """Holds all the information stored in the info.txt file. 
    Provides the required information about txNames, results and all the 
    meta-information needed for a single analysis object.
    
    :ivar _path: path to the info.txt file
    :ivar txname: if None, holds the information for all Txnames in the info.txt file.
                  If specified, only stores the information concerning the specific txname.    
    """
    
    def __init__(self, path):
        
        self._path = path
        logger.debug('Creating object based on info.txt: %s' %self._path)        
        self.globalInfo = None
        self.txNameInfoList = []        
        self._txnameFields = ['constraint', 'condition', 'fuzzycondition', \
        'unconstraint', 'exclusions', 'expectedexclusions', 'exclusionsp1', \
        'expectedexclusionsp1','exclusionsm1', 'expectedexclusionsm1', \
        'category', 'branchcondition']        
        self._readInfo()
    
    @property
    def _readInfo(self):
        """Reads the whole info.txt file, returns a tuple containing a 
        dictionary holding the meta-information (e.g. PAS, lumi, comment, ...) 
        and a list holding all the lines with keywords that show up more than once 
        (e.g. 'constraint', 'condition', 'exclusions', ...).    
        """
        
        infoFile = open(self._path)
        content = infoFile.readlines()
        infoFile.close()
        globalInfo = GlobalInfo(self._path)
        txObjects = {}
        
        #Get tags in info file:
        tags = [line.split(':', 1)[0].strip() for line in content]
        for i,tag in enumerate(tags):
            if not tag: continue
            line = content[i]
            value = line.split(':',1)[1].strip()            
            if tags.count(tag) == 1 and not tag in self._txnameFields:
                globalInfo.addInfo(tag,value)
            elif not tag in self._txnameFields:
                logger.info("Ignoring unknown field %s found in file %s" % (tag, self._path))
                continue
            elif "->" in value:
                txname = value.split('->')[0]
                newvalue = value.split('->')[1]
                if not txname in txObjects: txObjects[txname] = TxNameInfo(txname)
                txObjects[txname].addInfo(tag,newvalue)

        self.globalInfo = globalInfo
        self.txNameInfoList = txObjects.values()

        
    @property
    def getInfo(self,infoLabel,txname=None):
        """Returns the value of info field.
        :param infoLabel: label of the info field (string). It must be an attribute of
                          the GlobalInfo object or one of the TxNameInfo objects
        :param txname: If infoLabel belongs to a TxnameInfo object, the desired TxName must
                      be specified (string)
        """
        
        infoValue = '##not specified##'
        if hasattr(self.globalInfo,infoLabel):
            infoValue = getattr(self.globalInfo,infoLabel)
        else:
            for txNameInfo in self.txNameInfoList:
                if txNameInfo.name == txname and hasattr(txNameInfo,infoLabel):
                    infoValue = getattr(txNameInfo,infoLabel)
        
        if infoValue == '##not specified##':
            logger.error("Info field %s not found" %infoLabel)
            raise InvalidInfoFieldException("Invalid info field")
                
        return infoValue
  
    
    def getTxNames(self):
        """Returns a list of all the txnames contained in the info.txt file.
        
        :return: list of txnames (strings)        
        """
        
        txnames = [txNameInfo.name for txNameInfo in self.txNameInfoList]
        return txnames

