"""
.. module:: infoObjects
   :synopsis: Holds the classes and methods used to read and store the information in the
              info.txt files.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>


"""

import logging,os,sys
from smodels.tools.physicsUnits import GeV, fb, TeV, pb

FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.ERROR)


class TxNameInfo(object):
    """Holds the information related to one txname in the info.txt
    file (constraint, condition,...).
    Its attributes are generated according to the lines in the
    info.txt file which contain "info_tag: Txname -> value".
    
    :ivar name: Txname (string)    
    """
    
    def __init__(self, txname=None):
        
        self.txname = txname    
        
    def addInfo(self,tag,value):
        """
        Adds the info field labeled by tag with value value to the object.
        :param tag: information label (string)
        :param value: value for the field in string format 
        """

        setattr(self,tag,value)
    
    def getInfo(self, infoLabel):
        """Returns the value of info field.
        :param infoLabel: label of the info field (string). It must be an attribute of
                          the TxNameInfo object
        """
        
        if hasattr(self,infoLabel): return getattr(self,infoLabel)
        else: return False
    

class GlobalInfo(object):
    """
    Holds the global information contained in a info.txt file
    (luminosity, sqrts, experimentID,...).
    Its attributes are generated according to the lines in the
    info.txt file which contain "info_tag: value".
    
    :ivar infofile: path to the info.txt file
    :ivar analysisType: default type for analysis
    :ivar numericalAttr: list of properties which should be evaluated
                         when added to the object
    """
    
    def __init__(self, path):        
        self.infofile = path
        self.analysisType = 'UpperLimit'
        self._numericalAttr = ['sqrts','lumi']
        
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
                logger.error("The value for %s should be numerical" % tag)
                sys.exit()
        else: setattr(self,tag,value)
        
    def getInfo(self, infoLabel):
        """Returns the value of info field.
        :param infoLabel: label of the info field (string). It must be an attribute of
                          the GlobalInfo object
        """
        
        if hasattr(self,infoLabel): return getattr(self,infoLabel)
        else: return False
    
class InfoFile(object):
    """Holds all the information stored in the info.txt file. 
    Provides the required information about txNames, results and all the 
    meta-information needed for a single analysis object.
       
    :ivar infopath: path to the info.txt file
    :ivar globalInfo: GlobalInfo object. Containts all the global information in the file
                      (lum, sqrts, ID,...)
    :ivar txNameInfoList: a list of TxNameInfo objects constaining the information
                         which is txname-specific (constraint, condition,...)
    """
    
    def __init__(self, path):
        """
        Reads the whole info.txt file and builds the GlobalInfo and txNameInfo
        objects. The attributes of GlobalInfo and TxNameInfo are setting according
        to the info labels in the file.
        
        :param path: path to the info.txt file    
        """
        
        self.infopath = path
        logger.debug('Creating object based on info.txt: %s' %self.infopath)        
        self.globalInfo = None
        self.txNameInfoList = []        
 
        #Open the info file and get the information:
        if not os.path.isfile(path):
            logger.error("Info file %s not found" % path)
            sys.exit()            
        infoFile = open(self.infopath)
        content = infoFile.readlines()
        infoFile.close()
        globalInfo = GlobalInfo(self.infopath)
        txObjects = {}
        
        #Get tags in info file:
        tags = [line.split(':', 1)[0].strip() for line in content]
        for i,tag in enumerate(tags):
            if not tag: continue
            line = content[i]
            value = line.split(':',1)[1].strip()            
            if tags.count(tag) == 1 and not "->" in value:
                globalInfo.addInfo(tag,value)
            elif "->" in value:
                txname = value.split('->')[0].strip()
                newvalue = value.split('->')[1]
                if not txname in txObjects: txObjects[txname] = TxNameInfo(txname)
                txObjects[txname].addInfo(tag,newvalue)
            else:
                logger.info("Ignoring unknown field %s found in file %s" % (tag, self.infopath))
                continue
                

        self.globalInfo = globalInfo
        self.txNameInfoList = txObjects.values()


    def getInfo(self, infoLabel, txname=None):
        """Returns the value of info field.
        If the info field belongs to a Txname (TxNameInfo object) and txname
        is specified, returns the value for the respective TxName. Otherwise returns
        a dictionary with the values for all TxNames containing the info field.
        :param infoLabel: label of the info field (string). It must be an attribute of
                          the GlobalInfo object or one of the TxNameInfo objects
        :param txname: If infoLabel belongs to a TxnameInfo object, the desired TxName must
                      be specified (string)
        """
        
        globInfo = self.globalInfo.getInfo(infoLabel)
        if not globInfo is False: return globInfo
        txInfo = False
        if txname:
            for txNameInfo in self.txNameInfoList:           
                if txNameInfo.txname == txname:
                    txInfo = txNameInfo.getInfo(infoLabel)
                    break
            if not txInfo is False: return txInfo
        else:
            infoList = {}
            for txNameInfo in self.txNameInfoList:
                txInfo = txNameInfo.getInfo(infoLabel)
                if txInfo: infoList[txNameInfo.txname] = txInfo
            if infoList: return infoList        
        
        logger.error("Info field %s not found" % infoLabel)
        sys.exit()
  
    
    def getTxNames(self):
        """Returns a list of all the txnames contained in the info.txt file.
        
        :return: list of txnames (strings)        
        """
        
        txnames = [txNameInfo.txname for txNameInfo in self.txNameInfoList]
        return txnames

    def getTxInfoFor(self,txname):
        """Returns the TxNameInfo object for the corresponding txname
        :return: TxNameInfo object      
        """
        
        for txInfo in self.txNameInfoList:
            if txInfo.txname == txname: return txInfo
        logger.info("No object found for Txname = %s" % txname)
        return None