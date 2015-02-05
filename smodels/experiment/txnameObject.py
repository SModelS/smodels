"""
.. module:: dataObjects
   :synopsis: Holds the classes and methods used to read and store the information in the
              sms.py files. Also contains the interpolation methods

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>


"""

import logging,os,sys
from smodels.tools.physicsUnits import GeV, fb, TeV, pb
from smodels.theory.particleNames import elementsInStr
from smodels.theory.element import Element

FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.DEBUG)


class TxName(object):
    """Holds the information related to one txname in the Txname.txt
    file (constraint, condition,...) as well as the data.
    """
    
    def __init__(self, path):
        self.txnameFile = path
        self.txnameData = None
        self._elements = []
        
        logger.debug('Creating object based on txname file: %s' %self.txnameFile)        
 
        #Open the info file and get the information:
        if not os.path.isfile(path):
            logger.error("Txname file %s not found" % path)
            sys.exit()      
        txfile = open(self.txnameFile)
        content = txfile.readlines()
        txfile.close()
        
        #Get tags in info file:
        tags = [line.split(':', 1)[0].strip() for line in content]
        for i,tag in enumerate(tags):
            if not tag: continue
            line = content[i]
            value = line.split(':',1)[1].strip()            
            if tags.count(tag) == 1:
                if ';' in value: value = value.split(';')
                if tag == 'upperLimits' or tag == 'efficiencyMap':
                    self.txnameData = TxNameData(tag,value)
                else: self.addInfo(tag,value)
            else:
                logger.info("Ignoring unknown field %s found in file %s" % (tag, self.infopath))
                continue
        
        #Builds up a list of _elements appearing in constraints:        
        if hasattr(self,'constraint'):
            self._elements = [Element(el) for el in elementsInStr(self.constraint)]

    def __str__(self):
        return self.txname
        
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

    def getEfficiencyFor(self,element):
        """
        For upper limit results, checks if the input element appears in the constraints
        and falls inside the upper limit grid.
        If it does, returns efficiency = 1, else returns efficiency = 0.
        For efficiency map results, checks if the input element appears in the constraints.
        and falls inside the efficiency map grid.
        If it does, returns the corresponding efficiency value, else returns efficiency = 0.
        
        :param element: Element object
        :return: efficiency (float)
        """
        
        if self.txnameData.type == 'upperLimits':
            for el in self._elements:                
                if element.particlesMatch(el):
                    ul = self.txnameData.getValueFor(el.getMasses())
                    if type(ul) == type(fb): return 1.
            return 0.
        elif self.txnameData.type == 'efficiencyMap':
            for el in self._elements:
                if element.particlesMatch(el):
                    eff = self.txnameData.getValueFor(element.getMasses())
                    if eff: return eff                    
            return 0.
        else:
            logger.error("Unknown data type: %s" % self.txnameData.type)
            sys.exit()
            
        


class TxNameData(object):
    """Holds the data for the Txname object.  It holds Upper limit values or efficiencies."""
    
    def __init__(self,tag,value):
        self.type = tag
        self.data = eval(value, {'fb' : fb, 'pb' : pb, 'GeV' : GeV, 'TeV' : TeV})
       
    def getValueFor(self,massarray):
        """
        Interpolates the data and returns the UL or efficiency for the respective massarray
        :param massarray: mass array values (with units), i.e. [[100*GeV,10*GeV],[100*GeV,10*GeV]]
        """
        
        return 1.*fb
    
