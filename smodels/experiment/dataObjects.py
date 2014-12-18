"""
.. module:: dataObjects
   :synopsis: Holds the classes and methods used to read and store the information in the
              sms.py files. Also contains the interpolation methods

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>


"""

import logging,os,sys
from smodels.tools.physicsUnits import GeV, fb, TeV, pb
from smodels.theory import element
import numpy as np
from scipy.interpolate import griddata

FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.INFO)


class ULdata(object):
    """Holds a list of points with upper limits for the respective Txname
    and analysis ID: [[massarray_1,UL_1],[massarray_2,UL_2],...]

    :ivar analysisID: analysis ID (CMS-SUS-13-006,...)
    :ivar txname: the respective Txname (TChiWZ,...)
    :ivar data: list of data point
    """
    
    def __init__(self, ID, txname):
        
        self.analysisID = ID
        self.txname = txname
        self.data = None
        
    def getULFor(self,massarray):
        """
        Interpolates the data and returns the UL for the respective massarray
        :param massarray: mass array values (with units), i.e. [[100*GeV,10*GeV],[100*GeV,10*GeV]]
        """    
        
        if np.array(self.data[0][0]).shape() != np.array(massarray).shape():
            logger.error("Wrong mass input format")
            sys.exit()

        mass = np.array(massarray).flatten()
        xpts, ypts = np.array(), np.array()
        for pt in self.data:            
            xpts.append(np.array(pt[0]).flatten())
            ypts.append(pt[1])

        ul = griddata(xpts, ypts, [mass], method="linear")
        #Deal with nested result from griddata:
        while isinstance(ul[0],np.ndarray): ul = ul.flatten()
        ul = ul[0]
    
        if type(ul) != type(pb):
            logger.error("Masses out of efficiency map range (no extrapolation)")
            return 0.
        else: return ul            
            

class EMdata(object):
    """Holds a dictionary with elements as keys and the efficiency map grid
    as values: {Element : [list of points],...}
    
    :ivar analysisID: a unique analysis ID (CMS-SUS-13-006:SR1,...)
    :ivar dataDic: dictionary with  data points
    :ivar obsEvents: number of observed events
    :ivar expectedBG: expected number of BG events
    :ivar expBGerror: systematical error for the expected number of BG events
    """
    
    def __init__(self, ID):
        
        self.analysisID = ID
        self.dataDic = {}
        self.obsEvents = None
        self.expectedBG = None
        self.expBGerror = None
        
    def getEffFor(self,element):
        """
        Interpolates the data and returns the efficiency for the respective element
        :param element: a Element() object
        """    
        
        for el, effmap in enumerate(self.dataDic):
            if el == element:
                return "FIX: needs to implement interpolation of effmap"

    
class DataFile(object):
    """Holds all the information stored in the sms.py file. 
    Provides the required information about txNames and data.
       
    :ivar path: path to the sms.py file
    :ivar dataList: list of data objects (either
    :ivar txNameInfoList: a list of TxNameInfo objects constaining the information
                         which is txname-specific (constraint, condition,...)
    """
    
    def __init__(self, path, infoObject):
        """
        Loads the sms.py file and generate the data objects to hold the data.
        For UL analyses, generates a list of ULdata objects, where each object
        holds the UL data points for a single Txname.
        For EM analyses, generates a single EMdata object, which holds all the
        efficiency maps (for all the elements) for a single analysis/signal region.
        
        :param datapath: path to the sms.py file (string)
        :param infoObject: InfoData object containing the respective info.txt information    
        """
        
        self.datapath = path
        self.infoObject = infoObject
        logger.debug('Creating object based on sms.py: %s' %self.datapath)        
        self.dataList = []
 
        #Open the info file and get the information:
        if not os.path.isfile(path):
            logger.error("Data file %s not found" % path)
            sys.exit()
        smspy = {}
        execfile(path,smspy)
        if not 'Dict' in smspy:
            logger.error("Data Dict not found in %" % path)
            sys.exit()
        
        #Generates the data objects
        dataDict = smspy['Dict']
        analysisID = self.infoObject.getInfo('id')
        if self.infoObject.getInfo('analysisType') == 'UpperLimit':            
            for txname,data  in dataDict.items():
                uldata = ULdata(analysisID,txname)
                uldata.data = data
                self.dataList.append(uldata)
        elif self.infoObject.getInfo('analysisType') == 'EfficiencyMap':
            emdata = EMdata(analysisID)
            for txname,data in dataDict.items():
                txInfo = self.infoObject.getTxInfoFor(txname)
                el = element.Element(txInfo.constraint)
                emdata.dataDic[el] = data
            self.dataList.append(emdata)
        else:
            logger.warning("Unknown analysis type: %s " % self.infoObject.getInfo('analysisType'))
            return False
        
    
    def getTxNames(self):
        """Returns a list of all the txnames contained in the smspy file.
        
        :return: list of txnames (strings)        
        """
        
        smspy = {}
        execfile(self.datapath,smspy)
        txnames = smspy['Dict'].keys()
        return txnames

    def getData(self,txname=None):
        """For UL analyses returns the list of the ULdata object for the corresponding txname.
        For EM-type analyses, returns the EMdata object (txname is not used).
        
        :return: ULdata or EMdata object      
        """
        
        
        if self.infoObject.getInfo('analysisType') == 'EfficiencyMap': return self.dataList[0]
        elif self.infoObject.getInfo('analysisType') == 'UpperLimit':
            for uldata in self.dataList:
                if txname == uldata.txname: return uldata
       
        logger.info("No data found for Txname = %s " % txname)
        return None