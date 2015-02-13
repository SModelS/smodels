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

logger.setLevel(level=logging.DEBUG)

class Info(object):
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
        self.analysisType = 'upperLimit'

        logger.debug('Creating object based on info.txt: %s' %self.infofile)        
 
        #Open the info file and get the information:
        if not os.path.isfile(path):
            logger.error("Info file %s not found" % path)
            sys.exit()      
        from smodels.tools.stringTools import concatenateLines
        infoFile = open(self.infofile)
        content = concatenateLines ( infoFile.readlines() )
        infoFile.close()
        
        #Get tags in info file:
        tags = [line.split(':', 1)[0].strip() for line in content]
        for i,tag in enumerate(tags):
            if not tag: continue
            line = content[i]
            value = line.split(':',1)[1].strip()            
            if tags.count(tag) == 1:
                self.addInfo(tag,value)
            else:
                logger.info("Ignoring unknown field %s found in file %s" % (tag, self.infofile))
                continue

        
    def addInfo(self,tag,value):
        """
        Adds the info field labeled by tag with value value to the object.
        If the attribute contains
        :param tag: information label (string)
        :param value: value for the field in string format 
        """
                  
        try:
            setattr(self,tag,eval(value, {'fb' : fb, 'pb' : pb, 'GeV' : GeV, 'TeV' : TeV}))
        except SyntaxError:          
            setattr(self,tag,value)
        except NameError:
            setattr(self,tag,value)
        
    def getInfo(self, infoLabel):
        """Returns the value of info field.
        :param infoLabel: label of the info field (string). It must be an attribute of
                          the GlobalInfo object
        """
        
        if hasattr(self,infoLabel): return getattr(self,infoLabel)
        else: return False
