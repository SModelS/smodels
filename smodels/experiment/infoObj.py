"""
.. module:: infoObj
   :synopsis: Holds the classes and methods used to read and store the
              information in the globalInfo.txt or dataInfo.txt files.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>


"""

import os,sys
from smodels.tools.physicsUnits import GeV, fb, TeV, pb
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from smodels.tools.smodelsLogging import logger

class Info(object):
    """
    Holds the meta data information contained in a .txt file
    (luminosity, sqrts, experimentID,...). Its attributes are generated according to the lines in the
    .txt file which contain "info_tag: value".
    """

    def __init__(self, path=None):
        """
        :param path: path to the .txt file
        """

        self.path = path
        if path:
            logger.debug('Creating object based on  %s' %self.path)

            #Open the info file and get the information:
            if not os.path.isfile(path):
                logger.error("Info file %s not found" % path)
                raise SModelSError()
            from smodels.tools.stringTools import concatenateLines
            infoFile = open(self.path)
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
                    logger.info("Ignoring unknown field %s found in file %s"
                                % (tag, self.path))
                    continue

            self.cacheJsons()

    def __eq__ ( self, other ):
        if self.__dict__ != other.__dict__:
            return False
        return True

    def cacheJsons ( self ):
        """ if we have the "jsonFiles" attribute defined,
            we cache the corresponding jsons. Needed when pickling """
        if not hasattr ( self, "jsonFiles" ):
            return
        if hasattr ( self, "jsons" ): ## seems like we already have them
            return
        import json
        self.jsons = list()
        dirp = os.path.dirname ( self.path )
        jsonFiles = [os.path.join( dirp, js) for js in self.jsonFiles]
        for js in jsonFiles:
            with open(js, "r") as fi:
                self.jsons.append(json.load(fi))


    def dirName ( self, up=0 ):
        """ directory name of path. If up>0,
            we step up 'up' directory levels.
        """
        s_up = "/".join ( [ ".." ] * up )
        p = os.path.dirname ( self.path )
        return os.path.abspath ( os.path.join ( p, s_up ) )

    def __ne__ ( self, other ):
        return not self.__eq__ ( other )

    def addInfo(self,tag,value):
        """
        Adds the info field labeled by tag with value value to the object.

        :param tag: information label (string)
        :param value: value for the field in string format
        """
        if tag == "lastUpdate": # dont eval that!
            setattr ( self, "lastUpdate", str(value) )
            return
        try:
            setattr(self,tag,eval(value, {'fb':fb, 'pb':pb, 'GeV':GeV, 'TeV':TeV}))
        except SyntaxError:
            setattr(self,tag,value)
        except NameError:
            setattr(self,tag,value)
        except TypeError:
            setattr(self,tag,value)

    def getInfo(self, infoLabel):
        """
        Returns the value of info field.

        :param infoLabel: label of the info field (string). It must be an attribute
                          of the GlobalInfo object
        """

        if hasattr(self,infoLabel): return getattr(self,infoLabel)
        else: return False
