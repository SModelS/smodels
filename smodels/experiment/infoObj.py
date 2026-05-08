"""
.. module:: infoObj
   :synopsis: Holds the classes and methods used to read and store the
              information in the globalInfo.txt or dataInfo.txt files.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>


"""

import os
from smodels.base.physicsUnits import GeV, fb, TeV, pb
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from smodels.base.smodelsLogging import logger
from typing import Optional


class Info(object):
    """
    Holds the meta data information contained in a .txt file
    (luminosity, sqrts, experimentID,...). Its attributes are generated according to the lines in the
    .txt file which contain "info_tag: value".
    """

    def canonizeRegions ( self, regions : Optional[list] ) -> list:
        """ given a list of regions in globalInfo.txt in the
        srMappings field,
        return a canonical version of that list: strings in
        that list get transformed into dictionaries, if region type is
        missing, "SR" is assumed. if the "smodels" counterpart is not
        given for a region, we assume that there is None.
        If no pyhf name is given, we assume it to be the smodels name.
        if no onnx name is given, we assume it to be the pyhf name.
        if no label is given, we assume it to be the smodels name.
        :param regions: list of regions in srMappings globalInfo.txt
        :returns: canonical list of regions
        """
        if regions == None:
            return regions
        newregions = []
        for region in regions:
            if type(region)==str:
                region={"smodels": region}
            if not "type" in region:
                region["type"]="SR"
            if not "smodels" in region:
                region["smodels"]=None
            if not "pyhf" in region:
                region["pyhf"]=region["smodels"]
            if not "onnx" in region:
                region["onnx"]=region["pyhf"]
            if not "label" in region:
                # if there isnt a label, we take
                # the smodels as the label
                if region["smodels"] is None:
                    region["label"]=region["smodels"]
                else:
                    region["label"]=region["pyhf"]
            newregions.append ( region )
        return newregions

    def __init__(self, path=None):
        """
        :param path: path to the .txt file
        """

        self.path = path
        if path:
            logger.debug(f'Creating object based on {self.path}')

            # Open the info file and get the information:
            if not os.path.isfile(path):
                logger.error(f"Info file {path} not found")
                raise SModelSError()
            from smodels.experiment.expAuxiliaryFuncs import concatenateLines
            infoFile = open(self.path)
            content = concatenateLines(infoFile.readlines())
            infoFile.close()

            # Get tags in info file:
            tags = [line.split(':', 1)[0].strip() for line in content]
            for i, tag in enumerate(tags):
                if not tag:
                    continue
                if tag.startswith("#"):  # a comment!
                    continue
                line = content[i]
                value = line.split(':', 1)[1].strip()
                if tag in [ "srMappings" ]:
                    regions = eval(value)
                    regions = self.canonizeRegions ( regions )
                    value = str(regions)
                if tags.count(tag) == 1:
                    self.addInfo(tag, value)
                else:
                    logger.info(f"tag {tag} given multiple times in {self.path}" )
                    continue

            self.cacheStatsModels()

    def __eq__(self, other):
        if self.__dict__ != other.__dict__:
            return False
        return True

    def cacheStatsModels(self):
        """ if we have the "statModels" attribute defined,
            we cache the corresponding onnx files. Needed when pickling """
        if not hasattr(self, "statModels"):
            ## we dont have any stats models, nothing to cache
            return
        if hasattr(self, "cachedModels"):  
            # seems like we already have cached them
            return
        import json
        self.cachedModels = {}
        dirp = os.path.dirname(self.path)
        for setName, models in self.statModels.items():
            for model in models:
                fullPath = os.path.join(dirp, model )
                with open ( fullPath, "rb" ) as f:
                    if model.endswith ( ".json" ):
                        with open(fullPath,"rt") as f:
                            txt = json.load(f)
                    elif model.endswith ( ".onnx" ):
                        txt = f.read()
                    else:
                        logger.error ( f"{model} has unrecognized file extension: should be either json or onnx" )
                    self.cachedModels[model] = txt
                    f.close()

    def dirName(self, up=0):
        """ directory name of path. If up>0,
            we step up 'up' directory levels.
        """
        s_up = "/".join([".."] * up)
        p = os.path.dirname(self.path)
        return os.path.abspath(os.path.join(p, s_up))

    def __ne__(self, other):
        return not self.__eq__(other)

    def addInfo(self, tag, value):
        """
        Adds the info field labeled by tag with value value to the object.

        :param tag: information label (string)
        :param value: value for the field in string format
        """
        if tag == "lastUpdate":  # dont eval that!
            setattr(self, "lastUpdate", str(value))
            return
        try:
            setattr(self, tag, eval(value, {'fb': fb, 'pb': pb, 'GeV': GeV, 'TeV': TeV}))
        except (SyntaxError,NameError,TypeError):
            setattr(self, tag, value)

    def getInfo(self, infoLabel):
        """
        Returns the value of info field.

        :param infoLabel: label of the info field (string). It must be an attribute
                          of the GlobalInfo object
        """

        if hasattr(self, infoLabel):
            return getattr(self, infoLabel)
        else:
            return False
