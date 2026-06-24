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
from typing import Optional,Union


class Info(object):
    """
    Holds the meta data information contained in a .txt file
    (luminosity, sqrts, experimentID,...).
    Its attributes are generated according to the lines in the
    .txt file which contain "info_tag: value".
    """

    def canonizeRegions ( self, regions : Optional[dict] ) -> Optional[dict]:
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
        if regions is None:
            return regions
        newregions={}
        for label,region in regions.items():
            if type(region)==str:
                logger.debug ( f"exploding {region} to a dictionary" )
                region={"smodels": region}
            if "type" not in region:
                logger.debug ( f"no type specified for {region}, will assume SR" )
                region["type"]="SR"
            if "smodels" not in region:
                if region["type"]=="SR":
                    region["smodels"]=label
                    logger.debug ( f"no 'smodels' specified in {region}, will assume {label}" )
                else:
                    region["smodels"]=None
                    logger.debug ( f"no 'smodels' specified in {region}, will assume None" )
            if "pyhf" not in region:
                region["pyhf"]=region["smodels"]
            if "onnx" not in region:
                region["onnx"]=region["pyhf"]
            if "sl" not in region:
                region["sl"]=region["smodels"]
            newregions[label] = region
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
                    try:
                        regions = eval(value)
                    except Exception as e:
                        try:
                            from pygments import highlight
                            from pygments.lexers import PythonLexer, JsonLexer
                            from pygments.formatters import TerminalFormatter
                            print(highlight(value, JsonLexer(), TerminalFormatter()))
                        finally:
                            raise e
                    regions = self.canonizeRegions ( regions )
                    value = str(regions)
                if tags.count(tag) == 1:
                    self.addInfo(tag, value)
                else:
                    logger.info(f"tag {tag} given multiple times in {self.path}" )
                    continue
            # self.createSRMappingsDict()
            self.checkConsistencyOfStatsModels()
            self.cacheStatsModels()

    def isInSRMapping ( self, region ):
        if hasattr ( self, "srMappings" ):
            for label,m in self.srMappings.items():
                if region == label:
                    return True
            return False
        ## if there is no srMappings, we look directly
        ## in srSets
        for regions in self.srSets.values():
            if region in regions:
                return True
        return False

    def checkConsistencyOfStatsModels( self ):
        """ check that all SRs mentioned in srSets are indeed in srMappings """
        if not hasattr ( self, "srSets" ):
            return
        elif not isinstance ( self.srSets, dict ):
            raise SModelSError ( f"srSets has to be a dict, but is {type(self.srSets)}" )
        for regions in self.srSets.values():
            for region in regions:
                if not self.isInSRMapping ( region ):
                    raise SModelSError ( f"region {region} not mentioned in srMapping" )

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
        if type ( self.statModels ) == str:
            raise SModelSError ( f"in {self.id}: could not parse {self.statModels}" )
        for setName, model_tuples in self.statModels.items():
            for model_tuple in model_tuples:
                if type(model_tuple)==str:
                    raise SModelSError ( f"in {self.id}: statModels have to be tuples" )
                model_type = model_tuple[0]
                model = model_tuple[1]
                fullPath = os.path.join(dirp, model )
                with open ( fullPath, "rb" ) as f:
                    if model_type == "sl":
                        with open(fullPath,"rt") as f:
                            txt = eval ( f.read() )
                    elif model_type in [ "full_pyhf", "pyhf" ]:
                        with open(fullPath,"rt") as f:
                            txt = json.load(f)
                    elif model_type == "onnx":
                        txt = f.read()
                    else:
                        logger.error ( f"model_type {model_type} is unknown. should be of: onnx, pyhf, full_pyhf, sl" )
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
            setattr( self, tag, eval(value,
                     {'fb': fb, 'pb': pb, 'GeV': GeV, 'TeV': TeV}))
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
