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
from typing import Optional, Union
from smodels.base.types import PathType

class Info(object):
    """
    Holds the meta data information contained in a .txt file
    (luminosity, sqrts, experimentID,...).
    Its attributes are generated according to the lines in the
    .txt file which contain "info_tag: value".
    """

    def canonizeRegions ( self, regions : Optional[dict] = None ) \
            -> Optional[dict]:
        """ given a list of regions in globalInfo.txt in the
        regionMappings field,
        return a canonical version of that list: strings in
        that list get transformed into dictionaries, if region type is
        missing, "SR" is assumed. if the "smodels" counterpart is not
        given for a region, we assume that there is None.
        If no pyhf name is given, we assume it to be the smodels name.
        if no onnx name is given, we assume it to be the pyhf name.
        if no label is given, we assume it to be the smodels name.
        :param regions: list of regions in regionMappings globalInfo.txt
        :returns: canonical list of regions
        """
        if regions is None:
            return regions
        newregions={}
        for label,region in regions.items():
            regionDict = {}
            if type(region)==str:
                logger.debug ( f"exploding {region} to a dictionary in {self.path}" )
                regionDict["smodels"] = region
            elif type(region)==dict:
                regionDict.update(region)
            else:
                raise SModelSError ( f"region {region} is neither a string nor a dictionary in {self.path}" )
            for key in [ "type", "smodels", "pyhf", "sl", "onnx" ]:
                if key not in regionDict:
                    logger.debug ( f"region {label} has no {key} defined, setting to default value in {self.path}.")
            # Set defaults (if they were not defined)
            regionDict.setdefault("type", "SR")
            if regionDict["type"]=="SR":
                regionDict.setdefault("smodels", label)
            else:
                regionDict.setdefault("smodels", None)
            
            regionDict.setdefault("pyhf", regionDict["smodels"])
            regionDict.setdefault("sl", regionDict["smodels"])
            regionDict.setdefault("onnx", regionDict["pyhf"])
            
            newregions[label] = regionDict
        return newregions

    def __init__(self, path : Optional [ PathType ] = None ):
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
                if tag in [ "covariance", "jsonFiles", \
                            "jsonFiles_FullLikelihood", "datasetOrder" ]:
                    from smodels import installation
                    logger.warning ( f"Tag '{tag}' in {self.path} was used in the old format of the database. This version {installation.version()} of SModelS ignores this tag. Use 'regionMappings', 'regionSets', 'statModels' instead. See the documentation -> 'Detailed Guide to SModelS' -> 'Database of Experimental Results' for more info." )
                if tag == "regionMappings":
                    try:
                        regions = eval(value)
                    except Exception as e:
                        try:
                            from pygments import highlight
                            from pygments.lexers import JsonLexer
                            from pygments.formatters import TerminalFormatter
                            print(highlight(value, JsonLexer(), TerminalFormatter()))
                        finally:
                            raise e
                    regions = self.canonizeRegions( regions )
                    value = str(regions)
                if tags.count(tag) == 1:
                    self.addInfo(tag, value)
                else:
                    logger.info(f"tag {tag} given multiple times in {self.path}" )
                    continue
            
            self.checkConsistencyOfRegionSets()
            self.checkConsistencyOfStatModels()
            self.cacheStatsModels()

    def checkConsistencyOfRegionSets( self ):
        """ check that all SRs mentioned in regionSets are included in regionMappings """
        if not hasattr ( self, "regionMappings" ):
            return
        if not hasattr ( self, "regionSets" ):
            return
        elif not isinstance ( self.regionSets, dict ):
            raise SModelSError ( f"regionSets has to be a dict, but is {type(self.regionSets)}" )
        
        for regionList in self.regionSets.values():
            if not isinstance ( regionList, list ):
                raise SModelSError ( f"regionSets has to be a dict of lists, but its values are of type {type(regionList)}" )
            for region in regionList:
                if region not in  self.regionMappings:
                    raise SModelSError ( f"region label {region} not mentioned in regionMappings" )
                
    def checkConsistencyOfStatModels( self ):
        """ check that all stat models mentioned in statModels are included in regionSets """
        if not hasattr ( self, "statModels" ):
            return
        if not hasattr ( self, "regionSets" ):
            raise SModelSError ( f"In {self.path}: statModels is defined, but regionSets is not. This is inconsistent." )
        elif not isinstance ( self.statModels, dict ):
            raise SModelSError ( f"In {self.path}: statModels has to be a dict, but is {type(self.statModels)}" )
        
        for regionSet, model_tuples in self.statModels.items():
            if regionSet not in self.regionSets:
                raise SModelSError ( f"In {self.path}: region set {regionSet} mentioned in statModels is not defined in regionSets" )
            if not isinstance ( model_tuples, list ):
                raise SModelSError ( f"In {self.path}: statModels has to be a dict of lists, but its values are of type {type(model_tuples)}" )
            for model_tuple in model_tuples:
                if not isinstance ( model_tuple, tuple ):
                    raise SModelSError ( f"In {self.path}: statModels has to be a dict of lists of tuples, but its values are of type {type(model_tuple)}" )
                if len(model_tuple) != 2:
                    raise SModelSError ( f"In {self.path}: statModels has to be a dict of lists of tuples of length 2, but its values are of length {len(model_tuple)}" )
                model_type,model = model_tuple
                if model_type not in [ "sl", "full_pyhf", "pyhf", "onnx" ]:
                    raise SModelSError ( f"In {self.path}: model_type {model_type} is unknown. should be of: onnx, pyhf, full_pyhf, sl" )
                if not isinstance ( model, str ):
                    raise SModelSError ( f"In {self.path}: model has to be a string, but is {type(model)}" )
                if not os.path.isfile ( os.path.join(os.path.dirname(self.path), model) ):
                    raise SModelSError ( f"In {self.path}: model file {model} does not exist in the same directory as {self.path}" )

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
        for regionSetName, model_tuples in self.statModels.items():
            assert regionSetName in self.regionSets, f"{regionSetName} is not in {self.regionSets}"
            n_regionset = len( self.regionSets[regionSetName] )
            for model_tuple in model_tuples:
                model_type, model = model_tuple
                fullPath = os.path.join(dirp, model )
                with open ( fullPath, "rb" ) as f:
                    txt = None
                    if model_type == "sl":
                        with open(fullPath,"rt") as f:
                            txt = eval ( f.read() )
                            n_rows = len(txt)
                            if n_regionset != n_rows:
                                raise SModelSError ( f"dimension({n_regionset}) of region set {regionSetName} does not match number of rows({n_rows}) of covariance matrix {model} in {self.path}" )
                            n_cols = len(txt[0])
                            if n_regionset != n_cols:
                                raise SModelSError ( f"dimension({n_regionset}) of region set {regionSetName} does not match number of columns({n_columns}) of covariance matrix {model} in {self.path}" )
                    elif model_type in [ "full_pyhf", "pyhf" ]:
                        with open(fullPath,"rt") as f:
                            txt = json.load(f)
                    elif model_type == "onnx":
                        txt = f.read()
                    else:
                        logger.error ( f"model_type {model_type} is unknown. should be of: onnx, pyhf, full_pyhf, sl" )
                    self.cachedModels[model] = txt
                    f.close()

    def dirName(self, up : int = 0 ) -> str:
        """ directory name of path. If up>0,
            we step up 'up' directory levels.
        """
        s_up = "/".join([".."] * up)
        p = os.path.dirname(self.path)
        return os.path.abspath(os.path.join(p, s_up))

    def __ne__(self, other):
        return not self.__eq__(other)

    def addInfo(self, tag : str, value : str ):
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

    def getInfo(self, infoLabel : str ) -> bool | str:
        """
        Returns the value of info field.

        :param infoLabel: label of the info field (string). 
        It must be an attribute of the GlobalInfo object
        """

        if hasattr(self, infoLabel):
            return getattr(self, infoLabel)
        else:
            return False
