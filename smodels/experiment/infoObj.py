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


class Info(object):
    """
    Holds the meta data information contained in a .txt file
    (luminosity, sqrts, experimentID,...). Its attributes are generated according 
    to the lines in the .txt file which contain "info_tag: value".
    """

    def canonizeRegions ( self, regions : list, forNN : bool = False ) -> list:
        """ given a list of regions in globalInfo.txt in any of the 
        jsonFiles, jsonFiles_FullLikelihood, or mlModels fields,
        return a canonical version of that list: strings in
        that list get transformed into dictionaries, if region type is
        missing, "SR" is assumed. if the "smodels" counterpart is not
        given for a region, we assume that there is None.

        :param regions: list of regions in globalInfo.txt
        :param forNN: if true, then also possibly translate "pyhf" fields into
        "onnx" field
        :returns: canonical list of regions
        """
        newregions = []
        for region in regions:
            if type(region)==str:
                region={"smodels": region}
            if not "type" in region:
                region["type"]="SR"
            if not "smodels" in region:
                region["smodels"]=None
            if forNN:
                if not "onnx" in region and "pyhf" in region:
                    region["onnx"]=region["pyhf"]
                    region.pop("pyhf")
            newregions.append ( region )
        return newregions

    def __init__(self, path=None):
        """
        :param path: path to the .txt file
        """

        self.path = path
        if path:
            logger.debug('Creating object based on  %s' % self.path)

            # Open the info file and get the information:
            if not os.path.isfile(path):
                logger.error("Info file %s not found" % path)
                raise SModelSError()
            from smodels.experiment.expAuxiliaryFuncs import concatenateLines
            infoFile = open(self.path)
            content = concatenateLines(infoFile.readlines())
            infoFile.close()

            # Get tags in info file:
            tags = [line.split(':', 1)[0].strip() for line in content]
            modelsLine = None # the mlModels line needs to be parsed
            ## after the jsonFiles line
            for i, tag in enumerate(tags):
                if not tag:
                    continue
                if tag.startswith("#"):  # a comment!
                    continue
                line = content[i]
                value = line.split(':', 1)[1].strip()
                if tag == "mlModels":
                    modelsLine = value
                    continue
                if tag in [ "jsonFiles", "jsonFiles_FullLikelihood" ]:
                    jsonFiles = eval(value)
                    for jsonFileName,regions in jsonFiles.items():
                        newregions = self.canonizeRegions ( regions, forNN=False )
                        jsonFiles[jsonFileName] = newregions
                    value = str(jsonFiles)
                if tags.count(tag) == 1:
                    self.addInfo(tag, value)
                else:
                    logger.info(f"Ignoring unknown field {tag} found in file {self.path}" )
                    continue
            ## only now add the mlModels field
            if modelsLine != None:
                if not "'" in modelsLine and not '"' in modelsLine:
                    # did you write without qoutes?
                    modelsLine = f'"{modelsLine}"'
                mlModels = eval(modelsLine)
                if type(mlModels)==str:
                    if len(jsonFiles.values())>1:
                        logger.error ( f"mlModels {mlModels} is a single model, but we have several json files." )
                        sys.exit()
                    mlModels = { mlModels: list(jsonFiles.values())[0] }
                if type(mlModels)==dict:
                    for onnxFile,pointer in mlModels.items():
                        if type(pointer) == str:
                            pointer = jsonFiles[pointer]
                        newregions = self.canonizeRegions ( pointer, forNN=True )
                        mlModels[onnxFile]=newregions
                value = str(mlModels)
                self.addInfo("mlModels", value )

            self.cacheJsons()
            self.cacheOnnxes()

    def __eq__(self, other):
        if self.__dict__ != other.__dict__:
            return False
        return True

    def cacheOnnxes(self):
        """ if we have the "mlModels" attribute defined,
            we cache the corresponding onnx files. Needed when pickling """
        if not hasattr(self, "mlModels"):
            return
        if hasattr(self, "onnxes"):  # seems like we already have them
            return
        dirp = os.path.dirname(self.path)
        if type( self.mlModels ) in [ str ]:
            jsonFileNames = list ( self.jsonFiles.keys() )
            if len ( jsonFileNames ) == 1:
                jsonFileName = jsonFileNames[0]
                onnxFile = jsonFileName.replace(".json",".onnx")
                fullPath = os.path.join(dirp, onnxFile )
                if not os.path.exists ( fullPath ):
                    onnxFile = "model.onnx" ## fall back to standard name
                # allow shorthand notation for entries with only one json file
                self.mlModels = { onnxFile: jsonFileName }
            else:
                logger.error ( f"mlModels field in {dirp} is a string, but  {len(jsonFileNames)} json files are mentioned!" )
                sys.exit(-1)
        self.onnxes, self.smYields, self.inputMeans = {}, {}, {}
        self.inputErrors, self.nll_exp_mu0, self.nll_obs_mu0 = {}, {}, {}
        for onnxFile, jsonfilename in self.mlModels.items():
            fullPath = os.path.join(dirp, onnxFile )
            with open ( fullPath, "rb" ) as f:
                self.onnxes[onnxFile] = f.read()
                f.close()
            import onnx
            m = onnx.load ( fullPath )
            smYields, inputMeans, inputErrors = {}, [], []
            nll_exp_mu0, nll_obs_mu0 = None, None
            nllA_exp_mu0, nllA_obs_mu0 = None, None
            nll_exp_max, nll_obs_max = None, None
            nllA_exp_max, nllA_obs_max = None, None
            # smYields = []
            import json
            for em in m.metadata_props:
                if em.key == "bkg_yields":
                    st = eval(em.value)
                    for l in st: ## the sm yields are tuple of (name,value)
                        smYields[ l[0] ] = l[1]
                    #    smYields.append ( l[1] )
                if em.key == "standardization_mean":
                    inputMeans = eval(em.value)
                elif em.key == "standardization_std":
                    inputErrors = eval(em.value)
                elif em.key == 'nLL_exp_mu0':
                    nll_exp_mu0 = json.loads(em.value)
                elif em.key == 'nLL_exp_max':
                    nll_exp_max = json.loads(em.value)
                elif em.key == 'nLL_obs_max':
                    nll_obs_max = json.loads(em.value)
                elif em.key == 'nLLA_exp_max':
                    nllA_exp_max = json.loads(em.value)
                elif em.key == 'nLLA_obs_max':
                    nllA_obs_max = json.loads(em.value)
                elif em.key == 'nLL_obs_mu0':
                    nll_obs_mu0 = json.loads(em.value)
                elif em.key == 'nLLA_exp_mu0':
                    nllA_exp_mu0 = json.loads(em.value)
                elif em.key == 'nLLA_obs_mu0':
                    nllA_obs_mu0 = json.loads(em.value)
                elif em.key == 'y_min':
                    values = json.loads(em.value)
                    if True: ## not_override
                        nllA_obs_max = [None,values[-1]]
                        nllA_exp_max = [None,values[-3]]
                        nll_obs_max = [None,values[-5]]
                        nll_exp_max = [None,values[-7]]
            self.smYields[onnxFile] = smYields
            self.inputMeans[onnxFile] = inputMeans
            self.inputErrors[onnxFile] = inputErrors
            self.nll_exp_mu0[onnxFile] = nll_exp_mu0
            self.nll_obs_mu0[onnxFile] = nll_obs_mu0
            self.nllA_exp_mu0[onnxFile] = nllA_exp_mu0
            self.nllA_obs_mu0[onnxFile] = nllA_obs_mu0
            self.nll_exp_max[onnxFile] = nll_exp_max
            self.nll_obs_max[onnxFile] = nll_obs_max
            self.nllA_exp_max[onnxFile] = nllA_exp_max
            self.nllA_obs_max[onnxFile] = nllA_obs_max

    def cacheJsons(self):
        """ if we have the "jsonFiles" attribute defined,
            we cache the corresponding jsons. Needed when pickling """
        if not hasattr(self, "jsonFiles"):
            return
        if hasattr(self, "jsons"):  # seems like we already have them
            return
        import json
        self.jsons = list()
        dirp = os.path.dirname(self.path)
        jsonFiles = [os.path.join(dirp, js) for js in self.jsonFiles]
        for js in jsonFiles:
            with open(js, "r") as fi:
                try:
                    self.jsons.append(json.load(fi))
                except Exception as e:
                    logger.error(f"cannot load {js}: {e}")
                    raise(e)

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
        except SyntaxError:
            setattr(self, tag, value)
        except NameError:
            setattr(self, tag, value)
        except TypeError:
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
