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
    (luminosity, sqrts, experimentID,...). Its attributes are generated according to the lines in the
    .txt file which contain "info_tag: value".
    """

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
            for i, tag in enumerate(tags):
                if not tag:
                    continue
                if tag.startswith("# "):  # a comment!
                    continue
                line = content[i]
                value = line.split(':', 1)[1].strip()
                if tag in [ "jsonFiles", "jsonFiles_FullLikelihood" ]:
                    jsonFiles = eval(value)
                    for jsonFileName,regions in jsonFiles.items():
                        newregions = []
                        for region in regions:
                            if type(region)==str:
                                region={"smodels": region}
                            if not "type" in region:
                                region["type"]="SR"
                            if not "smodels" in region:
                                region["smodels"]=None
                            newregions.append ( region )
                        jsonFiles[jsonFileName] = newregions
                    value = str(jsonFiles)
                if tags.count(tag) == 1:
                    self.addInfo(tag, value)
                else:
                    logger.info("Ignoring unknown field %s found in file %s"
                                % (tag, self.path))
                    continue

            self.cacheJsons()
            self.cacheOnnx()

    def __eq__(self, other):
        if self.__dict__ != other.__dict__:
            return False
        return True

    def cacheOnnx(self):
        """ if we have the "mlModel" attribute defined,
            we cache the corresponding onnx. Needed when pickling """
        if not hasattr(self, "mlModel"):
            return
        if hasattr(self, "onnx"):  # seems like we already have them
            return
        dirp = os.path.dirname(self.path)
        mlModel = os.path.join(dirp, self.mlModel)
        with open ( mlModel, "rb" ) as f:
            self.onnx = f.read()
            f.close()
        import onnx
        m = onnx.load ( mlModel )
        smYields, inputMeans, inputErrors = {}, [], []
        nll_exp_mu0, nll_obs_mu0 = None, None
        # smYields = []
        import json
        for em in m.metadata_props:
            if em.key == "bkg_yields":
                st = eval(em.value)
                for l in st:
                    smYields[ l[0] ] = l[1]
                #    smYields.append ( l[1] )
            if em.key == "standardization_mean":
                inputMeans = eval(em.value)
            elif em.key == "standardization_std":
                inputErrors = eval(em.value)
            elif em.key == 'nLL_exp_mu0':
                nll_exp_mu0 = json.loads(em.value)
            elif em.key == 'nLL_obs_mu0':
                nll_obs_mu0 = json.loads(em.value)
        self.smYields = smYields
        self.inputMeans = inputMeans
        self.inputErrors = inputErrors
        self.nll_exp_mu0 = nll_exp_mu0
        self.nll_obs_mu0 = nll_obs_mu0

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
