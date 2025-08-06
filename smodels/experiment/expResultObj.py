"""
.. module:: expResultObj
   :synopsis: Contains class that encapsulates an experimental result

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import os
from smodels.experiment import infoObj
from smodels.experiment import datasetObj
from smodels.experiment import metaObj
from smodels.experiment.exceptions import SModelSExperimentError
from smodels.base.smodelsLogging import logger
from smodels.experiment.expAuxiliaryFuncs import getAttributesFrom, getValuesForObj, cleanWalk

try:
    import cPickle as serializer
except ImportError as e:
    import pickle as serializer


class ExpResult(object):
    """
    Object  containing the information and data corresponding to an
    experimental result (experimental conference note or publication).
    """

    def __init__(self, path=None, databaseParticles=None):
        """
        :param path: Path to the experimental result folder, None means
                     transient experimental result
        :param databaseParticles: the model, i.e. the particle content
        """

        if path in [None, "<transient>"]:
            self.path = "<transient>"
            return
        if not os.path.isdir(path):
            raise SModelSExperimentError(f"{path} is not a path")

        self.path = path
        if not os.path.isfile(os.path.join(path, "globalInfo.txt")):
            logger.error("globalInfo.txt file not found in " + path)
            raise TypeError
        self.globalInfo = infoObj.Info(os.path.join(path, "globalInfo.txt"))
        # Add type of experimental result (if not defined)
        if not hasattr(self.globalInfo, 'type'):
            self.globalInfo.type = 'prompt'

        datasets = {}
        folders = []
        for root, _, files in cleanWalk(path):
            folders.append((root, files))
        folders.sort()
        self.datasets = []
        hasJsons = hasattr(self.globalInfo, "jsonFiles" )
        if hasJsons:
            dsOrder = []
            for SRs in self.globalInfo.jsonFiles.values():
                if not isinstance(SRs,list):
                    raise SModelSExperimentError(f"The entries in jsonFiles keys should be a list and not {type(SRs)}")
                for SR in SRs:
                    if not isinstance(SR,dict):
                        raise SModelSExperimentError(f"The SRs in the jsonFiles keys should be a dictionary and not {type(SR)}")
                    # In case the smodels label is not defined, ignore SR
                    if not "smodels" in SR:
                        continue
                    smodelsLabel = SR["smodels"]
                    if smodelsLabel is None or smodelsLabel in dsOrder:
                        continue
                    # Store the dataset following the order defined in globalInfo.jsonFiles
                    dsOrder.append ( smodelsLabel )
            self.globalInfo.datasetOrder = dsOrder
        hasOrder = hasattr(self.globalInfo, "datasetOrder")
        for root, files in folders:
            if 'dataInfo.txt' in files:  # data folder found
                # Build data set
                try:
                    dataset = datasetObj.DataSet(root, self.globalInfo,
                                                 databaseParticles=databaseParticles)
                    if hasOrder:
                        datasets[dataset.dataInfo.dataId] = dataset
                    else:
                        self.datasets.append(dataset)
                except TypeError as e:
                    logger.warning(f"Error creating dataset from dir {root}:\n {e}")
                    continue
        if not hasOrder:
            return
        dsOrder = self.globalInfo.datasetOrder
        if type(dsOrder) == str:
            ## for debugging only, we allow a single dataset
            dsOrder = [dsOrder]
        for dsname in dsOrder:
            self.datasets.append(datasets[dsname])
        if type(self.globalInfo.datasetOrder)==tuple:
            self.globalInfo.datasetOrder = list ( self.globalInfo.datasetOrder )
        # now append the rest -- but only for json file case
        if hasJsons:
            for dsName,ds in datasets.items():
                if dsName not in dsOrder:
                    self.datasets.append ( ds )
                    self.globalInfo.datasetOrder.append ( dsName )
        if len(self.datasets) != len(dsOrder):
            raise SModelSExperimentError( f"lengths of datasets and datasetOrder mismatch in {self.globalInfo.id}")

    def writePickle(self, dbVersion):
        """ write the pickle file """

        meta = metaObj.Meta(self.path, databaseVersion=dbVersion)
        pclfile = f"{self.path}/.{meta.getPickleFileName()}"
        logger.debug(f"writing expRes pickle file {pclfile}, mtime={meta.cTime()}")
        f = open(pclfile, "wb")
        ptcl = min(4, serializer.HIGHEST_PROTOCOL)
        serializer.dump(meta, f, protocol=ptcl)
        serializer.dump(self, f, protocol=ptcl)
        f.close()

    def __eq__(self, other):
        if self.globalInfo != other.globalInfo:
            return False
        if len(self.datasets) != len(other.datasets):
            return False
        for (myds, otherds) in zip(self.datasets, other.datasets):
            if myds != otherds:
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def id(self):
        return self.globalInfo.getInfo('id')

    def __str__(self):
        label = self.globalInfo.getInfo('id') + ": "
        dataIDs = [dataset.getID() for dataset in self.datasets]
        ct_dids = 0
        if dataIDs:
            for dataid in dataIDs:
                if dataid:
                    ct_dids += 1
                    label += dataid + ","
        label = "%s(%d):" % (label[:-1], ct_dids)
        txnames = []
        for dataset in self.datasets:
            for txname in dataset.txnameList:
                tx = txname.txName
                if not tx in txnames:

                    txnames.append(tx)
        if isinstance(txnames, list):
            for txname in txnames:
                label += txname + ','
            label = f"{label[:-1]}({len(txnames)}),"
        else:
            label += txnames + ','
        return label[:-1]

    def getDataset(self, dataId):
        """
        retrieve dataset by dataId
        """
        for dataset in self.datasets:
            if dataset.getID() == dataId:
                return dataset
        return None

    def getTxNames(self):
        """
        Returns a list of all TxName objects appearing in all datasets.
        """
        txnames = []
        for dataset in self.datasets:
            txnames += dataset.txnameList
        return txnames

    def getEfficiencyFor(self, dataID=None, txname=None, sms=None, mass=None):
        """
        For an Efficiency Map type, returns the efficiency for the corresponding
        txname and dataset for the given dataSet ID (signal region).
        For an Upper Limit type, returns 1 or 0, depending on whether the SMS
        matches the Txname.
        If SMS is not defined, but mass is given, give the efficiency using only the mass array
        (no width reweighting is applied) and the mass format is assumed
        to follow the evaluationType by the data.

        :param dataID: dataset ID (string) (only for efficiency-map type results)
        :param txname: TxName object or txname string (only for UL-type results)
        :param sms: SMS object
        :param mass: Mass array

        :return: efficiency (float)
        """


        dataset = self.getDataset(dataID)
        if dataset:
            return dataset.getEfficiencyFor(txname, sms=sms, mass=mass)
        else:
            return None

    def hasCovarianceMatrix(self):
        return hasattr(self.globalInfo, "covariance")

    def hasJsonFile(self):
        return hasattr(self.globalInfo, "jsonFiles")

    def isCombinableWith ( self, other ):
        """ can this expResult be safely assumed to be approximately uncorrelated
        with "other"? "Other" is another expResult. Later, "other" should also be
        allowed to be a dataset """
        if self == other: return False
        from smodels.base.physicsUnits import TeV
        if abs ( self.globalInfo.sqrts.asNumber(TeV) - other.globalInfo.sqrts.asNumber(TeV) ) > 1e-5:
            # different runs!
            return True
        if "CMS" in self.globalInfo.id and "ATLAS" in other.globalInfo.id:
            return True # different experiments!
        if "CMS" in other.globalInfo.id and "ATLAS" in self.globalInfo.id:
            return True # different experiments!
        ## first check the combinationsmatrix
        if hasattr ( self.globalInfo, "_combinationsmatrix" ):
            matrix = self.globalInfo._combinationsmatrix
            if self.globalInfo.id in matrix:
                if other.globalInfo.id in matrix[self.globalInfo.id]:
                    return True # found it!
            if other.globalInfo.id in matrix:
                if self.globalInfo.id in matrix[other.globalInfo.id]:
                    return True # found it!
        # now check if maybe the info is in the database itself
        if hasattr ( self.globalInfo, "combinableWith" ):
            # if we have this field in the exp result, we use it
            if other.globalInfo.id in self.globalInfo.combinableWith:
                return True
        if hasattr ( other.globalInfo, "combinableWith" ):
            if self.globalInfo.id in other.globalInfo.combinableWith:
                return True
        # nothing found? default is: False
        return False

    def getUpperLimitFor(self, dataID=None, alpha=0.05, evaluationType=False,
                         txname=None, sms=None, compute=False, mass=None):
        """
        Computes the 95% upper limit (UL) on the signal cross section according
        to the type of result.
        For an Efficiency Map type, returns the UL for the signal*efficiency
        for the given dataSet ID (signal region). For an Upper Limit type,
        returns the UL for the signal*BR for the given mass array and
        Txname.
        If SMS is not defined, but mass is given, compute the UL using only the mass array
        (no width reweighting is applied) and the mass format is assumed
        to follow the evaluationType by the data.

        :param dataID: dataset ID (string) (only for efficiency-map type results)
        :param alpha: Can be used to change the C.L. value. The default value is 0.05
                      (= 95% C.L.) (only for  efficiency-map results)
        :param expected: Compute evaluationType limit, i.e. Nobserved = NexpectedBG
                         (only for efficiency-map results)
        :param txname: TxName object or txname string (only for UL-type results)
        :param sms: SMS object
        :param mass: Mass array
        :param compute: If True, the upper limit will be computed
                        from evaluationType and observed number of events.
                        If False, the value listed in the database will be used
                        instead.
        :return: upper limit (Unum object)
        """

        dataset = self.getDataset(dataID)
        if dataset:
            upperLimit = dataset.getUpperLimitFor(sms=sms, evaluationType=evaluationType,
                                                  txnames=txname,
                                                  compute=compute, alpha=alpha,
                                                  mass=mass)
            return upperLimit
        else:
            logger.error(f"Dataset ID {dataID} not found in experimental result {self}")
            return None

    def getValuesFor(self, attribute):
        """
        Returns a list for the possible values appearing in the ExpResult
        for the required attribute (sqrts,id,constraint,...).
        If there is a single value, returns the value itself.

        :param attribute: name of a field in the database (string).
        :return: list of unique values for the attribute

        """

        return getValuesForObj(self, attribute)

    def getAttributes(self, showPrivate=False):
        """
        Checks for all the fields/attributes it contains as well as the
        attributes of its objects if they belong to smodels.experiment.

        :param showPrivate: if True, also returns the protected fields (_field)
        :return: list of field names (strings)

        """

        attributes = getAttributesFrom(self)

        if not showPrivate:
            attributes = list(filter(lambda a: a[0] != '_', attributes))

        return attributes

    def getTxnameWith(self, restrDict={}):
        """
        Returns a list of TxName objects satisfying the restrictions.
        The restrictions specified as a dictionary.

        :param restrDict: dictionary containing the fields and their allowed values.
                          E.g. {'txname' : 'T1', 'axes' : ....}
                          The dictionary values can be single entries or a list
                          of values.  For the fields not listed, all values are
                          assumed to be allowed.
        :return: list of TxName objects if more than one txname matches the selection
                 criteria or a single TxName object, if only one matches the
                 selection.

        """
        txnameList = []
        for tag, value in restrDict.items():
            for txname in self.getTxNames():
                txval = txname.getInfo(tag)
                if txval is False:
                    continue
                elif txval == value:
                    txnameList.append(txname)

        if len(txnameList) == 1:
            txnameList = txnameList[0]

        return txnameList

    def __lt__(self, other):
        """ experimental results are sorted alphabetically according
        to their description strings """
        return str(self) < str(other)
