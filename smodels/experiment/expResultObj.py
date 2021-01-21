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
from smodels.tools.smodelsLogging import logger
from smodels.tools.stringTools import cleanWalk
from smodels.theory.auxiliaryFunctions import getAttributesFrom,getValuesForObj

try:
    import cPickle as serializer
except ImportError as e:
    import pickle as serializer

class ExpResult(object):
    """
    Object  containing the information and data corresponding to an
    experimental result (experimental conference note or publication).
    """

    def __init__(self, path = None, discard_zeroes = True, databaseParticles = None):
        """
        :param path: Path to the experimental result folder
        :param discard_zeroes: Discard maps with only zeroes
        :param databaseParticles: the model, i.e. the particle content
        """

        if not path:
            return
        if not os.path.isdir ( path ):
            raise SModelSExperimentError ( "%s is not a path" % path )

        self.discard_zeroes = discard_zeroes
        self.path = path
        if not os.path.isfile(os.path.join(path, "globalInfo.txt")):
            logger.error("globalInfo.txt file not found in " + path)
            raise TypeError
        self.globalInfo = infoObj.Info(os.path.join(path, "globalInfo.txt"))
        #Add type of experimental result (if not defined)
        if not hasattr(self.globalInfo,'type'):
            self.globalInfo.type = 'prompt'

        datasets = {}
        folders=[]
        for root, _, files in cleanWalk(path):
            folders.append ( (root, files) )
        folders.sort()
        self.datasets = []
        hasOrder = hasattr(self.globalInfo, "datasetOrder")
        for root, files in folders:
            if 'dataInfo.txt' in files:  # data folder found
                # Build data set
                try:
                    dataset = datasetObj.DataSet(root, self.globalInfo,
                            discard_zeroes = discard_zeroes,
                            databaseParticles = databaseParticles)
                    if hasOrder:
                        datasets[dataset.dataInfo.dataId]=dataset
                    else:
                        self.datasets.append ( dataset )
                except TypeError:
                    continue
        if not hasOrder:
            return
        dsOrder = self.globalInfo.datasetOrder
        if type ( dsOrder ) == str:
            ## for debugging only, we allow a single dataset
            dsOrder = [ dsOrder ]
        for dsname in dsOrder:
            self.datasets.append ( datasets[dsname] )
        if len(self.datasets) != len(dsOrder):
            raise SModelSExperimentError ( "lengths of datasets and datasetOrder mismatch" )


    def writePickle(self, dbVersion):
        """ write the pickle file """

        meta = metaObj.Meta ( self.path, self.discard_zeroes, databaseVersion=dbVersion )
        pclfile = "%s/.%s" % ( self.path, meta.getPickleFileName() )
        logger.debug ( "writing expRes pickle file %s, mtime=%s" % (pclfile, meta.cTime() ) )
        f=open( pclfile, "wb" )
        ptcl = min ( 4, serializer.HIGHEST_PROTOCOL )
        serializer.dump(meta, f, protocol=ptcl)
        serializer.dump(self, f, protocol=ptcl)
        f.close()

    def __eq__(self, other ):
        if self.globalInfo != other.globalInfo:
            return False
        if len(self.datasets) != len ( other.datasets ):
            return False
        for (myds,otherds) in zip ( self.datasets, other.datasets ):
            if myds != otherds:
                return False
        return True

    def __ne__(self, other ):
        return not self.__eq__ ( other )

    def id(self):
        return self.globalInfo.getInfo('id')

    def __str__(self):
        label = self.globalInfo.getInfo('id') + ": "
        dataIDs = [dataset.getID() for dataset in self.datasets]
        ct_dids=0
        if dataIDs:
            for dataid in dataIDs:
                if dataid:
                    ct_dids+=1
                    label += dataid + ","
        label = "%s(%d):" % ( label[:-1], ct_dids )
        txnames = []
        for dataset in self.datasets:
            for txname in dataset.txnameList:
                tx = txname.txName
                if not tx in txnames:

                    txnames.append(tx)
        if isinstance(txnames, list):
            for txname in txnames:
                label += txname + ','
            label = "%s(%d)," % (label[:-1], len(txnames) )
        else:
            label += txnames + ','
        return label[:-1]

    def getDataset(self, dataId ):
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

    def getEfficiencyFor(self, txname, mass, dataset=None):
        """
        Convenience function. Get the efficiency for
        a specific dataset for a a specific txname.
        Equivalent to:
        self.getDataset ( dataset ).getEfficiencyFor ( txname, mass )
        """

        dataset = self.getDataset(dataset)
        if dataset:
            return dataset.getEfficiencyFor(txname, mass)
        return None

    def hasCovarianceMatrix( self ):
        return hasattr(self.globalInfo, "covariance")

    def hasJsonFile( self ):
        return hasattr(self.globalInfo, "jsonFiles")


    """ this feature is not yet ready
    def isUncorrelatedWith ( self, other ):
        can this expResult be safely assumed to be approximately uncorrelated
        with "other"? "Other" can be another expResult, or a dataset of
        an expResult.
        if self == other: return False
        if other.globalInfo.dirName ( 1 ) != self.globalInfo.dirName ( 1 ):
            return True
        # print ( "%s combinable with %s?" % ( self.globalInfo.id, other.globalInfo.id ) )
        if hasattr ( self.globalInfo, "combinableWith" ):
            #print ( "check: %s, %s" % (other.globalInfo.id, self.globalInfo.combinableWith) )
            if other.globalInfo.id in self.globalInfo.combinableWith:
                return True
        if hasattr ( other.globalInfo, "combinableWith" ):
            if self.globalInfo.id in other.globalInfo.combinableWith:
                return True
        return None ## FIXME implement
    """



    def getUpperLimitFor(self, dataID=None, alpha=0.05, expected=False,
                          txname=None, mass=None, compute=False):
        """
        Computes the 95% upper limit (UL) on the signal cross section according
        to the type of result.
        For an Efficiency Map type, returns the UL for the signal*efficiency
        for the given dataSet ID (signal region). For an Upper Limit type,
        returns the UL for the signal*BR for the given mass array and
        Txname.

        :param dataID: dataset ID (string) (only for efficiency-map type results)
        :param alpha: Can be used to change the C.L. value. The default value is 0.05
                      (= 95% C.L.) (only for  efficiency-map results)
        :param expected: Compute expected limit, i.e. Nobserved = NexpectedBG
                         (only for efficiency-map results)
        :param txname: TxName object or txname string (only for UL-type results)
        :param mass: Mass array with units (only for UL-type results)
        :param compute: If True, the upper limit will be computed
                        from expected and observed number of events.
                        If False, the value listed in the database will be used
                        instead.
        :return: upper limit (Unum object)
        """

        dataset = self.getDataset(dataID)
        if dataset:
            upperLimit = dataset.getUpperLimitFor(element=mass,expected = expected,
                                                      txnames = txname,
                                                      compute=compute,alpha=alpha)
            return upperLimit
        else:
            logger.error("Dataset ID %s not found in experimental result %s" %(dataID,self))
            return None


    def getValuesFor(self,attribute):
        """
        Returns a list for the possible values appearing in the ExpResult
        for the required attribute (sqrts,id,constraint,...).
        If there is a single value, returns the value itself.

        :param attribute: name of a field in the database (string).
        :return: list of unique values for the attribute

        """

        return getValuesForObj(self,attribute)


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


    def __lt__ ( self, other ):
        """ experimental results are sorted alphabetically according
        to their description strings """
        return str(self) < str(other)
