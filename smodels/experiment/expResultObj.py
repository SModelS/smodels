"""
.. module:: expResultObj
   :synopsis: Contains class that encapsulates an experimental result

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import os
from smodels.experiment import infoObj
from smodels.experiment import txnameObj
from smodels.experiment import datasetObj
from smodels.experiment.exceptions import DatabaseNotFoundException
from smodels.tools.physicsUnits import fb

from smodels.tools.smodelsLogging import logger

class ExpResult(object):
    """
    Object  containing the information and data corresponding to an
    experimental result (experimental conference note or publication).
    
    :ivar path: path to the experimental result folder (i.e. ATLAS-CONF-2013-047)
    :ivar globalInfo: Info object holding the data in <path>/globalInfo.txt
    :ivar datasets: List of DataSet objects corresponding to the dataset folders 
                    in <path>    
    """
        
    def __init__(self, path=None, discard_zeroes = True ):
        """
        :param path: Path to the experimental result folder
        :param discard_zeroes: Discard maps with only zeroes
        """ 

        if path and os.path.isdir(path):
            self.path = path
            if not os.path.isfile(os.path.join(path, "globalInfo.txt")):
                logger.error("globalInfo.txt file not found in " + path)
                raise TypeError
            self.globalInfo = infoObj.Info(os.path.join(path, "globalInfo.txt"))
            self.datasets = []
            folders=[]
            for root, _, files in os.walk(path):
                folders.append ( (root, files) )
            folders.sort()
            for root, files in folders:
                if 'dataInfo.txt' in files:  # data folder found
                    # Build data set
                    try:
                        dataset = datasetObj.DataSet(root, self.globalInfo,
                                discard_zeroes = discard_zeroes )
                        self.datasets.append(dataset)
                    except TypeError:
                        continue

    def __eq__(self, other ):
        if self.globalInfo != other.globalInfo:
            return False
        if len(self.datasets) != len ( other.datasets ):
            return False
        for (myds,otherds) in zip ( self.datasets, other.datasets ):
            if myds != otherds:
                return False
        return True

    def __lt__ ( self, other ):
        """ experimental results are sorted alphabetically according
            to their description strings """
        return str(self) < str(other)


    def __ne__(self, other ):
        return not self.__eq__ ( other )

    def __str__(self):
        label = self.globalInfo.getInfo('id') + ": "
        dataIDs = [dataset.dataInfo.dataId for dataset in self.datasets]
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
            if dataset.dataInfo.dataId == dataId:
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

    def getEfficiencyFor ( self, txname, mass, dataset=None):
        """
        Convenience function. Get the efficiency for
        a specific dataset for a a specific txname.
        Equivalent to:
        self.getDataset ( dataset ).getEfficiencyFor ( txname, mass )
        """
        dataset = self.getDataset ( dataset )
        if dataset: return dataset.getEfficiencyFor ( txname, mass )
        return None

    def getUpperLimitFor(self, dataID=None, alpha=0.05, expected=False,
                          txname=None, mass=None, compute=False):
        """
        Computes the 95% upper limit (UL) on the signal cross section according
        to the type of result.
        For an Efficiency Map type, returns the UL for the signal*efficiency
        for the given dataSet ID (signal region). For an Upper Limit type,
        returns the UL for the signal*BR for for the given mass array and
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
        if self.datasets[0].dataInfo.dataType == 'efficiencyMap':
            if dataID == None:
                """ FIXME when no dataID is given, and we have a covariance
                    matrix, compute the combined CL UL """
                print ( "globalInfo covariance matrix? %s" % self.globalInfo.covariance )
                if hasCovarianceMatrix:
                    return combinedUL
            if not dataID or not isinstance(dataID, str):
                logger.error("The data set ID must be defined when computing ULs" \
                             " for efficiency-map results (as there is no covariance matrix).")
                return False
            
            useDataset = False
            for dataset in self.datasets:
                if dataset.dataInfo.dataId == dataID:
                    useDataset = dataset
                    break
            if useDataset is False:
                logger.error("The data set ID not found.")
                return False
                
            if compute:
                upperLimit = useDataset.getSRUpperLimit(alpha, expected, compute)
            else:
                if expected:
                    upperLimit = useDataset.dataInfo.expectedUpperLimit
                else:
                    upperLimit = useDataset.dataInfo.upperLimit
                if (upperLimit/fb).normalize()._unit:
                    logger.error("Upper limit defined with wrong units for %s and %s"
                                  %(dataset.globalInfo.id,dataset.dataInfo.dataId))
                    return False           
            
        elif self.datasets[0].dataInfo.dataType == 'upperLimit':
            if not txname or not mass:
                logger.error("A TxName and mass array must be defined when \
                             computing ULs for upper-limit results.")
                return False
            if not isinstance(txname, txnameObj.TxName) and \
            not isinstance(txname, str):
                logger.error("txname must be a TxName object or a string")
                return False
            if not isinstance(mass, list):
                logger.error("mass must be a mass array")
                return False
            for tx in self.getTxNames():
                if tx == txname or tx.txName == txname:
                    if expected:
                        if not tx.txnameDataExp: upperLimit = None
                        else: upperLimit = tx.txnameDataExp.getValueFor(mass)
                    else: upperLimit = tx.txnameData.getValueFor(mass)
        else:
            logger.warning("Unkown data type: %s. Data will be ignored.", 
                           self.datasets[0].dataInfo.dataType)

        
        return upperLimit
    

    def getValuesFor(self, attribute=None):
        """
        Returns a list for the possible values appearing in the ExpResult
        for the required attribute (sqrts,id,constraint,...).
        If there is a single value, returns the value itself.
        
        :param attribute: name of a field in the database (string). If not
                          defined it will return a dictionary with all fields 
                          and their respective values
        :return: list of values or value
        
        """
        fieldDict = list ( self.__dict__.items() )
        valuesDict = {}
        while fieldDict:
            for field, value in fieldDict:
                if not '<smodels.experiment' in str(value):
                    if not field in valuesDict:
                        valuesDict[field] = [value]
                    else: valuesDict[field].append(value)
                else:
                    if isinstance(value, list):
                        for entry in value:
                            fieldDict += entry.__dict__.items()
                    else: fieldDict += value.__dict__.items()
                fieldDict.remove((field, value))

        # Try to keep only the set of unique values
        for key, val in valuesDict.items():
            try:
                valuesDict[key] = list(set(val))
            except TypeError:
                pass
        if not attribute:
            return valuesDict
        elif not attribute in valuesDict.keys():
            logger.warning("Could not find field %s in %s", attribute, self.path)
            return False
        else:
            return valuesDict[attribute]


    def getAttributes(self, showPrivate=False):
        """
        Checks for all the fields/attributes it contains as well as the
        attributes of its objects if they belong to smodels.experiment.
        
        :param showPrivate: if True, also returns the protected fields (_field)
        :return: list of field names (strings)
        
        """
        fields = self.getValuesFor().keys()
        fields = list(set(fields))

        if not showPrivate:
            for field in fields[:]:
                if "_" == field[0]:
                    fields.remove(field)
        return fields
    

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
