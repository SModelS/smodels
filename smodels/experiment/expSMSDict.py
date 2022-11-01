#!/usr/bin/env python3

"""
.. module:: expSMSDict
   :synopsis: A two-way dictionary containing the links between the (TxName, smsLabel) and the unique ExpSMS
              objects. It can be used to quickly map the the TxName to the list of unique SMS and vice-versa.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.experiment.txnameObj import TxName
from smodels.experiment.expSMS import ExpSMS
from smodels.base.smodelsLogging import logger
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError

class ExpSMSDict(dict):
    """ 
    A two-way dictionary for storing the connections between unique ExpSMS and their
    corresponding TxNames.
    """

    def __init__(self, expResultList):
        """
        :param exptResultList: List of ExptResult objects used to build the map
        """

        self._smsDict = {}
        self._txDict = {}

        if len(expResultList) > 0:
            self.computeDicts(expResultList)

    def __str__(self):
        """
        Returns description of SMS dict with number
        of unique SMS and txname objects.
        """

        strD = "Dict with %i unique SMS for %i txname objects" %(len(self._smsDict),len(self._txDict))
        
        return strD
    
    def __repr__(self):

        return str(self)

    def __len__(self):
        """
        Returns the number of unique SMS
        """

        return len(self._smsDict)

    def __iter__(self):
        return iter(self._smsDict)

    def __getitem__(self,key):

        if isinstance(key,ExpSMS):
            useDict = self._smsDict
        elif isinstance(key,TxName):
            useDict = self._txDict
        else:
            logger.error("Can not assing a %s object to a ExpSMSDict" %str(type(key)))
            raise SModelSError()

        return useDict[key]

    def __setitem__(self, key, value):
        """
        Add item to self._smsDict if key is an ExpSMS
        or to self._txDict if key in a TxNameObj
        """

        if isinstance(key,ExpSMS):
            useDict = self._smsDict
        elif isinstance(key,TxName):
            useDict = self._txDict
        else:
            logger.error("Can not assing a %s object to a ExpSMSDict" %str(type(key)))
            raise SModelSError()

        useDict[key] =  value

    def __delitem__(self, key):

        if isinstance(key,ExpSMS):
            useDict = self._smsDict
        elif isinstance(key,TxName):
            useDict = self._txDict
        else:
            logger.error("Can not assing a %s object to a ExpSMSDict" %str(type(key)))
            raise SModelSError()

        dict.__delitem__(useDict,key)

    def getSMS(self):
        """
        Iterate over the unique ExpSMS stored in self.
        """

        return iter(self._smsDict)

    def getTx(self):
        """
        Iterate over the TxNames stored in self.
        """

        return iter(self._txDict)


    def computeDicts(self,expResultList):
        """
        Iterates over all (active) experimental results and build two dictionaries:
        one mapping indices to the unique SMS and another one
        a with the unique SMS indices as keys and a dictionary
        {ExpResult_index : {DataSet_index : {TxName_index : smsLabel}}} as values.

        :return: Index->SMS dictionary and SMS_index -> Exp_index dictionary
        """
                
        # Loop over active experimental results:
        # (Note that expResultList returns only the active results)
        for exp in expResultList:
            # Loop over txnames:
            for tx in exp.getTxNames():
                # Sort TxName (if it is already tagged as sorted, do nothing)
                tx.sortSMSMap()
                # Loop over ExpSMS in the txname:
                if tx in self.getTx():
                    logger.error("TxName %s fpr %s already present in the map (should never happen)!" %(tx,exp))
                    raise SModelSError()
                self[tx] = {}
                for sms,smsLabel in tx.smsMap.items():
                    smsMatch = None
                    # Loop for an identical SMS in smsDict
                    # (there can only be one match, since the SMS within
                    # a given txname must be unique)
                    for sms1 in self.getSMS():
                        if sms.identicalTo(sms1):
                            smsMatch = sms1
                            break
                                        
                    # Update dictionary
                    if smsMatch is None:
                        self[sms] = {tx : smsLabel}
                        self[tx][smsLabel] = sms
                    else:
                        self[smsMatch].update({tx : smsLabel})
                        self[tx][smsLabel] = smsMatch
