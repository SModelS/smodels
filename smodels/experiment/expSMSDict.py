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
from collections import OrderedDict

class ExpSMSDict(dict):
    """ 
    A two-way dictionary for storing the connections between unique ExpSMS and their
    corresponding TxNames.

    :ivar _smsDict: Dictionary mapping the unique SMS to the TxNames ({smsUnique : {TxName : smsLabel}})
    :ivar _txDict: Dictionary mapping the TxNames to the unique SMS ({TxName : {smsLabel : smsUnique}})
    :ivar _nodesDict: Dictionary mapping the node numbering in unique SMS to the
                      original numbering in the Txname SMS ({TxName : {smsLabel : nodesDict}})
    """

    def __init__(self, expResultList=[]):
        """
        :param exptResultList: List of ExptResult objects used to build the map
        """

        self._smsDict = {}
        self._txDict = {}
        self._nodesDict = {}

        if len(expResultList) > 0:
            self.computeDicts(expResultList)

    def __str__(self):
        """
        Returns description of SMS dict with number
        of unique SMS and txname objects.
        """

        strD = f"Dict with {len(self._smsDict)} unique SMS for {len(self._txDict)} txname objects"
        
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
            logger.error(f"Can not assing a {str(type(key))} object to a ExpSMSDict")
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
            logger.error(f"Can not assing a {str(type(key))} object to a ExpSMSDict")
            raise SModelSError()

        useDict[key] =  value

    def __delitem__(self, key):

        if isinstance(key,ExpSMS):
            useDict = self._smsDict
        elif isinstance(key,TxName):
            useDict = self._txDict
        else:
            logger.error(f"Can not assing a {str(type(key))} object to a ExpSMSDict")
            raise SModelSError()

        dict.__delitem__(useDict,key)

    def copy(self):
        """
        Create a copy of self.

        :return: the new ExpSMSDict object 
        """

        newDict = ExpSMSDict()
        newDict._smsDict = {key : {key2 : val2 for key2,val2 in val.items()} 
                                  for key,val in self._smsDict.items()}
        newDict._txDict = {key : {key2 : val2 for key2,val2 in val.items()} 
                                  for key,val in self._txDict.items()}
        newDict._nodesDict = {key : {key2 : val2 for key2,val2 in val.items()} 
                                  for key,val in self._nodesDict.items()}

        return newDict

    def filter(self,expResultList):
        """
        Returns a copy of self contanining only the mapping
        to the TxNames contained in expResultList.

        :param expResultList: List of experimental results (ExpResult obj)

        :return: A new ExpSMSDict with the selected TxNames.
        """

        selectedTx = []
        for exp in expResultList:
            selectedTx += exp.getTxNames()

        newDict = self.copy()
        for tx in list(newDict.getTx()):
            # If tx is in the list, do nothing
            if tx in selectedTx:
                continue
            # Get list of unique SMS matching
            # the txname:
            smsList = newDict[tx].values()
            # Remove tx from newDict:
            del newDict[tx]
            # Remove tx from the unique sms List:
            for sms in smsList:
                del newDict[sms][tx]

        # Finally, remove the unique sms
        # without any txnames attached to them
        for sms in list(newDict.getSMS()):
            if not newDict[sms]:
                del newDict[sms]

        return newDict
            
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

    def setTxNodeOrdering(self,sms,tx,smsLabel,reverse=False):
        """
        Relabel the node indices in sms accoding to the ExpSMS represented
        by the smsLabel in the TxName tx (unique ExpSMS indices - > tx sms indices).
        If reverse=True, do the the reverse labeling (tx sms indices -> unique ExpSMS indices).

        :param sms: SMS object to be relabeled
        :param tx: TxName object
        :param smsLabel: Label of the ExpSMS in the TxName object
        :param reverse: If True, do the reverse labeling

        :return: sms with indices relabeled
        """

        nodesDict = self._nodesDict[tx][smsLabel]
        # If no relabeling is needed, return original sms
        if nodesDict is None:
            return sms
        
        # Restrict the nodes dict to the indices present in the sms
        # (for the case of inclusive txnames, some nodes
        # might not be present in the sms)
        if not reverse:
            reducedNodesDict = {oldIndex : newIndex for oldIndex,newIndex in nodesDict.items()
                            if oldIndex in sms.nodeIndices}
        else:
            reducedNodesDict = {oldIndex : newIndex for newIndex,oldIndex in nodesDict.items()
                            if oldIndex in sms.nodeIndices}
        # If no relabeling is needed, return original sms
        if len(reducedNodesDict) == 0:
            return sms

        # Make sure all old nodes appear:
        # (for the case of inclusive names some nodes in
        # the sms might not appear in the dict)
        for oldIndex in sms.nodeIndices:
            if oldIndex not in reducedNodesDict:
                reducedNodesDict[oldIndex] = oldIndex
        sms.relabelNodeIndices(reducedNodesDict)
        
        return sms

    def computeDicts(self,expResultList):
        """
        Iterates over all (active) experimental results and build two dictionaries:
        one mapping TxNames and smsLabels to the unique SMS and another one
        a with the unique SMS indices as keys and a dictionary
        {TxName : smsLabel} as values. It also stores the mapping of the
        node numbering from the original Txname SMS to the unique (sorted) SMS.


        """
                
        # Loop over active experimental results:
        # (Note that expResultList returns only the active results)
        for exp in expResultList:
            # Loop over txnames:
            for tx in exp.getTxNames():
                # Loop over ExpSMS in the txname:
                if tx in self.getTx():
                    logger.error(f"TxName {tx} for {exp} already present in the map (should never happen)!")
                    raise SModelSError()
                self[tx] = {}
                self._nodesDict[tx] = {}
                for sms,smsLabel in tx.smsMap.items():
                    newSMS = sms.copy()
                    # Keep the mapping between original nodes
                    # and the sorted nodes (if the SMS has changed)
                    nodesDict = newSMS.sort()
                    if any(key != val for key,val in nodesDict.items()):
                        invertNodesDict = {newIndex : oldIndex 
                                          for oldIndex,newIndex 
                                          in nodesDict.items()}
                    else:
                        invertNodesDict = None
                    smsMatch = None
                    # Loop for an identical SMS in smsDict
                    # (there can only be one match, since the SMS within
                    # a given txname must be unique)
                    for sms1 in self.getSMS():
                        if newSMS.identicalTo(sms1):
                            smsMatch = sms1
                            break
                                        
                    # Update dictionary
                    if smsMatch is None:
                        self[newSMS] = {tx : smsLabel}
                        self[tx][smsLabel] = newSMS
                        self._nodesDict[tx][smsLabel] = invertNodesDict
                    else:
                        self[smsMatch].update({tx : smsLabel})
                        self[tx][smsLabel] = smsMatch
                        self._nodesDict[tx][smsLabel] = invertNodesDict

        # Sort txnames by dataMap 
        # (the one requiring the largest amount of data comes first)
        for txsms in self.getSMS():
            txDict = self[txsms]
            sortedTx = sorted(txDict.keys(), key = lambda tx: len(tx.dataMap), reverse=True)
            newDict = OrderedDict([(tx,txDict[tx]) for tx in sortedTx])
            self[txsms] = newDict


    def getMatchesFrom(self,smsTopDict):
        """
        Checks for all the matches between the SMS in smsTopDict (from decomposition)
        and the unique SMS in self. Returns a dictionary with the mapping:
        {unique SMS : [(matched SMS, original SMS),...]}

        :param smsTopDict: TopologyDict object with the TheorySMS from decomposition
        :return: Dictionary with unique SMS as keys and lists of matched SMS as values.
        """


        # Create dictionary with matches
        smsMatch = {expSMS : [] for expSMS in self.getSMS()}
        # Group expSMS by canon names:
        cNamesDict = {}
        for expSMS in smsMatch:
            if expSMS.canonName not in cNamesDict:
                cNamesDict[expSMS.canonName] = []
            cNamesDict[expSMS.canonName].append(expSMS)

        # Loop over theory SMS (from decomposition)
        # and compute matches:
        for sms in smsTopDict.getSMSList():
            canonName  = sms.canonName
            # Select txSMS with matching canon name:
            selectedSMSlist = []
            for cName,expSMSList in cNamesDict.items():
                if cName == canonName:
                    selectedSMSlist += expSMSList

            # Loop over selected SMS and check for matches
            # Store the (correctly ordered) match in smsMatch
            for expSMS in selectedSMSlist:
                # Use the first txname to compute the match:
                txname = list(self[expSMS].keys())[0]
                smsLabel = self[expSMS][txname]
                matchedSMS = txname.hasSMSas(sms,useLabel=smsLabel)
                if matchedSMS is None:
                    continue
                # Now set the node indices according to the
                # unique ExpSMS labels:
                self.setTxNodeOrdering(matchedSMS,txname,
                                        smsLabel,reverse=True)
                matchedSMS.smsID = sms.smsID
                smsMatch[expSMS].append((matchedSMS,sms))    

        return smsMatch