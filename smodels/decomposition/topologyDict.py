#!/usr/bin/env python3

"""
.. module:: topology
   :synopsis: Provides a Topology class and a TopologyList collection type.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>

"""

from smodels.decomposition.theorySMS import TheorySMS
from collections import OrderedDict
from typing import Any, List, Optional


class TopologyDict(OrderedDict):
    """
    An instance of this class represents an iterable collection of topologies.
    """

    def __init__(self):

        self.__dict__ = {}

    def numberOfSMS(self, canonName: Optional[int] = None) -> int:
        """
        Return the number of SMS in the dict. 
        If canonName is not None, return the number of SMS with the corresponding canonName.

        :param canonName: if None, return the total number of SMS, otherwise return only the number of SMS with the corresponding canonical name.

        :return: Number of SMS in the dict.
        """

        if canonName is not None:
            if canonName not in self:
                return 0
            else:
                return len(self[canonName])
        else:
            nSMS = 0
            for smsList in self.values():
                nSMS += len(smsList)
            return nSMS

    def addSMS(self, newSMS: TheorySMS) -> bool:
        """
        Add a new SMS to the dictionary. If an SMS with the same
        canonical name and topology already exists, merge them.

        :param newSMS: TheorySMS object to be added

        :return: True if the SMS was added, False otherwise
        """

        if isinstance(newSMS, TheorySMS):
            canonName = newSMS.canonName
            if canonName not in self:
                self[canonName] = [newSMS]
            else:
                smsList = self[canonName]
                # Find position to insert element
                # (using a bisection method)
                lo = 0
                hi = len(smsList)
                cmp = 1
                while lo < hi:
                    mid = (lo+hi)//2
                    cmp = smsList[mid].compareTo(newSMS)
                    if cmp < 0:
                        lo = mid+1
                    elif cmp > 0:
                        hi = mid
                    elif cmp == 0:
                        lo = mid
                        break
                index = lo
                if cmp == 0:
                    smsList[index] = smsList[index] + newSMS
                else:
                    smsList.insert(index, newSMS)

                self[canonName] = smsList[:]

            return True
        else:
            return False

    def getSMSList(self, canonName: Optional[int] = None) -> List[TheorySMS]:
        """
        Return a list with all the SMS appearing in the dict.
        If canonName is not None, return the SMS with the corresponding
        canonical name only.

        :param canonName: if None, return all SMS, otherwise return only the
                          the SMS with the corresponding canonical name.

        :return: List of TheorySMS objects.
        """

        allsmsList = []
        if canonName is None:
            cNames = sorted(list(self.keys()))
        else:
            cNames = [canonName]
        for cName in cNames:
            if cName not in self:
                continue
            allsmsList.extend(self[cName])
        return allsmsList

    def compress(self, doCompress: bool, doInvisible: bool, minmassgap: float, minmassgapISR: float):
        """
        Compress all SMS in the dictionary and include the compressed
        SMS in the topology list.

        :parameter doCompress: if True, perform mass compression
        :parameter doInvisible: if True, perform invisible compression
        :parameter minmassgap: value (in GeV) of the maximum
                               mass difference for compression
                               (if mass difference < minmassgap, perform mass compression)
        :parameter minmassgapISR: value (in GeV) for mass compression leading to pure
                                  ISR signature, i.e. PV > MET + MET + ... MET,
                                  (if all mass differences < minmassgapISR allow a pure ISR SMS)
        """

        for sms in self.getSMSList():
            newSMSList = sms.compress(doCompress, doInvisible, 
                                      minmassgap, minmassgapISR)
            if not newSMSList:
                continue
            for sms in newSMSList:
                self.addSMS(sms)

    def sort(self):
        """
        Sort the dictionary keys and store it in a new ordered dict in self.
        """

        newDict = OrderedDict()
        for key in sorted(self.keys()):
            newDict[key] = self.pop(key)

        for key,smsList in newDict.items():
            self[key] = smsList

    def setSMSIds(self):
        """
        Assign unique ID to each SMS in the Topology list
        """
        smsID = 1
        for sms in self.getSMSList():
            sms.smsID = smsID
            smsID += 1

    def setSMSAncestors(self):
        """
        Set the list of ancestors for each SMS in the Topology list
        keeping only the SMS which belongs to self (remove intermediate ancestors).
        """

        keepIDs = set([sms.smsID for sms in self.getSMSList()])
        for sms in self.getSMSList():
            sms.setAncestors(keepIDs)
        for sms in self.getSMSList():
            sms._ancestors = None # Clear ancestors to save memory


    def getTotalWeight(self, canonName: Optional[int] = None) -> Any:
        """
        Compute the summed cross-section over all the SMS.
        If canonName is not None, return the total cross-section
        for the SMS with the corresponding canonName.
        """

        if canonName is not None:
            smsList = self[canonName]
        else:
            smsList = self.getSMSList()

        totxsec = smsList[0].weightList
        for sms in smsList[1:]:
            totxsec = totxsec + sms.weightList

        return totxsec
