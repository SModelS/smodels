#!/usr/bin/env python3

"""
.. module:: topology
   :synopsis: Provides a Topology class and a TopologyList collection type.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>

"""

from smodels.theory.theorySMS import TheorySMS
from collections import OrderedDict


class TopologyDict(dict):
    """
    An instance of this class represents an iterable collection of topologies.
    """

    def __init__(self):

        self.__dict__ = OrderedDict()

    def addSMS(self, newSMS):

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

    def getSMSList(self):
        """
        Return a list with all the SMS appearing in the dict

        """
        allsmsList = []
        for smsList in self.values():
            allsmsList.extend(smsList)
        return allsmsList

    def compress(self, doCompress, doInvisible, minmassgap):
        """
        Compress all SMS in the dictionary and include the compressed
        SMS in the topology list.

        :parameter doCompress: if True, perform mass compression
        :parameter doInvisible: if True, perform invisible compression
        :parameter minmassgap: value (in GeV) of the maximum
                               mass difference for compression
                               (if mass difference < minmassgap, perform mass compression)

        """

        for sms in self.getSMSList():
            newSMSList = sms.compressElement(doCompress, doInvisible, minmassgap)
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
            newDict[key] = self[key]

        self.__dict__ = newDict

    def setSMSIds(self):
        """
        Assign unique ID to each SMS in the Topology list
        """
        smsID = 1
        for sms in self.getSMSList():
            sms.smsID = smsID
            smsID += 1

