#!/usr/bin/env python3

"""
.. module:: topology
   :synopsis: Provides a Topology class and a TopologyList collection type.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>

"""

from smodels.theory.element import Element


class TopologyDict(dict):
    """
    An instance of this class represents an iterable collection of topologies.

    :ivar topos: list of topologies (Topology objects)

    """

    def __init__(self):

        self.__dict__ = {}

    def addElement(self, newelement):

        if isinstance(newelement, Element):
            canonName = newelement.canonName
            if canonName not in self:
                self[canonName] = [newelement]
            else:
                elementList = self[canonName]
                # Find position to insert element
                # (using a bisection method)
                lo = 0
                hi = len(elementList)
                while lo < hi:
                    mid = (lo+hi)//2
                    cmp = elementList[mid].compareTo(newelement)
                    if cmp < 0:
                        lo = mid+1
                    elif cmp > 0:
                        hi = mid
                    elif cmp == 0:
                        lo = mid
                        break
                index = lo
                if cmp == 0:
                    elementList[index] += newelement
                else:
                    elementList.insert(index, newelement)

                self[canonName] = elementList[:]

            return True
        else:
            return False

    def getElements(self):
        """
        Return a list with all the elements in all the topologies.

        """
        elements = []
        for elementList in self.values():
            elements.extend(elementList)
        return elements

    def compressElements(self, doCompress, doInvisible, minmassgap):
        """
        Compress all elements in the dictionary and include the compressed
        elements in the topology list.

        :parameter doCompress: if True, perform mass compression
        :parameter doInvisible: if True, perform invisible compression
        :parameter minmassgap: value (in GeV) of the maximum
                               mass difference for compression
                               (if mass difference < minmassgap, perform mass compression)

        """

        for el in self.getElements():
            newElements = el.compressElement(doCompress, doInvisible, minmassgap)
            if not newElements:
                continue
            for newelement in newElements:
                self.addElement(newelement)

    def _setElementIds(self):
        """
        Assign unique ID to each element in the Topology list
        """
        elID = 1
        for element in self.getElements():
            element.elID = elID
            elID += 1
