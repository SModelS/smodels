#!/usr/bin/env python3

"""
.. module:: topology
   :synopsis: Provides a Topology class and a TopologyList collection type.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>

"""

from smodels.theory import crossSection
from smodels.theory.element import Element
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.theory.auxiliaryFunctions import index_bisect
from smodels.tools.smodelsLogging import logger
import itertools


class Topology(object):
    """
    An instance of this class represents a topology.
    """

    def __init__(self, elements=None):
        """
        Constructor.
        If elements is defined, create the topology from it. If elements it is
        a list, all elements must share a common global topology.

        :parameter elements: Element object or list of Element objects
        """

        self.vertnumb = []
        self.vertparts = []
        self.elementList = []

        if elements:
            if isinstance(elements,Element):
                self.addElement(elements)
            elif isinstance(elements,list):
                for element in elements:
                    self.addElement(element)

    def __str__(self):
        """
        Return string with numbers of particles per vertex, e.g.
        [1],[2,1]

        :returns: string with number of final states in each branch
        """
        ret = ""
        for p in self.vertparts:
            ret += "%s" % str(p).replace(" ", "")
        return ret

    def __repr__(self):
        return self.__str__()

    def __ne__(self,other):
        return not ( self.__eq__(other) )

    def __eq__(self,other):
        ret = (self.__cmp__(other)==0 )
        return ret

    def __lt__(self,other):
        return self.__cmp__(other)<0

    def __cmp__(self,other):
        """
        Compares the topology with other.
        The comparison is made on number of vertices and then on the
        total number of particles coming out of the vertices.
        :param other:  topology to be compared (Topology object)
        :return: -1 if self < other, 0 if self == other, +1, if self > other.
        """

        #Check for any permutation of branches:
        for v1 in itertools.permutations(self.vertparts):
            v1 = list(v1)
            if v1 == other.vertparts:
                return 0

        if sorted(self.vertnumb,reverse=True) != sorted(other.vertnumb,reverse=True):
            comp = sorted(self.vertnumb,reverse=True) > sorted(other.vertnumb,reverse=True)
            if comp: return 1
            else: return -1
        elif sorted(self.vertparts) != sorted(other.vertparts):
            comp = sorted(self.vertparts) > sorted(other.vertparts)
            if comp: return 1
            else: return -1
        else:
            return 0


    def checkConsistency(self):
        """
        Check if the all the elements in elementList are
        consistent with the topology (same number of vertices and final states)

        :returns: True if all the elements are consistent. Print error message
                  and exits otherwise.
        """

        for element in self.elementList:
            info = element.getEinfo()
            if self.vertnumb != info["vertnumb"]:
                logger.error("Inconsistent topology.")
                raise SModelSError()
            if self.vertparts != info["vertparts"]:
                logger.error("Inconsistent topology.")
                raise SModelSError()
        logger.info("Consistent topology.")
        return True


    def describe(self):
        """
        Create a detailed description of the topology.

        :returns: list of strings with a description of the topology
        """
        ret = ("number of vertices: %s, number of vertex particles: %s, "
               "number of elements: %d" % \
               (self.vertnumb, self.vertparts, len(self.elementList)))
        return ret


    def getElements(self):
        """
        Get list of elements of the topology.

        :return: elementList (list of Element objects)
        """
        return self.elementList


    def addElement(self, newelement):
        """
        Add an Element object to the elementList.

        For all the pre-existing elements, which match the new element, add
        weight. If no pre-existing elements match the new one, add it to the
        list. OBS: newelement MUST ALREADY BE SORTED (see element.sort())

        :parameter newelement: element to be added (Element object)
        :returns: True, if the element was added. False, otherwise
        """

        # If the topology info has not been set yet, set it using the element
        # info

        if not self.vertparts:
            self.vertparts = newelement.getEinfo()["vertparts"]
        if not self.vertnumb:
            self.vertnumb = newelement.getEinfo()["vertnumb"]

        #First check if element matches topology structure
        info = newelement.getEinfo()

        if info != self._getTinfo():
            logger.warning('Element to be added does not match topology')
            return False

        index = index_bisect(self.elementList,newelement)
        if index != len(self.elementList) and self.elementList[index] == newelement:
            self.elementList[index] += newelement
        else:
            self.elementList.insert(index,newelement)

        return True



    def _getTinfo(self):
        """
        Return a dictionary with the topology number of vertices and vertparts.

        :returns: dictionary with topology information
        """
        return {'vertnumb' : self.vertnumb, 'vertparts' : self.vertparts}


    def getTotalWeight(self):
        """
        Return the sum of all elements weights.

        :returns: sum of weights of all elements (XSection object)
        """
        if len(self.elementList) == 0:
            return None

        sumw = crossSection.XSectionList()
        for element in self.elementList:
            sumw += element.weight

        return sumw


class TopologyList(object):
    """
    An instance of this class represents an iterable collection of topologies.
    """
    
    def __init__(self, topologies=[]):
        """
        Add topologies sequentially, if provided.
        """

        self.topos = []
        for topo in topologies:
            self.add(topo)

    def __ne__(self,other):
        return not self.__eq__(other)

    def __eq__(self,other):
        return self.topos == other.topos

    def __len__(self):
        return len(self.topos)


    def __getitem__(self, index):
        return self.topos[index]


    def __iter__(self):
        return iter(self.topos)


    def __str__(self):
        s = "TopologyList:\n"
        for topo in self.topos:
            s += str(topo) + "\n"
        return s

    def __repr__(self):
        return self.__str__()

    def insert(self,index,topo):
        self.topos.insert(index,topo)


    def addList(self, topoList):
        """
        Adds topologies in topoList using the add method.

        """
        for topo in topoList:
            self.add(topo)


    def describe(self):
        """
        Returns string with basic information about the topology list.

        """
        s = "TopologyList:\n"
        for topo in self.topos:
            s += str(topo) + "\n"
        return s

    def index(self,topo):
        """
        Uses bisect to find the index where of topo in the list.
        If topo does not appear in the list, returns None.

        :param topo: Topology object
        :return: position of topo in the list. If topo does not
                appear in the list, return None.
        """

        i = index_bisect(self, topo)
        if i != len(self) and self[i] == topo:
            return i

        return None


    def hasTopology(self,topo):
        """
        Checks if topo appears in any of the topologies in the list.

        :param topo: Topology object
        :return: True if topo appears in the list, False otherwise.
        """

        for t in self:
            if t == topo:
                return True

        return False


    def add(self, newTopology):
        """
        Check if elements in newTopology matches an entry in self.topos.

        If it does, add weight. If the same topology exists, but not the same
        element, add element. If neither element nor topology exist, add the
        new topology and all its elements.

        :param newTopology: Topology object

        """

        index = index_bisect(self, newTopology)
        if index != len(self) and self[index] == newTopology:
            for newelement in newTopology.elementList:
                self.topos[index].addElement(newelement)
        else:
            self.insert(index,newTopology)


    def addElement(self, newelement):
        """
        Add an Element object to the corresponding topology in the list.
        If the element topology does not match any of the topologies in
        the list, create a new topology and insert it in the list.
        If the element topology already exists, add it to the respective
        topology.
        :parameter newelement: element to be added (Element object)
        :returns: True, if the element was added. False, otherwise
        """

        #First create a dummy topology from the element to check
        #if this topology already exists in the list:
        elInfo = newelement.getEinfo()
        topoDummy = Topology()
        topoDummy.elementList.append(newelement)
        topoDummy.vertnumb = elInfo["vertnumb"]
        topoDummy.vertparts = elInfo["vertparts"]

        index = index_bisect(self,topoDummy)
        if index != len(self) and self.topos[index] == topoDummy:
            self.topos[index].addElement(newelement)
        else:
            self.topos.insert(index,topoDummy)


    def getTotalWeight(self):
        """
        Return the sum of all topologies total weights.

        """
        sumw = crossSection.XSectionList()
        for topo in self:
            topoweight = topo.getTotalWeight()
            if topoweight:
                sumw += topoweight
        return sumw


    def getElements(self):
        """
        Return a list with all the elements in all the topologies.

        """
        elements = []
        for top in self.topos:
            elements.extend(top.elementList)
        return elements

    def compressElements(self,doCompress,doInvisible,minmassgap):
        """
        Compress all elements in the list and include the compressed
        elements in the topology list.

        :parameter doCompress: if True, perform mass compression
        :parameter doInvisible: if True, perform invisible compression
        :parameter minmassgap: value (in GeV) of the maximum
                               mass difference for compression
                               (if mass difference < minmassgap, perform mass compression)

        """

        for el in self.getElements():
            newElements = el.compressElement(doCompress,doInvisible,minmassgap)
            if not newElements:
                continue
            for newelement in newElements:
                newelement.sortBranches()  #Make sure elements are sorted BEFORE adding them
                self.addElement(newelement)

    def _setElementIds(self):
        """
        Assign unique ID to each element in the Topology list
        """
        elID = 1
        for element in self.getElements():
            element.elID = elID
            elID += 1
