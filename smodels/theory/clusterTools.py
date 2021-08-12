"""
.. module:: clusterTools
   :synopsis: Module holding the ElementCluster class and cluster methods used to combine similar elements according
      to the analysis.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.theory import crossSection
from smodels.theory.element import Element
from smodels.experiment.datasetObj import DataSet,CombinedDataSet
from smodels.tools.physicsUnits import fb
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.theory.auxiliaryFunctions import average
import numpy as np
from smodels.tools.smodelsLogging import logger


class AverageElement(Element):
    """
    Represents an element or list of elements containing only
    the basic attributes required for clustering or computing efficiencies/upper limits.
    Its properties are given by the average properties of the elements
    it represents and its weight is given by the total weight of
    all elements.
    """

    def __init__(self,elements=[]):
        if any(not isinstance(el,Element) for el in elements):
            raise SModelSError("An AverageElement must be created from a list of Element objects.")

        #Define relevant properties to be stored and averaged over:
        self.properties=['mass','totalwidth','txname']
        self.elements = elements[:]
        if self.elements:
            for attr in self.properties:
                setattr(self,attr,self.getAverage(attr))
            self.weight = self.elements[0].weight.copy()
            for el in self.elements[1:]:
                self.weight += el.weight

    def __str__(self):
        """
        Simply returns "averageElement", since the element
        has no well defined branches/final states (in general).

        :returns: averageElement (string)
        """

        return "averageElement"


    def __cmp__(self,other):
        """
        Compares the element with other. Only the properties
        defined in self.properties are used for comparison.
        :param other:  element to be compared (Element or AverageElement object)
        :return: -1 if self < other, 0 if self == other, +1, if self > other.
        """

        if not isinstance(other,(Element,AverageElement)):
            return -1

        otherProperties = [getattr(other,attr) for attr in self.properties]
        selfProperties = [getattr(self,attr) for attr in self.properties]
        comp = (selfProperties > otherProperties) - (otherProperties > selfProperties)

        return comp

    def __eq__(self,other):
        return self.__cmp__(other)==0

    def __lt__(self,other):
        return self.__cmp__(other)<0

    def __gt__(self,other):
        return self.__cmp__(other)>0

    def __neq__(self,other):
        return not self.__eq__(other)

    def __getattr__(self, attr):
        """
        Returns the attribute of self (necessary to overwrite the Element
        class method).

        :param attr: Attribute name

        :return: Attribute value
        """

        if not attr in self.__dict__:
            raise AttributeError ( "%s not in AverageElement" % attr )
        else:
            return self.__dict__[attr]

    def getAverage(self,attribute,weighted=True,nround=5):
        """
        Compute the average value for the attribute using
        the elements in self.elements.
        If weighted = True, compute the weighted average
        using the elements weights.

        :param attribute: Name of attribute to be averaged over (string)
        :param weighted: If True, uses the element weights to compute a weighted average
        :param nround: If greater than zero and the returning attibute is numeric, will round it
                      to this number of significant digits.

        :return: Average value of attribute.
        """

        if not self.elements:
            return None
        if len(self.elements) == 1:
            return getattr(self.elements[0],attribute)

        values = [getattr(el,attribute) for el in self.elements]
        if weighted:
            weights = [el.weight.getMaxXsec().asNumber(fb) for el in self.elements]
        else:
            weights = [1.]*len(self.elements)
        return average(values,weights,nround)

    def contains(self,element):
        """
        Check if the average element contains the element

        :param element: Element object

        :return: True/False
        """

        if any(el is element for el in self.elements):
            return True
        return False


class ElementCluster(object):
    """
    An instance of this class represents a cluster of elements.
    This class is used to store the relevant information about a cluster of
    elements and to manipulate this information.
    """

    def __init__(self, elements = [], dataset = None, distanceMatrix = None):

        self.elements = elements
        self.dataset = dataset
        self.maxInternalDist = 0.
        self._distanceMatrix = distanceMatrix
        #Compute maximal internal distance
        if self.elements and not self._distanceMatrix is None:
            self.maxInternalDist = max([self.getDistanceTo(el) for el in self])

    def __eq__(self, other):

        if type(self) != type(other):
            return False
        elif set(self.indices()) != set(other.indices()):
            return False
        else:
            return True

    def __iter__(self):
        return iter(self.elements)

    def __getitem__(self, iel):
        return self.elements[iel]

    def __len__(self):
        return len(self.elements)

    def __str__(self):
        return str(self.elements)

    def __repr__(self):
        return str(self.elements)

    def indices(self):
        """
        Return a list of element indices appearing in cluster
        """

        indices = [el._index for el in self]

        return indices

    def getTotalXSec(self):
        """
        Return the sum over the cross sections of all elements belonging to
        the cluster.

        :returns: sum of weights of all the elements in the cluster (XSectionList object)
        """
        totxsec = crossSection.XSectionList()
        for el in self:
            totxsec += el.weight
        if len(totxsec) != 1:
            logger.error("Cluster total cross section should have a single value")
            raise SModelSError()
        return totxsec[0]

    def getDataType(self):
        """
        Checks to which type of data (efficiency map or upper limit)
        the cluster refers to. It uses the cluster.dataset attribute.
        If not defined, returns None
        :return: upperLimits or efficiencyMap (string)
        """

        if self.dataset:
            dataType = self.dataset.getType()
        else:
            dataType = None

        return dataType

    def averageElement(self):
        """
        Computes the average element for the cluster.
        The average element is an empty element,
        but with its mass and width given by the average over the cluster elements.

        :return: Element object
        """

        avgEl = AverageElement(self.elements[:])
        if self.dataset:
            avgEl._upperLimit = self.dataset.getUpperLimitFor(avgEl,
                                                          txnames=avgEl.txname)

        avgEl._index = None
        return avgEl

    def copy(self):
        """
        Returns a copy of the index cluster (faster than deepcopy).
        """

        newcluster = ElementCluster(self.elements[:],self.dataset,self._distanceMatrix)
        newcluster.maxInternalDist = self.maxInternalDist
        return newcluster

    def add(self, elements):
        """
        Add an element or list of elements.

        :param elements: Element object or list of elements
        """

        if not isinstance(elements,list):
            elementList = [elements]
        else:
            elementList = elements

        for el in elementList:
            if el._index in self.indices():
                continue

            self.elements.append(el)
            #Update internal distance:
            self.maxInternalDist = max(self.maxInternalDist,self.getDistanceTo(el))

    def remove(self, elements):
        """
        Remove an element or a list of element  from the cluster.

        :param elements: Element object or list of elements
        """

        if not isinstance(elements,list):
            elementList = [elements]
        else:
            elementList = elements

        for el in elementList:
            indices = self.indices()
            if not el._index in indices:
                continue
            iel = indices.index(el._index)
            self.elements.pop(iel)
        #Update internal distance:
        self.maxInternalDist = max([self.getDistanceTo(elB) for elB in self])


    def getDistanceTo(self, element):
        """
        Return the maximum distance between any elements belonging to the
        cluster and element.

        :parameter element: Element object
        :return: maximum distance (float)
        """

        if not hasattr(element, '_upperLimit'):
            element._upperLimit = self.dataset.getUpperLimitFor(element,
                                                                txnames=element.txname)
        if element._upperLimit is None:
            return None

        #Use pre-computed distances for regular (non-averge) elements
        if not element._index is None:
            return max([self._distanceMatrix[element._index,el._index] for el in self])

        dmax = 0.
        for el in self:
            if not hasattr(el, '_upperLimit'):
                el._upperLimit = self.dataset.getUpperLimitFor(el,
                                                               txnames=el.txname)
            dmax = max(dmax,relativeDistance(element,el,self.dataset))

        return dmax

    def isConsistent(self,maxDist):
        """
        Checks if the cluster is consistent.
        Computes an average element in the cluster
        and checks if this average element belongs to the cluster
        according to the maximum allowed distance between cluster elements.

        :return: True/False if the cluster is/is not consistent.
        """

        avgElement = self.averageElement()
        if avgElement._upperLimit is None:
            return False

        dmax = self.getDistanceTo(avgElement)
        if dmax > maxDist:
            return False

        return True

def relativeDistance(el1, el2, dataset):
    """
    Defines the relative distance between two elements according to their
    upper limit values.
    The distance is defined as d = 2*|ul1-ul2|/(ul1+ul2).

    :parameter el1: Element object
    :parameter el2: Element object

    :returns: relative distance
    """

    if not hasattr(el1,'_upperLimit'):
        el1._upperLimit = dataset.getUpperLimitFor(el1,
                                                        txnames=el1.txname)
    if not hasattr(el2,'_upperLimit'):
        el2._upperLimit = dataset.getUpperLimitFor(el2,
                                                        txnames=el2.txname)

    ul1 = el1._upperLimit
    ul2 = el2._upperLimit

    if ul1 is None or ul2 is None:
        return None
    if (ul1+ul2).asNumber(fb) == 0.:
        return 0.
    ulDistance = 2.*abs(ul1 - ul2)/(ul1 + ul2)

    return ulDistance

def clusterElements(elements, maxDist, dataset):
    """
    Cluster the original elements according to their distance in upper limit space.

    :parameter elements: list of elements (Element objects)
    :parameter dataset: Dataset object to be used when computing distances in upper limit space
    :parameter maxDist: maximum distance for clustering two elements

    :returns: list of clusters (ElementCluster objects)
    """
    if len(elements) == 0:
        return []

    if any(not isinstance(el,Element) for el in elements):
        raise SModelSError("Asked to cluster non Element objects")
    if not isinstance(dataset,(DataSet,CombinedDataSet)):
        raise SModelSError("A dataset object must be defined for clustering")

    #Make sure only unique elements are clustered together (avoids double counting weights)
    #Sort element, so the ones with highest contribution (weight*eff) come first:
    elementList = sorted(elements, key = lambda el: el.weight.getMaxXsec()*el.eff, reverse=True)
    #Remove duplicated elements:
    elementsUnique = []
    for el in elementList:
        #Skip the element if it is related to any another element in the list
        if any(el.isRelatedTo(elB) for elB in elementsUnique):
            continue
        elementsUnique.append(el)

    #Get txname list only with the txnames from unique elements used for clustering
    txnames = list(set([el.txname for el in elementsUnique]))
    if dataset.getType() == 'upperLimit' and len(txnames) != 1 :
        logger.error("Clustering elements with different Txnames for an UL result.")
        raise SModelSError()


    if dataset.getType() == 'upperLimit': #Group according to upper limit values
        clusters = doCluster(elementsUnique, dataset, maxDist)
    elif dataset.getType() == 'efficiencyMap': #Group all elements together
        distanceMatrix = np.zeros((len(elementsUnique),len(elementsUnique)))
        cluster = ElementCluster(dataset=dataset,distanceMatrix=distanceMatrix)
        for iel,el in enumerate(elementsUnique):
            el._index = iel
        cluster.elements = elementsUnique
        clusters = [cluster]

    for cluster in clusters:
        cluster.txnames = txnames
    return clusters

def groupElements(elements,dataset):
    """
    Group elements into groups where the average element
    identical to all the elements in group.
    The groups contain all elements which share the same mass,width and upper limit
    and can be replaced by their average element when building clusters.

    :parameter elements: list of all elements to be grouped
    :parameter dataset: Dataset object to be used when computing distances in upper limit space

    :returns: a list of AverageElement objects
              which represents a group of elements with same mass, width and upper limit.
    """

    # First make sure all elements contain their upper limits
    for el in elements:
        if not hasattr(el,'._upperLimit'):
            el._upperLimit = dataset.getUpperLimitFor(el,txnames=el.txname)
        if el._upperLimit is None:
            raise SModelSError("Trying to cluster element outside the grid.")

    #Group elements if they have the same UL
    #and give the same average element (same mass and same width)
    avgElements = []
    for iA,elA in enumerate(elements):
        avgEl = AverageElement([elA])
        avgEl._upperLimit = elA._upperLimit
        for iB,elB in enumerate(elements):
            if iB <= iA:
                continue
            if elA._upperLimit != elB._upperLimit:
                continue
            if avgEl != elB:
                continue
            avgEl.elements.append(elB)
            avgEl.weight += elB.weight
        if not avgEl in avgElements:
            avgElements.append(avgEl)


    #Make sure each element belongs to a average element:
    for el in elements:
        nclusters = sum([avgEl.contains(el) for avgEl in avgElements])
        if nclusters != 1:
            raise SModelSError("Error computing average elements. Element %s belongs to %i average elements."
                               %(str(el),nclusters))
    return avgElements

def doCluster(elements, dataset, maxDist):
    """
    Cluster algorithm to cluster elements.

    :parameter elements: list of all elements to be clustered
    :parameter dataset: Dataset object to be used when computing distances in upper limit space
    :parameter maxDist: maximum distance for clustering two elements

    :returns: a list of ElementCluster objects containing the elements
              belonging to the cluster
    """

    #Get average elements:
    averageElements = groupElements(elements,dataset)

    #Index average elements:
    elementList = sorted(averageElements, key = lambda el: el._upperLimit)
    for iel,el in enumerate(elementList):
        el._index = iel

    #Pre-compute all necessary distances:
    distanceMatrix = np.zeros((len(elementList),len(elementList)))
    for iel,elA in enumerate(elementList):
        for jel,elB in enumerate(elementList):
            if jel <= iel:
                continue
            distanceMatrix[iel,jel] = relativeDistance(elA, elB, dataset)
    distanceMatrix = distanceMatrix + distanceMatrix.T

    #Start building maximal clusters
    clusterList = []
    for el in elementList:
        cluster = ElementCluster([],dataset,distanceMatrix)
        for elB in elementList:
            if distanceMatrix[el._index,elB._index] <= maxDist:
                cluster.add(elB)
        if not cluster.elements:
            continue
        if cluster.averageElement()._upperLimit is None:
            continue
        if not cluster in clusterList:
            clusterList.append(cluster)

    #Split the maximal clusters until all elements inside each cluster are
    #less than maxDist apart from each other and the cluster average position
    #is less than maxDist apart from all elements
    finalClusters = []
    while clusterList:
        newClusters = []
        for cluster in clusterList:
            #Check if maximal internal distance is below maxDist
            isConsistent = cluster.isConsistent(maxDist)
            if isConsistent and cluster.maxInternalDist < maxDist:
                if not cluster in finalClusters:
                    finalClusters.append(cluster)

            #Cluster violates maxDist:
            else:
                #Loop over cluster elements and if element distance
                #falls outside the cluster, remove element
                for el in cluster:
                    if cluster.getDistanceTo(el) > maxDist or not isConsistent:
                        newcluster = cluster.copy()
                        newcluster.remove(el)
                        if newcluster.averageElement()._upperLimit is None:
                            continue
                        if newcluster in newClusters:
                            continue
                        newClusters.append(newcluster)

        clusterList = newClusters
        # Check for oversized list of indexCluster (too time consuming)
        if len(clusterList) > 100:
            logger.warning("ElementCluster failed, using unclustered masses")
            finalClusters = []
            clusterList = []

    # finalClusters = finalClusters + clusterList
    # Add clusters of individual masses (just to be safe)
    for el in elementList:
        finalClusters.append(ElementCluster([el],dataset,distanceMatrix))

    # Clean up clusters (remove redundant clusters)
    for ic, clusterA in enumerate(finalClusters):
        if clusterA is None:
            continue
        for jc, clusterB in enumerate(finalClusters):
            if clusterB is None:
                continue
            if ic != jc and set(clusterB.indices()).issubset(set(clusterA.indices())):
                finalClusters[jc] = None
    while finalClusters.count(None) > 0:
        finalClusters.remove(None)

    #Replace average elements by the original elements:
    for cluster in finalClusters:
        originalElements = []
        for avgEl in cluster.elements[:]:
            originalElements += avgEl.elements[:]
        cluster.elements = originalElements[:]


    return finalClusters
