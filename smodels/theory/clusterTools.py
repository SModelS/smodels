"""
.. module:: clusterTools
   :synopsis: Module holding the ElementCluster class and cluster methods used to combine similar elements according
      to the analysis.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

from smodels.theory import crossSection
from smodels.theory.element import Element
from smodels.experiment.datasetObj import DataSet,CombinedDataSet
from smodels.tools.physicsUnits import fb, GeV
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
import numpy as np
from smodels.tools.smodelsLogging import logger

class ElementCluster(object):
    """
    An instance of this class represents a cluster of elements.
    This class is used to store the relevant information about a cluster of
    elements and to manipulate this information.
    
    :ivar elements: list of elements in the cluster (Element objects)    
    """

    def __init__(self, elements = [], maxDist = 0., dataset = None):

        self.elements = elements
        self.maxDist = maxDist
        self.dataset = dataset

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
    
    def relativeDistance(self,el1, el2):
        """
        Defines the relative distance between two elements according to their
        upper limit values.
        The distance is defined as d = 2*|ul1-ul2|/(ul1+ul2).
       
        :parameter el1: Element object
        :parameter el2: Element object
    
        :returns: relative distance
        """
    
        if not hasattr(el1,'_upperLimit'):
            el1._upperLimit = self.dataset.getUpperLimitFor(el1,
                                                            txnames=el1.txname)
        if not hasattr(el2,'_upperLimit'):
            el2._upperLimit = self.dataset.getUpperLimitFor(el2,
                                                            txnames=el1.txname)

        ul1 = el1._upperLimit
        ul2 = el2._upperLimit    
        
        if ul1 is None or ul2 is None:
            return None
        ulDistance = 2.*abs(ul1 - ul2)/(ul1 + ul2)
    
        return ulDistance
    

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

    def getAvgMass(self):
        """
        Return the average mass of all elements belonging to the cluster, weighted
        by their individual weights.
        If the masses have distinct shapes (e.g. for different txnames in efficiencyMap
        results), return None.
        
        :returns: average mass array appearing in the cluster     
        """                          
        
        massList = [el.mass for el in self]
        if any(np.array(m).shape != np.array(massList[0]).shape for m in massList):
            return None
        
        weights = [el.weight.getMaxXsec().asNumber(fb) for el in self]
        avgmass = [[0. for m in br] for br in massList[0]]
        for ib, branch in enumerate(massList[0]):
            for im,_ in enumerate(branch):
                vals = [mass[ib][im].asNumber(GeV) for mass in massList]
                weights = [ float(weight) for weight in weights ]
                avg = np.average(vals,weights=weights)
                avgmass[ib][im] = avg*GeV

        return avgmass

    def getAvgWidth(self):
        """
        Return the average width of all elements belonging to the cluster, weighted
        by their individual weights.
        If the widths have distinct shapes (e.g. for different txnames in efficiencyMap
        results), return None.

        :returns: average width array appearing in the cluster
        """                          

        widthList = [el.totalwidth for el in self]
        if any(np.array(w).shape != np.array(widthList[0]).shape for w in widthList):
            return None

        weights = [el.weight.getMaxXsec().asNumber(fb) for el in self]
        avgwidth = [[0. for w in br] for br in widthList[0]]
        for ib, branch in enumerate(widthList[0]):
            for im,_ in enumerate(branch):
                vals = [mass[ib][im].asNumber(GeV) for mass in widthList]
                weights = [ float(weight) for weight in weights ]
                avg = np.average(vals,weights=weights)
                avgwidth[ib][im] = avg*GeV

        return avgwidth
    
    def averageElement(self):
        """
        Computes the average element for the cluster.
        The average element is an empty element,
        but with its mass and width given by the average over the cluster elements.

        :return: Element object
        """

        avgEl = self.elements[0].__class__()
        avgEl.mass = self.getAvgMass()
        avgEl.totalwidth = self.getAvgWidth()
        avgEl.txname = self.elements[0].txname
        avgEl._upperLimit = self.dataset.getUpperLimitFor(avgEl,
                                                          txnames=avgEl.txname)
        return avgEl

    def getPIDs(self):
        """
        Return the list of all PIDs appearing in all elements in the cluster,
        i.e. [ [[pdg1, pdg2,...],[pdg3,pdg4,...]], [[pdg1', pdg2',...],[pdg3',pdg4',...]]

        :returns: list of PIDs
        """
        
        PIDs = []
        for el in self:
            for pidList in el.getMothers():
                if not pidList in PIDs:
                    PIDs.append(pidList)
            
        return PIDs
    
    def getIDs(self):
        """
        Return list of all element IDs appearing in the cluster
        :return: list of element IDs
        """
        IDs = []
        for el in self:
            if not el.elID in IDs: IDs.append(el.elID)
        return IDs

    def copy(self):
        """
        Returns a copy of the index cluster (faster than deepcopy).
        """

        newcluster = ElementCluster(self.elements,self.maxDist,self.dataset)
        newcluster.maxDist = self.maxDist
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

        dmax = 0.
        for el in self:
            if not hasattr(el, '_upperLimit'):
                el._upperLimit = self.dataset.getUpperLimitFor(el,
                                                               txnames=el.txname)
            dmax = max(dmax,self.relativeDistance(element,el))

        return dmax

    def getMaxInternalDist(self):
        """
        Return the maximum distance between any pair of elements belonging
        to the cluster.
        """

        dmax = max([self.getDistanceTo(el) for el in self])

        return dmax

    def isConsistent(self):
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
        if dmax > self.maxDist:
            return False

        return True



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

    txnames = list(set([el.txname for el in elements]))
    if dataset.getType() == 'upperLimit' and len(txnames) != 1 :
        logger.error("Clustering elements with different Txnames for an UL result.")
        raise SModelSError()

    # ElementCluster elements:
    clusters = _doCluster(elements, dataset, maxDist)
    for cluster in clusters:
        cluster.txnames = txnames
    return clusters


def _doCluster(elements, dataset, maxDist):
    """
    Cluster algorithm to cluster elements.

    :parameter elements: list of all elements to be clustered
    :parameter dataset: Dataset object to be used when computing distances in upper limit space
    :parameter maxDist: maximum distance for clustering two elements

    :returns: a list of ElementCluster objects containing the elements
              belonging to the cluster
    """

    # First make sure all elements contain their upper limits
    for el in elements:
        if not hasattr(el,'._upperLimit'):
            el._upperLimit = dataset.getUpperLimitFor(el,txnames=el.txname)
        if el._upperLimit is None:
            raise SModelSError("Trying to cluster element outside the grid.")

    #Index elements:
    elementList = sorted(elements, key = lambda el: el._upperLimit)
    for iel,el in enumerate(elementList):
        el._index = iel

    #Start building maximal clusters
    clusterList = []
    for el in elementList:
        cluster = ElementCluster([],maxDist,dataset)
        for elB in elementList:
            if cluster.relativeDistance(el, elB) <= maxDist:
                cluster.add(elB)
        if not cluster.isConsistent():
            continue
        clusterList.append(cluster)

    #Split the maximal clusters until all elements inside each cluster are
    #less than maxDist apart from each other and the cluster average position
    #is less than maxDist apart from all elements
    finalClusters = []
    while clusterList:
        newClusters = []
        for cluster in clusterList:
            #Check if maximal internal distance is below maxDist
            if cluster.getMaxInternalDist() < maxDist:
                if not cluster in finalClusters:
                    finalClusters.append(cluster)

            #Cluster violates maxDist:
            else:
                #Loop over cluster elements and if element distance
                #falls outside the cluster, remove element
                for el in cluster:
                    if cluster.getDistanceTo(el) > maxDist:
                        newcluster = cluster.copy()
                        newcluster.remove(el)
                        if not newcluster.isConsistent():
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
        finalClusters.append(ElementCluster([el],maxDist,dataset))

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

    return finalClusters
