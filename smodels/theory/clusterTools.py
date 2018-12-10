"""
.. module:: clusterTools
   :synopsis: Module holding the ElementCluster class and cluster methods used to combine similar elements according
      to the analysis.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

from smodels.theory import crossSection
from smodels.theory.auxiliaryFunctions import massAvg, distance
from smodels.tools.physicsUnits import fb
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
import numpy as np

from smodels.tools.smodelsLogging import logger


class ElementCluster(object):
    """
    An instance of this class represents a cluster.    
    This class is used to store the relevant information about a cluster of
    elements and to manipulate this information.
    
    :ivar elements: list of elements in the cluster (Element objects)    
    """
    def __init__(self):
        self.elements = []

    def __iter__(self):
        return iter(self.elements)

    def __getitem__(self, iel):
        return self.elements[iel]

    def getTotalXSec(self):
        """
        Return the sum over the cross sections of all elements belonging to
        the cluster.
        
        :returns: sum of weights of all the elements in the cluster (XSectionList object)
        """
        totxsec = crossSection.XSectionList()
        for el in self.elements:
            totxsec.combineWith(el.weight)
        if len(totxsec) != 1:
            logger.error("Cluster total cross section should have a single value")
            raise SModelSError()
        return totxsec[0]

    def getAvgMass(self):
        """
        Return the average mass of all elements belonging to the cluster, weighted
        by their individual weights.
        If the masses have distinct shapes (e.g. for different txnames in efficiencyMap
        results), return None.
        
        :returns: average mass array appearing in the cluster     
        """                          
        
        massList = [el.getMasses() for el in self.elements]
        if any(np.array(m).shape != np.array(massList[0]).shape for m in massList):
            return None
        
        weights = [el.weight.getMaxXsec().asNumber(fb) for el in self.elements]

        return massAvg(massList,weights=weights)


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

    def getDataType(self):
        """
        Checks to which type of data (efficiency map or upper limit)
        the cluster refers to. It uses the cluster.txnames attribute.
        If not defined, returns None
        :return: upperLimits or efficiencyMap (string)
        """

        if not hasattr(self, 'txnames') or not self.txnames:
            return None
        else:
            #Check the data types
            dataType = list(set([txname.txnameData.dataType for txname in self.txnames]))            
            if len(dataType) != 1:
                logger.error("A single cluster contain mixed data types!")
                raise SModelSError()
            elif 'upperLimit' in dataType[0]:
                return 'upperLimit'
            elif 'efficiencyMap' in dataType[0]:
                return 'efficiencyMap'
            else:
                logger.error("Unknown data type %s" % (dataType[0]))
                raise SModelSError()


class IndexCluster(object):
    """
    An instance of this class represents a cluster storing element indices.
    This auxiliary class is used to store element indices and positions in
    upper limit space. It is only used by the clustering algorithm.

    :ivar indices: list of integers mapping the cluster elements to their position in the list
                   (1st element -> index 0, 2nd element -> index 1,...)
    :ivar elementMap: dictionary with indices as keys and the corresponding element as values
    :ivar ulMap: dictionary with indices as keys and the corresponding element upper limit as values
    :ivar weightMap: dictionary with indices as keys and the corresponding element weight
                     as values
    :ivar txdata: TxNameData object to be used for computing distances in UL space
    """
    def __init__(self, elementMap=None, ulMap=None, wMap=None, indices=set([]), txdata = None):
        self.indices = indices
        self.elementMap = elementMap
        self.positionMap = ulMap
        self.weightMap = wMap
        self.txdata = txdata


    def __eq__(self, other):
        if type(self) != type(other):
            return False
        elif set(self.indices) != set(other.indices):
            return False
        else:
            return True


    def __iter__(self):
        return iter(list(self.indices))


    def __getitem__(self, iel):
        return list(self.indices)[iel]

    def copy(self):
        """
        Returns a copy of the index cluster (faster than deepcopy).
        """

        newcluster = IndexCluster()
        newcluster.indices = set(list(self.indices)[:])
        if type(self.positionMap) == type(dict()):
            newcluster.positionMap = dict(self.positionMap.items())
        else: newcluster.positionMap = None
        if type(self.elementMap) == type(dict()):
            newcluster.maselementMapsMap = dict(self.elementMap.items())
        else: newcluster.elementMap = None
        if type(self.weightMap) == type(dict()):
            newcluster.weightMap = dict(self.weightMap.items())
        else: newcluster.weightMap = None
        newcluster.txdata = self.txdata

        return newcluster


    def add(self, iels):
        """
        Add an index or a list of indices to the list of indices and update
        the avgPosition value.

        """
        if type(iels) == type(int()):
            ilist = [iels]
        else:
            ilist = iels

        indices = list(self.indices).extend(ilist)
        self.indices = set(indices)

    def remove(self, iels):
        """
        Remove an index or a list of indices to the list of indices and
        update the avgPosition value.

        """
        if type(iels) == type(int()):
            ilist = [iels]
        else:
            ilist = iels

        indices = list(self.indices)
        for iel in ilist:
            indices.remove(iel)
        self.indices = set(indices)


    def _getDistanceTo(self, elIndex):
        """
        Return the maximum distance between any elements belonging to the
        cluster and the element object corresponding to the element index elIndex.

        :parameter elIndex: Element index in elementMap.
        :return: maximum distance (float)
        """

        if not isinstance(elIndex,int):
            logger.error("Unknown object type (must be an element index")
            raise SModelSError()
        elif not elIndex in self.positionMap:
            logger.error("Element index %i not found in postion map" %elIndex)
            raise SModelSError()

        pos = self.positionMap[elIndex]
        el = self.elementMap[elIndex]
        dmax = max([distance(pos, self.positionMap[iel],el,self.elementMap[iel]) for iel in self])

        return dmax


    def _getMaxInternalDist(self):
        """
        Return the maximum distance between any pair of elements belonging
        to the cluster.
        """

        dmax = max([self._getDistanceTo(iel) for iel in self])

        return dmax


def groupAll(elements):
    """
    Create a single cluster containing all the elements.
    Skip mother elements which contain the daughter in the list (avoids double counting).

    :param elements: List of Element objects
    :return: ElementCluster object containing all unique elements
    """


    cluster = ElementCluster()
    cluster.elements = []
    allmothers = []
    #Collect the list of all mothers:
    for el in elements:
        allmothers += [elMom[1].elID for elMom in el.motherElements if not elMom[0]=='original']

    for el in elements:
        #Skip the element if it is a mother of another element in the list
        if any((elMom is el.elID) for elMom in allmothers):
            continue
        cluster.elements.append(el)

    #Collect the txnames appearing in the cluster
    cluster.txnames = list(set([el.txname for el in cluster.elements]))
    cluster.txnames.sort()
    return cluster


def clusterElements(elements, maxDist):
    """
    Cluster the original elements according to their distance in mass and lifetime space.

    :parameter elements: list of elements (Element objects)
    :parameter txname: TxName object to be used for computing distances in UL space
    :parameter maxDist: maximum distance for clustering two elements

    :returns: list of clusters (ElementCluster objects)
    """
    if len(elements) == 0:
        return []
    txnames = list(set([el.txname for el in elements]))
    if len(txnames) != 1:
        logger.error("Clustering elements with different Txnames!")
        raise SModelSError()
    txdata = txnames[0].txnameData
    # ElementCluster elements:
    clusters = _doCluster(elements, txdata, maxDist)
    for cluster in clusters:
        cluster.txnames = txnames
    return clusters


def _doCluster(elements, txdata, maxDist):
    """
    Cluster algorithm to cluster elements.

    :parameter elements: list of all elements to be clustered
    :parameter txdata: TxNameData object to be used for computing distances in UL space
    :parameter maxDist: maximum distance for clustering two elements

    :returns: a list of ElementCluster objects containing the elements
              belonging to the cluster
    """
    # First build the index:element, index:upperLimit space
    # and element:maxWeight (in fb) dictionaries
    #(Combine elements with identical masses)
    elementMap = {}
    ulMap = {}
    weightMap = {}
    for iel, el in enumerate(elements):
        elementMap[iel] = el
        ulMap[iel] = txdata.getValueFor(el).asNumber(fb)
        weightMap[iel] = el.weight.getMaxXsec().asNumber(fb)

    # Start with maximal clusters
    clusterList = []
    for iel in ulMap:
        indices = [iel]
        for jel in ulMap:
            if distance(ulMap[iel], ulMap[jel], elementMap[iel], elementMap[jel]) <= maxDist:
                indices.append(jel)
        indexCluster = IndexCluster(elementMap, ulMap, weightMap, set(indices),txdata)
        clusterList.append(indexCluster)

    #Split the maximal clusters until all elements inside each cluster are
    #less than maxDist apart from each other
    finalClusters = []
    newClusters = True
    while newClusters:
        newClusters = []
        for indexCluster in clusterList:
            # cluster is good
            if indexCluster._getMaxInternalDist() < maxDist:
                if not indexCluster in finalClusters:
                    finalClusters.append(indexCluster)
                continue
            #Loop over cluster elements and if element distance is greater than maxDist
            #remove element
            for iel in indexCluster:
                dist = indexCluster._getDistanceTo(iel)
                if dist > maxDist:
                    newcluster = indexCluster.copy()
                    newcluster.remove(iel)
                    if not newcluster in newClusters:
                        newClusters.append(newcluster)

        clusterList = newClusters
        # Check for oversized list of indexCluster (too time consuming)
        if len(clusterList) > 100:
            logger.warning("ElementCluster failed, using unclustered masses")
            finalClusters = []
            clusterList = []

    # finalClusters = finalClusters + clusterList
    # Add clusters of individual masses (just to be safe)
    for iel in elementMap:
        finalClusters.append(IndexCluster(elementMap, ulMap, weightMap,
                                           set([iel])))

    # Clean up clusters (remove redundant clusters)
    for ic, clusterA in enumerate(finalClusters):
        if clusterA is None:
            continue
        for jc, clusterB in enumerate(finalClusters):
            if clusterB is None:
                continue
            if ic != jc and clusterB.indices.issubset(clusterA.indices):
                finalClusters[jc] = None
    while finalClusters.count(None) > 0:
        finalClusters.remove(None)

    # Transform index clusters to element clusters:
    clusterList = []
    for indexCluster in finalClusters:
        cluster = ElementCluster()
        cluster.elements = [elementMap[iel] for iel in indexCluster]
        clusterList.append(cluster)

    return clusterList
