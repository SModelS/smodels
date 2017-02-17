"""
.. module:: clusterTools
   :synopsis: Module holding the ElementCluster class and cluster methods used to combine similar elements according
      to the analysis.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

import sys
from smodels.theory import crossSection
from smodels.theory.auxiliaryFunctions import massAvg, massPosition, distance
from smodels.tools.physicsUnits import fb, MeV
from smodels.theory.exceptions import SModelSTheoryError as SModelSError

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
        Return the average mass of all elements belonging to the cluster.
        If the cluster does not refer to a TxName (i.e. in efficiency map results)
        AND the cluster contains more than one element (assuming they differ in
        the masses), returns None.
        
        :returns: average mass array         
        """

        def similar ( a1, a2 ):
            if len(a1) != len(a2):
                return False
            for l1,l2 in zip ( a1,a2 ):
                if len(l1) != len(l2):
                    return False
                for e1, e2 in zip (l1,l2 ):
                    d=abs ( (e1-e2).asNumber(MeV) )
                    if d>.1:
                        return False
            return True
        
        if self.getDataType() == 'efficiencyMap':
            if len(self.elements) > 1: 
                ret = self.elements[0].getMasses()
                for i in self.elements[1:]:
                    if not similar ( i.getMasses(), ret ):
                        return None
                return ret
            else: return self.elements[0].getMasses()
        elif self.getDataType() == 'upperLimit':
            massList = [el.getMasses() for el in self.elements]
            weights = [el.weight.getMaxXsec() / fb for el in self.elements]        
            return massAvg(massList,weights=weights)


    def getPIDs(self):
        """
        Return the list of all PIDs appearing in all elements in the cluster,
        i.e. [ [[pdg1, pdg2,...],[pdg3,pdg4,...]], [[pdg1', pdg2',...],[pdg3',pdg4',...]] 
        
        :returns: list of PIDs
        """
        
        PIDs = []
        for el in self:
            for pidList in el.getPIDs():
                if not pidList in PIDs: PIDs.append(pidList)
            
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
            #for txname in self.txnames:                
            #    if not txname.txnameData._data:
            #        txname.txnameData.loadData()  #Make sure the _data is loaded
                    
            dataTag = list(set([txname.txnameData.dataTag for txname in self.txnames]))            
            if len(dataTag) != 1:
                logger.error("A single cluster contain mixed data types!")
                raise SModelSError()
            elif 'upperLimit' in dataTag[0]:
                return 'upperLimit'
            elif 'efficiencyMap' in dataTag[0]:
                return 'efficiencyMap'
            else:
                logger.error("Unknown data type %s" % (dataTag[0]))
                raise SModelSError()


class IndexCluster(object):
    """
    An instance of this class represents a cluster storing element indices.    
    This auxiliary class is used to store element indices and positions in
    upper limit space. It is only used by the clustering algorithm.
    
    :ivar indices: list of integers mapping the cluster elements to their position in the list
                   (1st element -> index 0, 2nd element -> index 1,...)
    :ivar avgPosition: position in upper limit space for the cluster average mass
    :ivar massMap: dictionary with indices as keys and the corresponding element mass as values
    :ivar positionMap: dictionary with indices as keys and the corresponding element position
                        in upper limit space as values
    :ivar weightMap: dictionary with indices as keys and the corresponding element weight
                     as values
    :ivar txdata: TxNameData object to be used for computing distances in UL space
    """
    def __init__(self, massMap=None, posMap=None, wMap=None, indices=set([]), txdata = None):
        self.indices = indices
        self.avgPosition = None
        self.massMap = massMap
        self.positionMap = posMap
        self.weightMap = wMap
        self.txdata = txdata


        if massMap and posMap and wMap and len(self.indices) > 0:
            self.avgPosition = self._getAvgPosition()


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
        newcluster.avgPosition = self.avgPosition
        if type(self.positionMap) == type(dict()):
            newcluster.positionMap = dict(self.positionMap.items())
        else: newcluster.positionMap = None
        if type(self.massMap) == type(dict()):
            newcluster.massMap = dict(self.massMap.items())
        else: newcluster.massMap = None
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
        self.avgPosition = self._getAvgPosition()


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
        self.avgPosition = self._getAvgPosition()


    def _getAvgPosition(self):
        """
        Return the average position in upper limit space for all indices
        belonging to the cluster.
        
        """
        if len(list(self.indices)) == 1:
            return self.positionMap[self[0]]
        masses = [self.massMap[iel] for iel in self]
        weights = [self.weightMap[iel] for iel in self]
        clusterMass = massAvg(masses,weights=weights)
        avgPos = self.txdata.getValueFor(clusterMass)
        if avgPos is None:
            return False
        else:
            return (avgPos/fb).asNumber()

    def _getDistanceTo(self, obj):
        """
        Return the maximum distance between any elements belonging to the
        cluster and the object obj.
        
        obj can be a position in upper limit space or an element index.
        
        """

        dmax = 0.
        if type(obj) == type(int()) and obj >= 0:
            pos = self.positionMap[obj]
        elif type(obj) == type(1.):
            pos = obj
        else:
            logger.error("Unknown object type (must be an element index or "
                         "position)")
            raise SModelSError()

        for jel in self:
            dmax = max(dmax, distance(pos, self.positionMap[jel]))
        return dmax


    def _getMaxInternalDist(self):
        """
        Return the maximum distance between any pair of elements belonging
        to the cluster as well as the cluster center and any element.
        
        """
        dmax = 0.
        if self.avgPosition is None:
            self.avgPosition = self._getAvgPosition()
        for iel in self:
            dmax = max(dmax, distance(self.positionMap[iel], self.avgPosition))
            dmax = max(dmax, self._getDistanceTo(iel))
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
        allmothers += [elMom for tp,elMom in el.motherElements]
        
    for el in elements:
        #Skip the element if it is a mother of another element in the list
        if any((elMom is el) for elMom in allmothers):
            continue
        cluster.elements.append(el) 
    
    #Collect the txnames appearing in the cluster
    cluster.txnames = list(set([el.txname for el in cluster.elements]))
    cluster.txnames.sort()
    return cluster


def clusterElements(elements, maxDist):
    """
    Cluster the original elements according to their mass distance.
    
    :parameter elements: list of elements (Element objects)
    :parameter txname: TxName object to be used for computing distances in UL space
    :parameter maxDist: maximum mass distance for clustering two elements
    
    :returns: list of clusters (ElementCluster objects)    
    """
    if len(elements) == 0:  return []
    txnames = list(set([el.txname for el in elements]))
    if len(txnames) != 1:
        logger.error("Clustering elements with different Txnames!")
        raise SModelSError()
    txdata = txnames[0].txnameData
    # ElementCluster elements by their mass:
    clusters = _doCluster(elements, txdata, maxDist)
    for cluster in clusters:
        cluster.txnames = txnames
    return clusters


def _doCluster(elements, txdata, maxDist):
    """
    Cluster algorithm to cluster elements.
    
    :parameter elements: list of all elements to be clustered
    :parameter txdata: TxNameData object to be used for computing distances in UL space
    :parameter maxDist: maximum mass distance for clustering two elements
    
    :returns: a list of ElementCluster objects containing the elements
    belonging to the cluster    
    """
    # First build the element:mass, element:position in UL space
    # and element:maxWeight (in fb) dictionaries
    #(Combine elements with identical masses)
    massMap = {}
    posMap = {}
    weightMap = {}
    for iel, el in enumerate(elements):
        if not el.getMasses() in massMap.values():
            massMap[iel] = el.getMasses()
            posMap[iel] = massPosition(massMap[iel], txdata)
            weightMap[iel] = el.weight.getMaxXsec() / fb
        else:
            j = list(massMap.keys())[list(massMap.values()).index(el.getMasses())] 
            weightMap[j] += el.weight.getMaxXsec() / fb

    # Start with maximal clusters
    clusterList = []
    for iel in posMap:
        indices = [iel]
        for jel in posMap:            
            if distance(posMap[iel], posMap[jel]) <= maxDist:
                indices.append(jel)        
        indexCluster = IndexCluster(massMap, posMap, weightMap, set(indices),txdata)
        #Ignore cluster which average mass falls oustide the grid:
        if indexCluster.avgPosition:
            clusterList.append(indexCluster)

    #Split the maximal clusters until all elements inside each cluster are
    #less than maxDist apart from each other and the cluster average position
    #is less than maxDist apart from all elements
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
            # Distance to cluster center (average)
            distAvg = indexCluster._getDistanceTo(indexCluster.avgPosition)
            #Loop over cluster elements and if element distance or cluster
            #average distance falls outside the cluster, remove element
            for iel in indexCluster:
                dist = indexCluster._getDistanceTo(iel)
                if max(dist, distAvg) > maxDist:
                    newcluster = indexCluster.copy()
                    newcluster.remove(iel)
                    if not newcluster in newClusters:
                        #Ignore cluster which average mass falls oustide the grid:
                        if newcluster.avgPosition:
                            newClusters.append(newcluster)

        clusterList = newClusters
        # Check for oversized list of indexCluster (too time consuming)
        if len(clusterList) > 100:
            logger.warning("ElementCluster failed, using unclustered masses")
            finalClusters = []
            clusterList = []

    # finalClusters = finalClusters + clusterList
    # Add clusters of individual masses (just to be safe)
    for iel in massMap:
        finalClusters.append(IndexCluster(massMap, posMap, weightMap,
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
        masses = [massMap[iel] for iel in indexCluster]
        for el in elements:
            if el.getMasses() in masses:
                cluster.elements.append(el)
        clusterList.append(cluster)

    return clusterList
