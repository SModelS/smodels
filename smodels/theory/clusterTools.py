"""
.. module:: theory.clusterTools
   :synopsis: Module holding the ElementCluster class and cluster methods used to combine similar elements according
   to the analysis.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

from smodels.theory import crossSection
from smodels.theory.auxiliaryFunctions import massAvg, massPosition, distance
from smodels.tools.physicsUnits import fb
import logging

logger = logging.getLogger(__name__)


class ElementCluster(object):
    """
    An instance of this class represents a cluster.
    
    This class is used to store the relevant information about a cluster of
    elements and to manipulate this information.
    
    """
    def __init__(self):
        self.elements = []


    def __iter__(self):
        return iter(self.elements)


    def __getitem__(self, iel):
        return self.elements[iel]


    def getTotalXSec(self):
        """
        Return the sum over the cross-sections of all elements belonging to
        the cluster.
        
        """
        totxsec = crossSection.XSectionList()
        for el in self.elements:
            totxsec.combineWith(el.weight)
        return totxsec

    def getAvgMass(self):
        """
        Return the average mass of all elements belonging to the cluster.
        
        :returns: average mass
         
        """
        massList = [el.getMasses() for el in self.elements]
        weights = [el.weight.getMaxXsec() / fb for el in self.elements]
        return massAvg(massList,weights=weights)


class IndexCluster(object):
    """
    An instance of this class represents a cluster storing element indices.
    
    This auxiliary class is used to store element indices and positions in
    upper limit space. It is used by the clustering algorithm.
    
    """
    def __init__(self, massMap=None, posMap=None, wMap=None, indices=set([]),
                 analysis=None):
        self.indices = indices
        self.avgPosition = None
        # Store the element mass list (for all elements)
        self.massMap = massMap
        # Store the element mass position in upper limit space list (for all
        # elements)
        self.positionMap = posMap
        self.weightMap = wMap
        self.analysis = analysis


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
        newcluster.analysis = self.analysis
      
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
        avgPos = massPosition(clusterMass, self.analysis)
        return avgPos


    def _getDistanceTo(self, obj):
        """
        Return the maximum distance between any elements belonging to the
        cluster and the object obj.
        
        obj can be a position in upper limit space or an element index.
        
        """
        dmax = 0.
        if type(obj) == type(int()) and obj >= 0:
            pos = self.positionMap[obj]
        elif type(obj) == type(float()):
            pos = obj
        else:
            logger.error("Unknown object type (must be an element index or "
                         "position)")
            import sys
            sys.exit()

        for jel in self:
            dmax = max(dmax, distance(pos, self.positionMap[jel]))
        return dmax


    def _getMaxInternalDist(self):
        """
        Return the maximum distance between any pair of elements belonging
        to the cluster as well as the cluster center and any element.
        
        """
        dmax = 0.
        if self.avgPosition == None:
            self.avgPosition = self._getAvgPosition()
        for iel in self:
            dmax = max(dmax, distance(self.positionMap[iel], self.avgPosition))
            dmax = max(dmax, self._getDistanceTo(iel))
        return dmax


def groupAll(elements):
    """
    Create a single cluster containing all the elements.
    
    """
    cluster = ElementCluster()
    cluster.elements = elements
    return cluster


def clusterElements(elements, analysis, maxDist):
    """
    Cluster the original elements according to their mass distance.
    
    :parameter elements: list of elements (Element objects)
    :parameter analysis: analysis to be considered (must be a ULanalysis object)
    :returns: list of clusters (ElementCluster objects)    
    """
    # Get the list of elements with good masses (with the masses replaced by
    # their 'good' value):
    goodElements = _getGoodElements(elements, analysis, maxDist)
    if len(goodElements) == 0:
        return []
    # ElementCluster elements by their mass:
    clusters = _doCluster(goodElements, analysis, maxDist)
    return clusters


def _doCluster(elements, analysis, maxDist):
    """
    Cluster algorithm to cluster elements.
    
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
            posMap[iel] = massPosition(massMap[iel], analysis)
            weightMap[iel] = el.weight.getMaxXsec() / fb
        else:
            j = massMap.keys()[massMap.values().index(el.getMasses())] 
            weightMap[j] += el.weight.getMaxXsec() / fb

    # Start with maximal clusters
    clusterList = []
    for iel in posMap:
        indices = [iel]
        for jel in posMap:            
            if distance(posMap[iel], posMap[jel]) <= maxDist:
                indices.append(jel)        
        indexCluster = IndexCluster(massMap, posMap, weightMap, set(indices), analysis)
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
                                           set([iel]),analysis))

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


def _getGoodElements(elements, analysis, maxDist):
    """
    Get the list of good masses appearing elements according to the analysis
    distance.
    
    :returns: list of elements with good masses with their masses replaced by
    the branch average. A mass is good if the mass distance between the
    branches is less than maxDist and if the element mass (or mass avg) falls
    inside the upper limit plane.
    
    """
    goodElements = []
    for element in elements:
        mass = element.getMasses()
        goodmass = None
        if mass[0] == mass[1]:
            p=massPosition(mass, analysis)
            if type(p)==type(fb):
                goodmass = mass
        else:
            mass1 = [mass[0], mass[0]]
            mass2 = [mass[1], mass[1]]
            mP1 = massPosition(mass1, analysis)
            mP2 = massPosition(mass2, analysis)
            if type(mP1)==type(None) or type(mP2)==type(None):
                continue
            if distance(mP1, mP2) < maxDist:
                goodmass = massAvg([mass1, mass2], method='harmonic')
        if goodmass and type(massPosition(goodmass, analysis))==type(fb):
            goodElements.append(element.copy())
            goodElements[-1].setMasses(goodmass)

    return goodElements
