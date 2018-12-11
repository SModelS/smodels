"""
.. module:: clusterTools
   :synopsis: Module holding the ElementCluster class and cluster methods used to combine similar elements according
      to the analysis.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

from smodels.theory import crossSection
from smodels.theory.auxiliaryFunctions import distance, removeUnits, addUnit
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

    def __init__(self, elements = [], maxDist = 0., txdata = None):

        self.elements = elements
        self.maxDist = maxDist
        self.txdata = txdata

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
        
        massList = [el.getMasses() for el in self]
        if any(np.array(m).shape != np.array(massList[0]).shape for m in massList):
            return None
        
        weights = [el.weight.getMaxXsec().asNumber(fb) for el in self]
        massList = np.array(removeUnits(massList, [GeV]))
        mAvg = np.average(massList,weights=weights,axis=0)
        mAvg = addUnit(mAvg,unit=GeV)

        return mAvg

    def getAvgWidth(self):
        """
        Return the average width of all elements belonging to the cluster, weighted
        by their individual weights.
        If the widths have distinct shapes (e.g. for different txnames in efficiencyMap
        results), return None.

        :returns: average width array appearing in the cluster
        """                          

        widthList = [el.getWidths() for el in self]
        if any(np.array(w).shape != np.array(widthList[0]).shape for w in widthList):
            return None

        weights = [el.weight.getMaxXsec().asNumber(fb) for el in self]
        widthList = np.array(removeUnits(widthList, [GeV]))
        wAvg = np.average(widthList,weights=weights,axis=0)
        wAvg = addUnit(wAvg,unit=GeV)

        return wAvg

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

        newcluster = ElementCluster(self.elements,self.indices,self.maxDist)
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
            element._upperLimit = self.txdata.getValueFor(element)
        if element._upperLimit is None:
            return None

        dmax = 0.
        for el in self:
            if not hasattr(el, '_upperLimit'):
                el._upperLimit = self.txdata.getValueFor(el)
            dmax = max(dmax,distance(element._upperLimit,el._upperLimit))

        return dmax

    def _getMaxInternalDist(self):
        """
        Return the maximum distance between any pair of elements belonging
        to the cluster.
        """

        dmax = max([self._getDistanceTo(el) for el in self])

        return dmax

    def averageElement(self):
        """
        Computes the average element for the cluster.
        The average element is identical to the first element of the cluster,
        except that the masses and widths of the oddParticles are replaced by
        their (weighted) average values over the cluster.

        :return: Element object
        """

        avgEl = self.elements[0].copy()
        avgMass = np.array(self.getAvgMass())
        avgWidth = np.array(self.getAvgWidth())
        for i,br in enumerate(avgEl.branches):
            for j,ptc in enumerate(br.oddParticles):
                avgParticle = ptc.copy()
                avgParticle.__setattr__('mass',avgMass[i,j])
                avgParticle.__setattr__('totalwidth',avgWidth[i,j])
                avgEl.branches[i].oddParticle[j]  = avgParticle

        avgEl._upperLimit = self.txdata.getValueFor(avgEl)

        return avgEl

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

    # First make sure all elements contain their upper limits
    for el in elements:
        if not hasattr(el,'._upperLimit'):
            el._upperLimit = txdata.getValueFor(el)
        if el._upperLimit is None:
            raise SModelSError("Trying to cluster element outside the grid.")

    #Index elements:
    elementList = sorted(elements, key = lambda el: el._upperLimit)
    for iel,el in enumerate(elementList):
        el._index = iel

    #Start building maximal clusters
    clusterList = []
    for el in elementList:
        maxels = [elB for elB in elementList if distance(elB,el) <= maxDist]
        cluster = ElementCluster(elements=maxels,maxDist=maxDist,txdata=txdata)
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
            if cluster._getMaxInternalDist() < maxDist:
                if not cluster in finalClusters:
                    finalClusters.append(cluster)

            #Cluster violates maxDist:
            else:
                #Loop over cluster elements and if element distance
                #falls outside the cluster, remove element
                for el in cluster:
                    if cluster._getDistanceTo(iel) > maxDist:
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
        finalClusters.append(ElementCluster(elements=[el],
                                            maxDist=maxDist,
                                            txdata=txdata))

    # Clean up clusters (remove redundant clusters)
    for ic, clusterA in enumerate(finalClusters):
        if clusterA is None:
            continue
        for jc, clusterB in enumerate(finalClusters):
            if clusterB is None:
                continue
            if ic != jc and clusterB.indices().issubset(clusterA.indices()):
                finalClusters[jc] = None
    while finalClusters.count(None) > 0:
        finalClusters.remove(None)

    return finalClusters
