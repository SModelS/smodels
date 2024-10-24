"""
.. module:: clusterTools
   :synopsis: Module holding the SMSCluster class and cluster methods used to combine similar SMS according
      to the analysis.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.base.particleNode import ParticleNode
from smodels.base.particle import MultiParticle
from smodels.base.physicsUnits import fb
from smodels.decomposition.theorySMS import TheorySMS
from smodels.experiment.datasetObj import DataSet, CombinedDataSet
from smodels.matching.exceptions import SModelSMatcherError as SModelSError
from smodels.matching.matcherAuxiliaryFuncs import average
from smodels.base.smodelsLogging import logger
import numpy as np



class AverageSMS(TheorySMS):
    """
    Represents an SMS or list of SMS containing only
    the basic attributes required for clustering or computing efficiencies/upper limits.
    Its properties are given by the average properties of the SMS
    it represents and its weight is given by the total weight of
    all SMS.
    """

    def __init__(self, smsList=[]):
        if any(not isinstance(sms, TheorySMS) for sms in smsList):
            raise SModelSError("An AverageSMS must be created from a list of TheorySMS objects.")

        TheorySMS.__init__(self)

        # Get the relevant properties needed by the txnames
        # (in addition to mass, totalwidth and isSM)
        attrList = ['mass', 'totalwidth', 'isSM']
        for sms in smsList:
            dataMap = sms.txname.dataMap
            attrList += [attr for _, attr,_ in dataMap.values()]
        attrList = list(set(attrList))

        # Define base SMS to copy (common) global properties from
        smsBase = smsList[0]

        #  Define relevant properties to be stored and averaged over:
        self.properties = attrList
        self.txname = smsBase.txname
        #If the average corresponds to a single SMS, we
        # can directly use its upperLimit value (if defined).
        if (len(smsList) == 1) and hasattr(smsBase,'_upperLimit'):
            self._upperLimit = smsBase._upperLimit

        #  Consistency checks:
        if any(sms.canonName != smsBase.canonName for sms in smsList):
            logger.error("Can not build average SMS from SMS with distinct topologies.")
            raise SModelSError()
        if any(sms.txname != self.txname for sms in smsList):
            logger.error("Can not build average SMS from SMS for distinct txnames.")
            raise SModelSError()


        # Create new node objects holding particles
        # with the average attributes
        avgNodesDict = {}
        if len(smsList) == 1:
            avgNodesDict = {n : node for n,node in zip(smsBase.nodeIndices,smsBase.nodes)}
        else:
            weights = [sms.weight.asNumber(fb) for sms in smsList]
            for nodeIndex in smsBase.nodeIndices:
                if nodeIndex == smsBase.rootIndex:
                    # For the root node make a dummy copy:
                    avgNode = smsBase.indexToNode(nodeIndex).copy()
                else:
                    # For all the other nodes compute the average
                    allParticles = [sms.indexToNode(nodeIndex).particle for sms in smsList]
                    attrDict = {}
                    for attr in self.properties:
                        values = [getattr(ptc, attr) for ptc in allParticles]
                        avgAttr = average(values, weights, nround=5)
                        attrDict[attr] = avgAttr
                    mp = allParticles[0]
                    for ptc in allParticles[1:]:
                        mp = mp + ptc
                    if isinstance(mp,MultiParticle):
                        avgParticle = MultiParticle(particles=mp.particles,
                                                attributesDict=attrDict)
                    else:
                        avgParticle = MultiParticle(particles=[mp],
                                                attributesDict=attrDict)
                    avgNode = ParticleNode(particle=avgParticle)
                avgNodesDict[nodeIndex] = avgNode

        # Copy the tree structure from smsBase,
        # but replacing the node objects by the average nodes
        self.copyTreeFrom(smsBase,nodesObjDict=avgNodesDict)
        # Compute the weight
        self.weight = smsBase.weight
        for sms in smsList[1:]:
            self.weight = self.weight + sms.weight

    def __eq__(self, other):
        """
        Compares the SMS with other. Only the properties
        defined in self.properties are used for comparison.
        :param other:  SMS to be compared (SMS object)
        :return: True/False.
        """

        if not isinstance(other, TheorySMS):
            return False
        
        if self.canonName != other.canonName:
            return False

        for nodeIndex in self.nodeIndices:
            nodeA = self.indexToNode(nodeIndex)
            nodeB = other.indexToNode(nodeIndex)
            if nodeA.isSM and nodeB.isSM:
                continue
            for attr in self.properties:
                attrA = getattr(nodeA,attr)
                attrB = getattr(nodeB,attr)
                if attrA != attrB:
                    return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

class SMSCluster(object):
    """
    An instance of this class represents a cluster of SMS,
    which holds the averageSMS.
    This class is used to store the relevant information about a cluster of
    SMS and to manipulate this information.
    """

    def __init__(self, smsList=[], dataset=None):

        self.smsList = smsList[:]
        self.dataset = dataset
        # Compute average SMS
        self.averageSMS = self.computeAverageSMS()
        
    def __str__(self):
        return str(self.averageSMS)

    def __repr__(self):
        return str(self.averageSMS)

    def getTotalXSec(self):
        """
        Return the sum over the cross sections of all SMS belonging to
        the cluster.

        :returns: sum of weights of all the SMS in the cluster (XSectionList object)
        """
        totxsec = 0.0*fb
        for sms in self.smsList:
            totxsec += sms.weight
        return totxsec

    def computeAverageSMS(self):
        """
        Computes the average SMS for the cluster.
        The average SMS has generic ParticleNodes
        with the attributes set to the average values of self.smsList.
        It can only be defined if all SMS share the same canonical name
        and the same txname. Otherwise, returns None

        :return: AverageSMS object or None (if it can not be defined)
        """

        # Check if an average SMS can be defined:
        cName = self.smsList[0].canonName
        tx = self.smsList[0].txname
        if any(sms.canonName != cName for sms in self.smsList):
            avgSMS = None
        elif any(sms.txname != tx for sms in self.smsList):
            avgSMS = None
        else:
            # Define the average SMS with the required properties averaged over:
            avgSMS = AverageSMS(self.smsList[:])

        # If the averageSMS can be defined, 
        # compute the averageSMS upper limit 
        if avgSMS is not None and self.dataset is not None:
            ul = self.dataset.getUpperLimitFor(avgSMS,txnames=tx)
            avgSMS._upperLimit = ul

        return avgSMS
    
    def distanceTo(self,sms):
        """
        Defines the relative distance between the cluster
        and SMS object or another cluster.
        The distance is computed using the upper limits for
        the averageSMS (for cluster objects) or the SMS upper limit.
        The distance is defined as d = 2*|ul1-ul2|/(ul1+ul2).

        :parameter sms: SMS object or SMSCluster object

        :returns: relative distance
        """

        sms1 = self.averageSMS

        if isinstance(sms,SMSCluster):
            sms2 = sms.averageSMS
        else:
            sms2 = sms

        if not hasattr(sms2, '_upperLimit'):
            sms2._upperLimit = self.dataset.getUpperLimitFor(sms2,
                                                    txnames=sms2.txname)

        ul1 = sms1._upperLimit
        ul2 = sms2._upperLimit

        if ul1 is None or ul2 is None:
            return None
        if (ul1+ul2).asNumber(fb) == 0.:
            return 0.
        ulDistance = 2.*abs(ul1 - ul2)/(ul1 + ul2)

        return ulDistance

    def isValid(self,maxDist):
        """
        Checks if the SMSCluster is a valid cluster,
        i.e. if its AverageSMS has a well defined UL
        and the distance between any SMS in the cluster
        and the cluster AverageSMS  is smaller than maxDist.

        :param dataset: Dataset object used to check the UL
        :parameter maxDist: maximum distance between the averageSMS and
                            all the SMS in the cluster

        :return: True/False 
        """

        if self.averageSMS is None:
            return False
        if self.averageSMS._upperLimit is None:
            return False
      
        # Compute distances between averageSMS and SMS in cluster
        dists = [self.distanceTo(sms) for sms in self.smsList]
        
        # Check if any distance could not be computed:
        if any(d is None for d in dists):
            return False
        
        # Compute maximal distance between averageSMS and SMS
        dist = max(dists,default=0.0)
        if dist > maxDist:
            return False
        
        return True



def clusterSMS(smsList, maxDist, dataset):
    """
    Cluster the original SMS according to their distance in upper limit space.

    :parameter smsList: list of sms (TheorySMS objects)
    :parameter dataset: Dataset object to be used when computing distances in upper limit space
    :parameter maxDist: maximum distance for clustering two SMS

    :returns: list of clusters (SMSCluster objects)
    """
    if len(smsList) == 0:
        return []

    if any(not isinstance(sms, TheorySMS) for sms in smsList):
        raise SModelSError("Asked to cluster non TheorySMS objects")
    if not isinstance(dataset, (DataSet, CombinedDataSet)):
        raise SModelSError("A dataset object must be defined for clustering")

    # Make sure only unique SMS are clustered together (avoids double counting weights)
    # Sort SMS, so the ones with highest contribution (weight) come first:
    smsList = sorted(smsList, key=lambda sms: sms.weight, reverse=True)

    # Remove duplicated SMS:
    smsUnique = []
    for sms in smsList:
        # Skip the SMS if it is related to any another SMS in the list
        if any(sms.isRelatedTo(smsB) for smsB in smsUnique):
            continue
        smsUnique.append(sms)

    # Get txname list only with the txnames from unique SMS used for clustering
    txnames = list(set([sms.txname for sms in smsUnique]))
    if dataset.getType() == 'upperLimit' and len(txnames) != 1:
        logger.error("Clustering SMS with different Txnames for an UL result.")
        raise SModelSError()

    if dataset.getType() == 'upperLimit':  # Group according to upper limit values
        clusters = doCluster(smsUnique, dataset, maxDist)
    elif dataset.getType() == 'efficiencyMap':  # Group all SMS together
        cluster = SMSCluster(smsUnique)
        clusters = [cluster]

    for cluster in clusters:
        cluster.txnames = txnames
    return clusters


def doCluster(smsList, dataset, maxDist):
    """
    Cluster algorithm to cluster SMS using a modified minimal spanning tree method.

    :parameter smsList: list of all SMS to be clustered
    :parameter dataset: Dataset object to be used when computing distances in upper limit space
    :parameter maxDist: maximum distance for clustering two SMS

    :returns: a list of SMSCluster objects containing the SMS
              belonging to the cluster
    """


    #  First make sure all SMS contain their upper limits
    for sms in smsList:
        if not hasattr(sms, '._upperLimit'):
            sms._upperLimit = dataset.getUpperLimitFor(sms, txnames=sms.txname)
        if sms._upperLimit is None:
            raise SModelSError("Trying to cluster SMS (id = %i, txname = %s) outside the grid for dataset %s." 
                                %(sms.smsID,sms.txname,dataset.longStr()))

    # Sort SMS by upperLimit 
    sortedSMSList = sorted(smsList, key = lambda sms: sms._upperLimit)

    # Cluster all identical SMS into clusters
    clusterList = groupSMS(sortedSMSList,dataset)

    # Sort clusters by proximity in upperLimit 
    clusterList = sorted(clusterList, 
                         key=lambda cluster: cluster.averageSMS._upperLimit)

    # Compute the distance matrix for the clusters:
    # (diagonal entries are set to None)
    dMatrix = np.full((len(clusterList),len(clusterList)),fill_value=None)
    for iA,clusterA in enumerate(clusterList):
        for iB,clusterB in enumerate(clusterList):
            if iA >= iB: continue
            dMatrix[iA,iB] = clusterA.distanceTo(clusterB)
            dMatrix[iB,iA] = dMatrix[iA,iB] # the matrix is symmetric
    
    # Check minimal distance between two distinct clusters
    # (set default to 2*maxDist in case len(clusterList) == 1, i.e. dMatrix is empty)
    minDist = min([d for d in dMatrix.flatten() if d is not None], 
                  default=2*maxDist)
    # Merge closest cluster up to the point when
    # no more merges are possible (all clusters have distances larger than maxDist)
    while (len(clusterList) > 1) and (minDist < maxDist):
        # Get indices of first pair of clusters with minimum distance
        mergeIndices = np.argwhere(dMatrix == minDist)[0]
        # Merge clusters with indices in mergeIndices
        newCluster = mergeClusters(np.take(clusterList,mergeIndices))
        # If new cluster is valid (averageSMS is close to the cluster's SMS),
        # replace the two orignal clusters by it. 
        if newCluster.isValid(maxDist):
            # Update cluster list:
            # Remove cluster with largest index
            # and replace cluster with smallest index by new cluster
            delete_index = max(mergeIndices)
            replace_index = min(mergeIndices)
            clusterList.pop(delete_index)
            clusterList[replace_index] = newCluster 
            # Remove row and column for cluster with largest index
            dMatrix = np.delete(dMatrix,delete_index,0)
            dMatrix = np.delete(dMatrix,delete_index,1)
            # Compute the distances between the new cluster and all the
            # other ones.
            # Recompute distances using the new cluster:
            for iA,clusterA in enumerate(clusterList):
                if iA >= replace_index: continue
                dMatrix[iA,replace_index] = clusterA.distanceTo(newCluster)
                dMatrix[replace_index,iA] = dMatrix[iA,replace_index] # the matrix is symmetric
        
        # Otherwise keep the original clusters, but set 
        # their distance above maxDist, so they will no longer be allowed to merge
        else:
            dMatrix[mergeIndices] = 2*maxDist
            dMatrix[np.flip(mergeIndices)] = 2*maxDist


        minDist = min([d for d in dMatrix.flatten() if d is not None], 
                  default=2*maxDist)
        

    # Finally sort clusters by total xsection and length
    clusters = sorted(clusterList, 
                      key = lambda c: (c.getTotalXSec(),len(c.smsList)),
                      reverse=True)

    return clusters


def groupSMS(smsList,dataset):
    """
    Group SMS into clusters where the average SMS is
    identical to all the SMS in cluster.
    Each cluster contains all SMS which share the same mass,width and upper limit
    and can be replaced by their average SMS.

    :parameter smsList: list of all SMS to be grouped

    :returns: a list of SMSCluster objects
              which represents a group of SMS with same mass, width and upper limit.
    """

    # Group SMS if they have the same UL
    # and give the same average SMS (same BSM attributes)
    smsClusterList = []
    nonClusteredSMS = smsList[:]
    while nonClusteredSMS:
        smsA = nonClusteredSMS.pop(0)
        smsCluster = SMSCluster([smsA],dataset=dataset)
        # Keep track of the SMS which could not
        # be clustered with smsA to use in the next iteration
        notClustered = []
        while nonClusteredSMS:
            smsB = nonClusteredSMS.pop(0)
            # If it can not be clustered with smsA, add to notClustered list
            if (smsA._upperLimit != smsB._upperLimit) or (smsCluster.averageSMS != smsB):
                notClustered.append(smsB)
                continue
            smsCluster.smsList.append(smsB)
            smsCluster.averageSMS.weight += smsB.weight
        
        smsClusterList.append(smsCluster)
        nonClusteredSMS = notClustered[:]

    # Make sure each SMS belongs to a cluster:
    nclustered = sum([len(c.smsList) for c in smsClusterList])
    nsms = len(smsList)
    if nclustered != nsms:
        raise SModelSError(f"Error grouping SMS. {nclustered}/{nsms} SMS were clustered")
    
    return smsClusterList

def mergeClusters(clusterList,useAverage=False):
    """
    Merge a list of SMSCluster objects using
    the averageSMS of each cluster if useAverage is True.
    Otherwise cluster the list of all SMS from all clusters.

    :param clusterList: List of SMSCluster objects
    :param useAverage: True/False. If True, 
                       will cluster the averageSMS in each cluster.
    """
        
    if len(clusterList) == 0:
        return None
    
    smsList = []
    for cluster in clusterList:
        smsList += cluster.smsList[:]
        
    if not useAverage:
        newCluster = SMSCluster(smsList,
                                clusterList[0].dataset)
    else:
        avgSMSList = [cluster.averageSMS 
                   for cluster in clusterList]
        newCluster = SMSCluster(avgSMSList,
                                clusterList[0].dataset)
        # Make sure the cluster SMS list contains the original SMS
        newCluster.smsList = smsList[:]

    return newCluster
