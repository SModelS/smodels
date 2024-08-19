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
            attrList += [attr for node, attr, unit in dataMap.values()]
        attrList = list(set(attrList))

        # Define base SMS to copy (common) global properties from
        smsBase = smsList[0]

        #  Define relevant properties to be stored and averaged over:
        self.properties = attrList
        self.smsList = smsList[:]
        self.txname = smsBase.txname

        #  Consistency checks:
        if any(sms.canonName != smsBase.canonName for sms in self.smsList):
            logger.error("Can not build average SMS from SMS with distinct topologies.")
            raise SModelSError()
        if any(sms.txname != self.txname for sms in self.smsList):
            logger.error("Can not build average SMS from SMS for distinct txnames.")
            raise SModelSError()


        # Create new node objects holding particles
        # with the average attributes
        avgNodesDict = {}
        if len(smsList) == 1:
            avgNodesDict = {n : node for n,node in zip(smsBase.nodeIndices,smsBase.nodes)}
        else:
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
                        avgAttr = self.getAverage(values)
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
        for sms in self.smsList[1:]:
            self.weight = self.weight + sms.weight

    def __cmp__(self, other):
        """
        Compares the SMS with other. Only the properties
        defined in self.properties are used for comparison.
        :param other:  SMS to be compared (SMS object)
        :return: -1 if self < other, 0 if self == other, +1, if self > other.
        """

        if not isinstance(other, TheorySMS):
            return -1

        for nodeIndex in self.nodeIndices:
            nodeA = self.indexToNode(nodeIndex)
            nodeB = other.indexToNode(nodeIndex)
            if nodeA.isSM and nodeB.isSM:
                continue
            for attr in self.properties:
                attrA = getattr(nodeA,attr)
                attrB = getattr(nodeB,attr)
                if attrA != attrB:
                    if attrA > attrB:
                        return 1
                    else:
                        return -1
        return 0

    def __lt__(self, other):
        return self.__cmp__(other) == -1

    def __gt__(self, other):
        return self.__cmp__(other) == 1

    def __eq__(self, other):
        return self.__cmp__(other) == 0

    def __ne__(self, other):
        return self.__cmp__(other) != 0

    def getAverage(self, values, weighted=True, nround=5):
        """
        Compute the average value for the attribute using
        the SMS list in self.smsList.
        If weighted = True, compute the weighted average
        using the SMS weights.

        :param values: List of values to be averaged over
        :param weighted: If True, uses the SMS weights to compute a weighted average
        :param nround: If greater than zero and the returning attibute is numeric, will round it
                      to this number of significant digits.

        :return: Average value of attribute.
        """

        if weighted:
            weights = [sms.weight.asNumber(fb) for sms in self.smsList]
        else:
            weights = [1.]*len(self.smsList)

        return average(values, weights, nround)

    def contains(self, sms):
        """
        Check if the average SMS contains the SMS

        :param sms: TheorySMS object

        :return: True/False
        """

        if any(sms is selfSMS for selfSMS in self.smsList):
            return True
        return False


class SMSCluster(object):
    """
    An instance of this class represents a cluster of SMS.
    This class is used to store the relevant information about a cluster of
    SMS and to manipulate this information.
    """

    def __init__(self, smsList=[]):

        self.smsList = smsList

    def __eq__(self, other):

        if type(self) != type(other):
            return False
        elif set(self.indices()) != set(other.indices()):
            return False
        else:
            return True

    def __iter__(self):
        return iter(self.smsList)

    def __getitem__(self, isms):
        return self.smsList[isms]

    def __len__(self):
        return len(self.smsList)

    def __str__(self):
        return str(self.smsList)

    def __repr__(self):
        return str(self.smsList)

    def indices(self):
        """
        Return a list of SMS indices appearing in cluster
        """

        indices = [sms._index for sms in self]

        return indices

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

    @property
    def averageSMS(self):
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
            return None
        if any(sms.txname != tx for sms in self.smsList):
            return None
        # Define the average SMS with the required properties averaged over:
        avgSMS = AverageSMS(self.smsList[:])

        avgSMS._index = None

        return avgSMS

    def isValid(self,dataset):
        """
        Checks if the SMSCluster is a valid cluster,
        i.e. if its AverageSMS has a well defined UL

        :param dataset: Dataset object used to check the UL

        :return: True/False 
        """

        ul = dataset.getUpperLimitFor(self.averageSMS, 
                                      txnames=self.averageSMS.txname)
        self.averageSMS._upperLimit = ul

        isvalid = (ul is not None)
        
        return isvalid

def relativeDistance(sms1, sms2, dataset):
    """
    Defines the relative distance between two SMS according to their
    upper limit values.
    The distance is defined as d = 2*|ul1-ul2|/(ul1+ul2).

    :parameter sms1: SMS object
    :parameter sms2: SMS object

    :returns: relative distance
    """

    if not hasattr(sms1, '_upperLimit'):
        sms1._upperLimit = dataset.getUpperLimitFor(sms1,
                                                   txnames=sms1.txname)
    if not hasattr(sms2, '_upperLimit'):
        sms2._upperLimit = dataset.getUpperLimitFor(sms2,
                                                   txnames=sms2.txname)

    ul1 = sms1._upperLimit
    ul2 = sms2._upperLimit

    if ul1 is None or ul2 is None:
        return None
    if (ul1+ul2).asNumber(fb) == 0.:
        return 0.
    ulDistance = 2.*abs(ul1 - ul2)/(ul1 + ul2)

    return ulDistance

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
        cluster = SMSCluster()
        for isms, sms in enumerate(smsUnique):
            sms._index = isms
        cluster.smsList = smsUnique
        clusters = [cluster]

    for cluster in clusters:
        cluster.txnames = txnames
    return clusters

def groupSMS(smsList, dataset):
    """
    Group SMS into groups where the average SMS
    identical to all the SMS in group.
    The groups contain all SMS which share the same mass,width and upper limit
    and can be replaced by their average SMS when building clusters.

    :parameter smsList: list of all SMS to be grouped
    :parameter dataset: Dataset object to be used when computing distances in upper limit space

    :returns: a list of AverageSMS objects
              which represents a group of SMS with same mass, width and upper limit.
    """


    #  First make sure all SMS contain their upper limits
    for sms in smsList:
        if not hasattr(sms, '._upperLimit'):
            sms._upperLimit = dataset.getUpperLimitFor(sms, txnames=sms.txname)
        if sms._upperLimit is None:
            raise SModelSError("Trying to cluster SMS (id = %i, txname = %s) outside the grid for dataset %s." 
                                %(sms.smsID,sms.txname,dataset.longStr()))

    # Group SMS if they have the same UL
    # and give the same average SMS (same BSM attributes)
    avgSMSList = []
    for iA, smsA in enumerate(smsList):
        avgSMS = AverageSMS([smsA])
        avgSMS._upperLimit = smsA._upperLimit
        for iB, smsB in enumerate(smsList):
            if iB <= iA:
                continue
            if smsA._upperLimit != smsB._upperLimit:
                continue
            if avgSMS != smsB:
                continue
            avgSMS.smsList.append(smsB)
            avgSMS.weight += smsB.weight
        if avgSMS not in avgSMSList:
            avgSMSList.append(avgSMS)

    # Make sure each SMS belongs to a average SMS:
    for sms in smsList:
        nclusters = sum([avgSMS.contains(sms) for avgSMS in avgSMSList])
        if nclusters != 1:
            raise SModelSError("Error computing average SMS. SMS %s belongs to %i average SMS."
                               % (str(sms), nclusters))
    return avgSMSList

def clusterTo(centroids,smsList,dataset,maxDist):
    """
    Assign a SMS from smsList to one of the centroids
    """
    
    dMatrix = np.full((len(smsList),len(centroids)),fill_value=maxDist)
    for isms,sms in enumerate(smsList):
        for ic,c in enumerate(centroids):
            dMatrix[isms,ic] = relativeDistance(c,sms,dataset)

    clusters = [[] for ic in range(len(centroids))]
    notClustered = []
    for isms,sms in enumerate(smsList):
        ic = np.argmin(dMatrix[isms])
        d = dMatrix[isms,ic]
        if d < maxDist:
            clusters[ic].append(isms)
        else:
            notClustered.append(isms)
    
    smsArray = np.array(smsList)
    clusterObjs = []
    for indexList in clusters:
        smsCluster = SMSCluster(smsArray[indexList].tolist())
        # Check if the cluster is valid:
        is_valid = smsCluster.isValid(dataset)        
        if is_valid:
            clusterObjs.append(smsCluster)
        else:
            # If cluster is not valid, keep only the first SMS of the
            # list and move all the other SMS to the notClustered list
            # (note that a cluster with a single SMS is always valid)
            logger.debug('AverageSMS in SMSCluster has no valid UL, splitting cluster.')
            smsCluster_single = SMSCluster(smsArray[indexList[:1]].tolist())
            clusterObjs.append(smsCluster_single)
            notClustered += indexList[1:]

    clusterObjs = sorted(clusterObjs, key = lambda c: sorted([sms.smsID for sms in c.smsList]))

    # Sort not clustered by largest distance first
    notClustered = sorted(notClustered, key = lambda isms: min(dMatrix[isms,:]),reverse=True)
    notClustered = np.array(smsList)[notClustered].tolist()

    return clusterObjs,notClustered   

def kmeansCluster(initialCentroids,sortedSMSList,dataset,maxDist):

    # First cluster around the initial centroids
    clusters,_ = clusterTo(initialCentroids,sortedSMSList,dataset,maxDist)
    
    stop = False
    # Now iterate until the clusters converge
    # (at this stage do not impose a distance limit for clustering)
    while (not stop):
        # Compute new centroids:
        centroids = [cluster.averageSMS for cluster in clusters]
        # Compute new clusters
        newClusters,_ = clusterTo(centroids,sortedSMSList,dataset,maxDist=float('inf'))
        # Stop when clusters no longer change
        if all(clusters[ic] == c for ic,c in enumerate(newClusters)):
            stop = True
        else:
            clusters = newClusters[:]

    # Using the converged centroids remove the SMS with d(centroid,sms) > maxDist
    newClusters,notClustered = clusterTo(centroids,sortedSMSList,dataset,maxDist=maxDist)
    return  newClusters,notClustered

def doCluster(smsList, dataset, maxDist,nmax=100):
    """
    Cluster algorithm to cluster SMS using a modified K-Means method.

    :parameter smsList: list of all SMS to be clustered
    :parameter dataset: Dataset object to be used when computing distances in upper limit space
    :parameter maxDist: maximum distance for clustering two SMS
    :parameter nmax: maximum number of iterations

    :returns: a list of SMSCluster objects containing the SMS
              belonging to the cluster
    """

    # Get average SMS:
    averageSMSList = groupSMS(smsList, dataset)

    # Index average SMS:
    sortedSMSList = sorted(averageSMSList, key=lambda sms: sms._upperLimit)
    for isms, sms in enumerate(sortedSMSList):
        sms._index = isms


    # Choose initial centroids as the first SMS in a chain of
    # SMS all with dist < maxDist:
    centroids = [sortedSMSList[0]]
    for sms in sortedSMSList:
        if relativeDistance(sms,centroids[-1],dataset) > maxDist:
            centroids.append(sms)

    clusters,notClustered = kmeansCluster(centroids,sortedSMSList,dataset,maxDist)

    niter = 1
    while (len(notClustered) !=0) and (niter < nmax):
        centroids = [c.averageSMS for c in clusters] + [notClustered[0]]
        clusters,notClustered = kmeansCluster(centroids,sortedSMSList,dataset,maxDist)
        niter += 1

    if (niter == nmax) and len(notClustered) != 0:
        logger.warning("SMSCluster failed, using unclustered topologies")
        clusters = [SMSCluster([sms]) for sms in sortedSMSList]

    # Replace average SMS by the original SMS:
    for cluster in clusters:
        originalSMS = []
        for avgSMS in cluster.smsList[:]:
            originalSMS += avgSMS.smsList[:]
        cluster.smsList = originalSMS[:]        
    
    # Finally sort clusters by total xsection and length
    clusters = sorted(clusters, key = lambda c: (c.getTotalXSec(),len(c)),reverse=True)

    return clusters
