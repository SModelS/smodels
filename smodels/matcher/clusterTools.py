"""
.. module:: clusterTools
   :synopsis: Module holding the SMSCluster class and cluster methods used to combine similar SMS according
      to the analysis.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.theory.theorySMS import TheorySMS
from smodels.theory.particleNode import ParticleNode
from smodels.theory.particle import MultiParticle
from smodels.experiment.datasetObj import DataSet, CombinedDataSet
from smodels.tools.physicsUnits import fb
from smodels.matcher.exceptions import SModelSMatcherError as SModelSError
from smodels.matcher.matcherAuxiliaryFuncs import average
import numpy as np
from smodels.tools.smodelsLogging import logger


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

        otherProperties = [getattr(other, attr) for attr in self.properties]
        selfProperties = [getattr(self, attr) for attr in self.properties]
        comp = (selfProperties > otherProperties) - (otherProperties > selfProperties)

        return comp

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

    def __init__(self, smsList=[], dataset=None, distanceMatrix=None):

        self.smsList = smsList
        self.dataset = dataset
        self.maxInternalDist = 0.
        self._distanceMatrix = distanceMatrix
        # Compute maximal internal distance
        if self.smsList and self._distanceMatrix is not None:
            self.maxInternalDist = max([self.getDistanceTo(sms) for sms in self])

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
        totxsec = sum([sms.weight for sms in self.smsList])
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
        if self.dataset:
            avgSMS._upperLimit = self.dataset.getUpperLimitFor(avgSMS,
                                                              txnames=avgSMS.txname)

        avgSMS._index = None

        return avgSMS

    def copy(self):
        """
        Returns a copy of the index cluster (faster than deepcopy).
        """

        newcluster = SMSCluster(self.smsList[:], self.dataset, self._distanceMatrix)
        newcluster.maxInternalDist = self.maxInternalDist
        return newcluster

    def add(self, smsList):
        """
        Add an SMS or list of SMS.

        :param smsList: TheorySMS object or list of SMS
        """

        if not isinstance(smsList, list):
            smsList = [smsList]
        else:
            smsList = smsList[:]

        for sms in smsList:
            if sms._index in self.indices():
                continue

            self.smsList.append(sms)
            # Update internal distance:
            self.maxInternalDist = max(self.maxInternalDist, self.getDistanceTo(sms))

    def remove(self, smsList):
        """
        Remove an SMS or a list of SMS from the cluster.

        :param smsList: TheorySMS object or list of SMS
        """

        if not isinstance(smsList, list):
            smsList = [smsList]
        else:
            smsList = smsList[:]

        for sms in smsList:
            indices = self.indices()
            if sms._index not in indices:
                continue
            isms = indices.index(sms._index)
            self.smsList.pop(isms)
        # Update internal distance:
        self.maxInternalDist = max([self.getDistanceTo(elB) for elB in self])

    def getDistanceTo(self, sms):
        """
        Return the maximum distance between any of the SMS belonging to the
        cluster and sms.

        :parameter sms: SMS object
        :return: maximum distance (float)
        """

        if not hasattr(sms, '_upperLimit'):
            sms._upperLimit = self.dataset.getUpperLimitFor(sms,
                                                            txnames=sms.txname)
        if sms._upperLimit is None:
            return None

        # Use pre-computed distances for regular (non-averge) SMS
        if sms._index is not None:
            return max([self._distanceMatrix[sms._index, el._index] for el in self])

        dmax = 0.
        for smsSelf in self:
            if not hasattr(smsSelf, '_upperLimit'):
                smsSelf._upperLimit = self.dataset.getUpperLimitFor(smsSelf,
                                                               txnames=smsSelf.txname)
            dmax = max(dmax, relativeDistance(sms, smsSelf, self.dataset))

        return dmax

    def isConsistent(self, maxDist):
        """
        Checks if the cluster is consistent.
        Computes an average SMS in the cluster
        and checks if this average SMS belongs to the cluster
        according to the maximum allowed distance between cluster SMS.

        :return: True/False if the cluster is/is not consistent.
        """

        avgSMS = self.averageSMS
        if avgSMS._upperLimit is None:
            return False

        dmax = self.getDistanceTo(avgSMS)
        if dmax > maxDist:
            return False

        return True


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
    # Sort SMS, so the ones with highest contribution (weight*eff) come first:
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
        distanceMatrix = np.zeros((len(smsUnique), len(smsUnique)))
        cluster = SMSCluster(dataset=dataset, distanceMatrix=distanceMatrix)
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
            raise SModelSError("Trying to cluster SMS outside the grid.")

    # Group SMS if they have the same UL
    # and give the same average SMS (same mass and same width)
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


def doCluster(smsList, dataset, maxDist):
    """
    Cluster algorithm to cluster SMS.

    :parameter smsList: list of all SMS to be clustered
    :parameter dataset: Dataset object to be used when computing distances in upper limit space
    :parameter maxDist: maximum distance for clustering two SMS

    :returns: a list of SMSCluster objects containing the SMS
              belonging to the cluster
    """

    # Get average SMS:
    averageSMSList = groupSMS(smsList, dataset)

    # Index average SMS:
    sortedSMSList = sorted(averageSMSList, key=lambda sms: sms._upperLimit)
    for isms, sms in enumerate(sortedSMSList):
        sms._index = isms

    # Pre-compute all necessary distances:
    distanceMatrix = np.zeros((len(sortedSMSList), len(sortedSMSList)))
    for isms, smsA in enumerate(sortedSMSList):
        for jsms, smsB in enumerate(sortedSMSList):
            if jsms <= isms:
                continue
            distanceMatrix[isms, jsms] = relativeDistance(smsA, smsB, dataset)
    distanceMatrix = distanceMatrix + distanceMatrix.T

    # Start building maximal clusters
    clusterList = []
    for sms in sortedSMSList:
        cluster = SMSCluster([], dataset, distanceMatrix)
        for smsB in sortedSMSList:
            if distanceMatrix[sms._index, smsB._index] <= maxDist:
                cluster.add(smsB)
        if not cluster.smsList:
            continue
        if cluster.averageSMS._upperLimit is None:
            continue
        if cluster not in clusterList:
            clusterList.append(cluster)

    # Split the maximal clusters until all SMS inside each cluster are
    # less than maxDist apart from each other and the cluster average position
    # is less than maxDist apart from all SMS
    finalClusters = []
    while clusterList:
        newClusters = []
        for cluster in clusterList:
            # Check if maximal internal distance is below maxDist
            isConsistent = cluster.isConsistent(maxDist)
            if isConsistent and cluster.maxInternalDist < maxDist:
                if cluster not in finalClusters:
                    finalClusters.append(cluster)

            # Cluster violates maxDist:
            else:
                # Loop over cluster SMS and if the SMS distance
                # falls outside the cluster, remove SMS
                for sms in cluster:
                    if cluster.getDistanceTo(sms) > maxDist or not isConsistent:
                        newcluster = cluster.copy()
                        newcluster.remove(sms)
                        if newcluster.averageSMS._upperLimit is None:
                            continue
                        if newcluster in newClusters:
                            continue
                        newClusters.append(newcluster)

        clusterList = newClusters
        #  Check for oversized list of indexCluster (too time consuming)
        if len(clusterList) > 100:
            logger.warning("SMSCluster failed, using unclustered masses")
            finalClusters = []
            clusterList = []

    #  finalClusters = finalClusters + clusterList
    #  Add clusters of individual masses (just to be safe)
    for sms in sortedSMSList:
        finalClusters.append(SMSCluster([sms], dataset, distanceMatrix))

    #  Clean up clusters (remove redundant clusters)
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

    # Replace average SMS by the original SMS:
    for cluster in finalClusters:
        originalSMS = []
        for avgSMS in cluster.smsList[:]:
            originalSMS += avgSMS.smsList[:]
        cluster.smsList = originalSMS[:]

    return finalClusters
