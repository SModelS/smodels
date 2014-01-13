from Experiment import LimitGetter
import crossSection
import numpy as np
from scipy import stats
import logging,copy
from Tools.PhysicsUnits import addunit,rmvunit
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class Cluster(object):
    """Main class to store the relevant information about a cluster of elements and to manipulate it"""
    
    def __init__(self):
        self.elements = []
        
    def getTotalXSec(self):
        """Returns the sum over the cross-sections of all elements belonging to the cluster""" 
               
        totxsec = crossSection.XSectionList()
        for el in self.elements: totxsec.combineWith(el.weight)
        return totxsec
    
    def getAvgMass(self,method='harmonic'):
        """Returns the average mass of all elements belonging to the cluster using the defined method.
         :param method: the method employed: "harmonic" = harmonic means, "mean" = algebaric (standard) mean
         :returns: the average mass """
          
        massList = [el.getMasses() for el in self.elements]
        return massAvg(massList)
    

def massPosition(mass,Analysis,nounit=True):
    """ gives the mass position in upper limit space, using the analysis experimental limit data.
    If nounit=True, the result is given as number assuming fb units """

    xmass = LimitGetter.GetPlotLimit(mass,Analysis,complain=False)
    if type(xmass) != type(addunit(1.,'pb')): return None
    if nounit: xmass = rmvunit(xmass,'fb')
    return xmass


def massDistance(xmass1,xmass2):
    """ Definition of distance between two mass positions."""

    if xmass1 is None or xmass2 is None: return None
    distance = 2.*abs(xmass1-xmass2)/(xmass1+xmass2)
    if distance < 0.: return None #Skip masses without an upper limit
    return distance


def massAvg(massList,method='harmonic'):
    """Computes the average mass of massList, according to method (harmonic or mean).
    If massList contains a zero mass, switch method to mean"""
        
    flatList = [rmvunit(mass,'GeV') for mass in np.hstack(massList)]
    if method == 'harmonic' and 0. in flatList: method = 'mean'
    
    for mass in massList:
        if len(mass) != len(massList[0]) or len(mass[0]) != len(massList[0]) or len(mass[1]) != len(massList[1]):  
            logger.error('[getAvgMass]: mass shapes do not match in mass list:\n'+str(massList))
            return False
    
    avgmass = copy.deepcopy(massList[0])
    for ib,branch in enumerate(massList[0]):
        for ival,v in enumerate(branch):
            vals = [rmvunit(mass[ib][ival],'GeV') for mass in massList]
            if method == 'mean': avg = np.mean(vals)
            elif method == 'harmonic': avg = stats.hmean(vals)
            avgmass[ib][ival] = addunit(float(avg),'GeV')        
    return avgmass    
    
    
def getGoodMasses(elements,massPos,maxDist):
    """ Get the list of good masses appearing elements according to the Analysis distance.
    Return a dictionary with keys = old mass, values = equivalent good mass"""
    goodMasses = {}
    for element in elements:
        mass = element.getMasses()
        if mass[0] == mass[1] and not massPos[mass]: continue   #skip masses out of bounds
        elif mass[0] == mass[1]: goodMasses[mass] = mass
        else:
            mass1 = [mass[0],mass[0]]
            mass2 = [mass[1],mass[1]]
            if massDistance(massPos[mass1],massPos[mass2]) < maxDist: goodMasses[mass] = massAvg([mass1,mass2])
    return goodMasses    


def groupAll(elements):
    """Creates a single cluster containing all the elements"""
    
    cluster = Cluster()
    cluster.elements = elements
    return cluster


def clusterElements(elements,Analysis,maxDist,keepMassInfo=False):
    """ Cluster the original elements according to their mass distance and return the list of clusters.
        If keepMassInfo, saves the original masses and their cluster value in massDict """
        
#First compute all relevant mass positions:
    massPos = {}
    for element in elements:
        mass = element.getMasses()
        if mass[0] == mass[1]: massPos[mass] = massPosition(mass,Analysis)
        else:
            massPos[[mass[0],mass[0]]] = massPosition([mass[0],mass[0]],Analysis)  #Compute positions for asymmetric masses
            massPos[[mass[1],mass[1]]] = massPosition([mass[1],mass[1]],Analysis)
       
#Get elements with good masses and replace their masses by their 'good' value:    
    goodMasses = getGoodMasses(elements,massPos,maxDist)
    for mass in goodMasses:
        if not mass in massPos: massPos[mass] = massPosition(mass,Analysis)  #Include possible new masses
    goodElements = []
    for el in elements:
        if el.getMasses() in goodMasses:
            element = el.copy()
            element.setMasses(goodMasses[el.getMasses()])
            goodElements.append(element)             
       
#Cluster elements by their mass:
    clusters = doCluster(goodElements,massPos,maxDist)
    return clusters



def doCluster(elements,massPos,maxDist):
    """Cluster algorithm to cluster elements
    :returns: a list of indexes with the clustered masses according to the order\
    in elements
    """

    #Store all distances     
    MD = [] #Mass distance matrix in upper limit space
    massList = [el.getMasses() for el in elements]
    for i,mass1 in massList:
        MD.append([])
        for j,mass2 in massList:
            if j == i: MD[i].append(0.)
            elif j < i: MD[i].append(MD[j][i])
            else:                
                MD[i].append(massDistance(massPos[mass1],massPos[mass2]))  
                                   
    #Begin clustering
    ClusterList = []
    for imass,mass1 in enumerate(massList):
        cluster = set([])
        for jmass,mass2 in enumerate(massList):
            if MD[imass][jmass] == None: continue
            if MD[imass][jmass] <= maxDist: cluster.add(jmass)
        if not cluster in ClusterList: ClusterList.append(cluster)   #Zero level (maximal) clusters


    FinalCluster = []
    newClusters = [0]
    while len(newClusters) > 0:
        newClusters = []
        for cluster in ClusterList:
            split = False
            if len(cluster) > 1:
                clusterMass = massAvg([massList[imass] for imass in cluster])          
                avgPos = massPosition(clusterMass)
                distAvg = max([massDistance(massPos[imass],avgPos) for imass in cluster])
            for imass in cluster:
                clusterDist = ClusterDist(set([imass]),cluster,MD)
                if  clusterDist == None or clusterDist > maxDist or distAvg > maxDist:    #If object or cluster average falls outside the cluster, remove object
                    newcluster = copy.deepcopy(cluster)
                    newcluster.remove(imass)
                    split = True
                    if not newcluster in newClusters: newClusters.append(newcluster)

            if not split and not cluster in FinalCluster: FinalCluster.append(cluster)

        ClusterList = newClusters
        if len(ClusterList) > 500:  #Check for oversized list of cluster (too time consuming)
            logger.warning("[clusterElements] Cluster failed, using unclustered masses")
            FinalCluster = []  
            ClusterList = []
            

    FinalCluster = FinalCluster + ClusterList
    #Add clusters of individual masses (just to be safe)
    for imass,mass in enumerate(massList): FinalCluster.append(set([imass]))

    #Clean up clusters (remove redundant clusters)
    i = 0
    for i,clusterA in enumerate(FinalCluster):
        for j,clusterB in enumerate(FinalCluster):
            if i != j and clusterB.issubset(clusterA): FinalCluster[j] = set([])
    while FinalCluster.count(set([])) > 0: FinalCluster.remove(set([]))

    clusterList = []
    for icluster in FinalCluster:
        cluster = Cluster()
        for iel in icluster: cluster.elements.append(elements[iel])
        clusterList.append(cluster)
    
    return clusterList



def ClusterDist(cluster1,cluster2,massDist):
    """Returns the distance between two clusters = maximum distance between two elements of each cluster.
    massDist = square matrix of distances."""
    
    d = 0.
    if type(cluster1) != type(set()) or type(cluster2) != type(set()):
        logger.warning("[ClusterDist]: unknown format input")
        return False

    for ic in cluster1:
        for jc in cluster2:
            if massDist[ic][jc] == None: return None
            d = max(d,massDist[ic][jc])
            return d



