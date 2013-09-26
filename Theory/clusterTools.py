from Experiment import LimitGetter

clusterAnalysis = None  #Stores the analysis to be used for the calculation of mass distances (only accessed by massPosition)

def massPosition(mass,nounit=True):
    """ gives the mass position in upper limit space, using the analysis experimental limit data.
    If nounit=True, the result is given as number assuming fb units """
    global clusterAnalysis

    xmass = LimitGetter.GetPlotLimit(mass,clusterAnalysis,complain=False)
    if type(xmass) != type(addunit(1.,'pb')): return None
    if nounit: xmass = rmvunit(xmass,'fb')
    return xmass


def massDistance(mass1,mass2):
    """ definition of distance between two mass arrays. The function is defined so it uses the analysis
    experimental limits to define distance """
#Get mass positions for each mass if input are not numbers:
    if type(mass1) != type(1.) or type(mass2) != type(1.):
#If masses differ by more than 100%, do not define distance
        if abs(mass1[0][0]-mass2[0][0])/(mass1[0][0]+mass2[0][0]) > 0.5: return None
        xmass1 = massPosition(mass1)
        xmass2 = massPosition(mass2)
    else:
        xmass1 = mass1
        xmass2 = mass2

    if xmass1 is None or xmass2 is None: return None
    distance = 2.*abs(xmass1-xmass2)/(xmass1+xmass2)
    if distance < 0.: return None         #Skip masses without an upper limit
    return distance

def getMassPositions(massList):
    """ Computes the positions for all masses appearing in massList"""
    massPos = []
    for mass in massList: massPos.append(massPosition(mass))            
    return massPos

def getMassDistances(massList,massPositions=None):
    """ Computes the mass distance matrix for all masses appearing in massList"""
    
    if massPositions:
        massPos = massPositions
    else:
        massPos = getMassPositions(massList)
    massDist = [[None]*len(massList)]*len(massList)
    
    for iel,el1 in enumerate(massList):
        for jel,el2 in enumerate(massList):                
            massDist[iel][jel] = massDistance(massPos[iel],massPos[jel])
    
    return massDist

def ClusterDist(cluster1,cluster2,massDist):
  """Definition of distance two clusters, massDist = square matrix of distances"""
  d = 0.
  if type(cluster1) != type(set()) or type(cluster2) != type(set()):
    print "ClusterDist: unknown format input"
    return False

  for ic in cluster1:
    for jc in cluster2:
      if massDist[ic][jc] == None: return None
      d = max(d,massDist[ic][jc])
  return d


def doCluster(massList,dmin):
    """Cluster algorithm to cluster the masses in massList
    :returns: a list of indexes with the clustered masses according to the order\
    in massList
    """
        
    massPos = getMassPositions(massList) # Mass positions in upper limit space
    massDist = getMassDistances(massList,massPos)  #Mass distance matrix in upper limit space
    #Begin clustering
    ClusterList = []
    for imass,mass1 in enumerate(massList):
        cluster = set([])
        for jmass,mass2 in enumerate(massList):
            if MD[imass][jmass] == None: continue
            if MD[imass][jmass] <= dmin: cluster.add(jmass)
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
                clusterDist = clusterDist(set([imass]),cluster,massDist)
                if  clusterDist == None or clusterDist > dmin or DistAvg > dmin:    #If object or cluster average falls outside the cluster, remove object
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

    return FinalCluster



def goodMass(mass,dmin):
  """Test if a mass array is "good"
     = have similar branch masses if branch topologies are equal
     = have similar mother and LSP masses if branch topologies are different
     AND has an experimental limit
     If it is, return an equivalent array with equal masses (= mass avg)"""

  if mass[0] == mass[1] and massDistance(mass,mass) == 0.: return mass
  if len(mass[0]) == len(mass[1]):
    mass1 = [mass[0],mass[0]]
    mass2 = [mass[1],mass[1]]
    massDist = massDistance(mass1,mass2)
    if massDist == None or massDist > dmin:
      return False
    else:
      return massAvg([mass1,mass2],"harmonic")
  else:
    mass1 = mass
    mass2 = mass
    mass1[1][0] = mass1[0][0]   #Force mothers and daughters to be equal in each branch
    mass1[1][len(mass1)-1] = mass1[0][len(mass1)-1]
    mass2[0][0] = mass2[1][0]
    mass2[0][len(mass2)-1] = mass2[1][len(mass2)-1]
    massDist = massDistance(mass1,mass2)
    if massDist == None or massDist > dmin:
      return False
    else:
      return massAvg([mass1,mass2],"harmonic")
  
  

  

def massAvg(massList, method = "harmonic"):
  """For a list of equivalent masses, compute an average mass (or mass array)
     using the defined method.
     :param method: the method employed: "harmonic" = harmonic means, "mean" = algebaric (standard) mean

     :returns: the average mass
  """
  import numpy
  from Tools.PhysicsUnits import rmvunit

  N = len(massList)
  if N == 0:
    print "MassAvg: Empty array"
    return False
  if N == 1: return massList[0]

  if type(massList[0]) != type(list()):
    equivinBr = [massList]
#In case the input has 2 branches of different sizes, average
#each one individually
  elif len(massList[0]) == 2 and type(massList[0][0]) == type(list()):
    if len(massList[0][0]) != len(massList[0][1]):
      equivinBr = [[],[]]
      for mass in massList:
        equivinBr[0].append(mass[0])
        equivinBr[1].append(mass[1])
    else:
      equivinBr = [massList]

  massout = []
  for ib in range(len(equivinBr)):
    equivmasses = numpy.array(equivinBr[ib])  #Generate numpy array

#Sanity checks:
    for mass in equivmasses.flat:
      if rmvunit(mass,'GeV') == 0.:
        print "MassAvg: Zero mass!"
        return False
      if rmvunit(mass,'GeV') < 0.:
        print "MassAvg: Negative mass!",mass
        return False

    if method == "mean":
      massavg = equivmasses[0]
    elif method == "harmonic":
      massavg = 1./equivmasses[0]
    else:
      print "MassAvg: Unknown method"
      return False

    for imass in range(1,N):
      mass = equivmasses[imass]
      if mass.shape != massavg.shape:    #Sanity check
        print "MassAvg: Wrong input"
        return False
      if method == "mean":
        massavg = massavg + mass
      elif method == "harmonic":
        massavg = massavg + 1./mass

    if method == "mean":
      massavg = massavg/float(N)
    elif method == "harmonic":
      massavg = float(N)/massavg

    if massavg.shape != equivmasses[0].shape:
      print "MassAvg: Error computing average"
      return False

    massout.append(massavg.tolist())

  if len(massout) == 1:
    return massout[0]
  else:
    return massout
