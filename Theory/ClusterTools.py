#!/usr/bin/env python

"""
.. module:: ClusterTools
   :synopsis: methods that deal with clustering masses and finding average masses

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

def DoCluster(objlist,Distfunc,dmin,AvgFunc=None,PosFunc=None):
  """Cluster algorithm (generic for any type of object, as long as the distance function is given)
    :returns: a list of indexes with the clustered objects according to the order\
    in objlist
    If AvgFunc is defined, this function is used to compute the cluster center and check if all objects are close to the center
    If PosFunc is defined, uses it to store obj positions (saves time if computing distances is expensive)
  """
  from Tools.PhysicsUnits import addunit
  import copy
  MD = [[None]*len(objlist)]*len(objlist)  #object distances
  MX = None  #object positions
#Compute distance matrix
  if PosFunc: MX = [PosFunc(obj) for obj in objlist]

  for iob,obj1 in enumerate(objlist):
    for job,obj2 in enumerate(objlist):
      if MX:
        MD[iob][job] = Distfunc(MX[iob],MX[job])
      else:
        MD[iob][job] = Distfunc(obj1,obj2)


#Begin clustering
  ClusterList = []
  for iob,obj1 in enumerate(objlist):
    cluster = set([])
    for job,obj2 in enumerate(objlist):
      if MD[iob][job] == None: continue
      if MD[iob][job] <= dmin: cluster.add(job)
    if not cluster in ClusterList: ClusterList.append(cluster)   #Zero level (maximal) clusters


  FinalCluster = []
  newClusters = [0]
  while len(newClusters) > 0:
    newClusters = []
    for cluster in ClusterList:
      split = False
      if len(cluster) > 1:
        if AvgFunc:       #Optional check to see if the cluster average falls inside the cluster
          obj_cluster = [objlist[ic] for ic in cluster]
          obj_avg = AvgFunc(obj_cluster)
          if PosFunc:
            x_avg = PosFunc(obj_avg)
            DistAvg = max([Distfunc(MX[ic],x_avg) for ic in cluster])
          else:
            DistAvg = max([Distfunc(obj,obj_avg) for obj in obj_cluster])
        else:
          DistAvg = 0.

        for i in cluster:
          ClDist = ClusterDist(set([i]),cluster,MD)
          if  ClDist == None or ClDist > dmin or DistAvg > dmin:    #If object or cluster average falls outside the cluster, remove object
            newcluster = copy.deepcopy(cluster)
            newcluster.remove(i)
            split = True
            if not newcluster in newClusters:
              newClusters.append(newcluster)

      if not split and not cluster in FinalCluster: FinalCluster.append(cluster)

    ClusterList = newClusters
    if len(ClusterList) > 1000 or (AvgFunc and len(ClusterList) > 250):  #Check for oversized list of cluster (too time consuming)
      print "DoCluster: Error clustering. ClusterList >",len(ClusterList)
      return None


  FinalCluster = FinalCluster + ClusterList
#Add clusters of individual masses (just to be safe)
  for ic in range(len(objlist)): FinalCluster.append(set([ic]))

#Clean up clusters (remove redundant clusters)
  i = 0
  for i,clusterA in enumerate(FinalCluster):
    for j,clusterB in enumerate(FinalCluster):
      if i != j and clusterB.issubset(clusterA):
        FinalCluster[j] = set([])
  while FinalCluster.count(set([])) > 0: FinalCluster.remove(set([]))

  return FinalCluster

def GoodMass(mass,Distfunc,dmin):
  """Test if a mass array is "good"
     = have similar branch masses if branch topologies are equal
     = have similar mother and LSP masses if branch topologies are different
     AND has an experimental limit
     If it is, return an equivalent array with equal masses (= mass avg)"""

  if mass[0] == mass[1] and Distfunc(mass,mass) == 0.: return mass
  if len(mass[0]) == len(mass[1]):
    mass1 = [mass[0],mass[0]]
    mass2 = [mass[1],mass[1]]
    MD = Distfunc(mass1,mass2)
    if MD == None or MD > dmin:
      return False
    else:
      return MassAvg([mass1,mass2],"harmonic")
  else:
    mass1 = mass
    mass2 = mass
    mass1[1][0] = mass1[0][0]   #Force mothers and daughters to be equal in each branch
    mass1[1][len(mass1)-1] = mass1[0][len(mass1)-1]
    mass2[0][0] = mass2[1][0]
    mass2[0][len(mass2)-1] = mass2[1][len(mass2)-1]
    MD = Distfunc(mass1,mass2)
    if MD == None or MD > dmin:
      return False
    else:
      return MassAvg([mass1,mass2],"harmonic")

def MassAvg(equivin, method = "harmonic"):
  """For a list of equivalent masses, compute an average mass (or mass array)
     using the defined method.
     :param method: the method employed: "harmonic" = harmonic means, "mean" = algebaric (standard) mean

     :returns: the average mass
  """
  import numpy
  from Tools.PhysicsUnits import rmvunit

  N = len(equivin)
  if N == 0:
    print "MassAvg: Empty array"
    return False
  if N == 1: return equivin[0]

  if type(equivin[0]) != type(list()):
    equivinBr = [equivin]
#In case the input has 2 branches of different sizes, average
#each one individually
  elif len(equivin[0]) == 2 and type(equivin[0][0]) == type(list()):
    if len(equivin[0][0]) != len(equivin[0][1]):
      equivinBr = [[],[]]
      for mass in equivin:
        equivinBr[0].append(mass[0])
        equivinBr[1].append(mass[1])
    else:
      equivinBr = [equivin]

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

def sumweights(wlist):
  """Sum a list of weights"""
  from Tools.PhysicsUnits import addunit
  neweight = {}
  for wk in wlist[0].keys(): neweight[wk]=addunit(0.,'fb')# .update({wk : addunit(0.,'fb')})
  for wk in wlist[0].keys():
    wsum = addunit(0.,'fb')
    for weight in wlist: wsum = wsum + weight[wk]
    # neweight.update({wk : wsum})
    neweight[wk]=wsum
  return neweight


def ClusterDist(cluster1,cluster2,MD):
  """Definition of distance two clusters, MD = square matrix of distances"""
  d = 0.
  if type(cluster1) != type(set()) or type(cluster2) != type(set()):
    print "ClusterDist: unknown format input"
    return False

  for ic in cluster1:
    for jc in cluster2:
      if MD[ic][jc] == None: return None
      d = max(d,MD[ic][jc])
  return d

