#!/usr/bin/env python

"""
.. module:: TheoryPrediction
    :synopsis: Classes encapsulating the results of the computation of reference cross sections
    
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
    
"""

import copy
from ClusterTools import DoCluster, GoodMass, MassAvg
from Experiment import LimitGetter
from SMSDataObjects import CTop



class XSecPredictionForCluster:
  """ a wrapper to store the theory predictions from evaluateCluster """

  def __init__ (self):
    self.result_dic = {} ## dictionary of the results, keys are the constraints, values are  TheoryPredicitions
    self.conditions_dic = {} ## dictionary of the conditions, keys are conditions, values are  TheoryPredicitions
    self.mass = None ## average cluster mass
    self.allmasses = None ## array to store original masses inside the cluster (if required by the user)
    self.explimit = None ## experimental limit on production cross section for the average mass


  def getConditionValues ( self, wlabel = ''):
    """ get the list of condition values for the cross-section label wlabel. Return False if not available """
    conds = self.conditions_dic.values()
    if not conds[0].has_key ( wlabel ):
      return False
    else:
      res = []
      for cond in conds:
        res.append(cond[wlabel])
    return res
  
  def getMaxCondition(self,wlabel = None):
    """ get the maximum condition value, False if wlabel does not exist. If wlabel = None, get the maximum condition
        over all weight labels"""

    if wlabel is None:
      all_conds = []
      conds = self.conditions_dic.values()
      for cond in conds:
        all_conds.extend(cond.values())  #Has all condition values, for all weight labels

      conds = list(set(all_conds))   #Remove repeated entries
    else:
      conds = self.getConditionValues(wlabel)

    if conds is False: return False
    if conds is None or conds == [None]: return None
    if 'N/A' in conds: return 'N/A'
    maxcond = 0.
    for tvalue in conds:
      if type(tvalue) == type(1.): maxcond = max(maxcond,tvalue)
    return maxcond

  def prediction ( self, wlabel = None ):
    """ get the theory prediction for all cross-sections if
    wlabel = None or for the specific cross-section label if wlabel = label. Assumes a single result per analysis.
      :returns: cross section prediction, with units
    """
    if len(self.result_dic.values()) > 1:
      return "[TheoryPrediction]: unknown result_dic format"
    if wlabel is None:
      return self.result_dic.values()[0]
    elif not self.result_dic.values()[0].has_key(wlabel):
      return None
    else:
      return self.result_dic.values()[0][wlabel]
  


class TheoryPrediction:
  """ Holds the theory prediction for a single analysis """
  def __init__ (self):
    self.analysis = None ## holds the corresponding Analysis object
    self.clusterResults = []  #holds a list of XSecPredictionForCluster objects
    self.goodmasses = {} ## a dictionary to store masses and their equivalent goodmass value (if required by the user)


  def computeTheoryPredictions(self,Analysis,SMSTopList,keepMassInfo=False):
    """ Main method for evaluating the theoretical predictions to the Analysis \
    from the decomposition results (SMSTopList). It first generates a topology with all the
    elements from the Analysis and add to them the matching elements from SMSTopList.
    Then it combines equivalent masses using the Cluster tools and calls \
    evaluateCluster to compute the theoretical predictions for each cluster. \
    Produces self.clusterResults.
      :returns: True if successful, None if list is empty.
    """

#Adds the matching elements from SMSTopList to the analysis elements
    resultsTop = copy.deepcopy(Analysis.Top)
    self.analysis = Analysis

    for eventTop in SMSTopList:
#Check if topologies match:
      if not eventTop == resultsTop: continue   
#Loop over (event) element list:
      for eventElement in eventTop.ElList: resultsTop.addEventElement(eventElement,Analysis.sqrts)

    dmin = Analysis.masscomp
  #Create a mass list with all masses appearing in the analysis elements which have similar branch masses:
    Goodmasses = []
    for El in resultsTop.ElList:
      for massweight in El.MassWeightList:
        gmass = GoodMass(massweight.mass,Analysis.MassDist,dmin)
        if keepMassInfo: self.goodmasses[str(massweight.mass)] = gmass
        if gmass:
           massweight.mass = gmass
           if not gmass in Goodmasses: Goodmasses.append(gmass)

  #Cluster masses:
    MCluster = DoCluster(Goodmasses,Analysis.MassDist,dmin,MassAvg,Analysis.MassPosition)

    if MCluster == None or MCluster == False:
      MCluster = []
      for i in range(len(Goodmasses)): MCluster.append(set([i]))
      print "[SMSAnalysis.py] Cluster failed, using unclustered masses"
      
  #Loop over clusters to evaluate constraints and conditions inside each cluster
    self.clusterResults = []  # Clear out results
    for cluster in MCluster:
  #Get masses in cluster
      masscluster = []
      for ic in cluster: masscluster.append(Goodmasses[ic])
  #Create clustered topology from the analysis topology elements belonging to cluster:
      NewTop = CTop(resultsTop,masscluster)

  #Now NewTop contains only elements with a common mass (replaced by the average mass)
  #Evaluate theoretical predictions for analysis inside cluster
      ClusterResult = NewTop.evaluateCluster(Analysis.results)
      ClusterResult.mass = NewTop.clustermass
      ClusterResult.explimit = LimitGetter.GetPlotLimit(ClusterResult.mass,Analysis,complain=False)
      if keepMassInfo: ClusterResult.allmasses = masscluster

  #Check if average mass is inside the cluster (exp. limit for average mass ~ exp. limit for individual masses):
      mavg = NewTop.clustermass
      davg = -1.
      for mass in masscluster:
        davg = max(davg,Analysis.MassDist(mass,mavg))
      if davg == -1. or davg > dmin:
        log = logging.getLogger(__name__)
        log.warning ( "Upper limit for average mass is not similar to upper limit for constituent masses (davg=%s, dmin=%s). Wont cluster." % ( davg, dmin ) )
        continue

      self.clusterResults.append(ClusterResult)

    if not self.clusterResults: return None      
    return True

