#!/usr/bin/env python

"""
.. module:: TheoryPrediction
    :synopsis: Classes encapsulating the results of the computation of reference cross sections
    
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

class XSecPredictionForCluster:
  """ a wrapper to store the theory predictions from evaluateCluster """

  def __init__ (self):
    self.result_dic = {} ## dictionary of the results, keys are the constraints, values are  TheoryPredicitions
    self.conditions_dic = {} ## dictionary of the conditions, keys are conditions, values are  TheoryPredicitions
    self.mass = None ## average cluster mass
    self.allmasses = None ## array to store original masses inside the cluster (if required by the user)
    self.explimit = None ## experimental limit on production cross section for the average mass
  def __str__ ( self ):
    from Tools.PhysicsUnits import rmvunit
    s=""
    for (key,value) in self.result_dic.items():
      m="None"
      if self.mass and len(self.mass)>0:
        m=str( [ rmvunit( x, "GeV" ) for x in self.mass[0] ] )+" GeV"
      onevalue=value
      sovs= [ "7 TeV (LO)", "8 TeV (LO)", "7 TeV (NLO)", "8 TeV (NLO)" ]
      for sov in sovs:
        if value.has_key ( sov ): onevalue=value[sov]
      s+="sigma_ref=%s {%s} for %s: " % ( onevalue, sov, key)
      s+="m=%s. {exp ul=%s}\n" % ( m, self.explimit )
    return s

  def describe ( self ):
    s="TheoryPrediction for m=%s:\n" % self.mass
    s+="================================="
    s+="experimental ul=%s\n" % self.explimit
    for (key,value) in self.result_dic.items():
      s+=">> %s: %s\n" % (key, value)
    s+"\n"
    for (key,value) in self.conditions_dic.items():
      s+="+>> %s: %s\n" % (key, value)
    return s
    
  def oldformat(self):
    """ Returns a dictionary with the old output format """
    import logging
    log = logging.getLogger(__name__)
    log.warning('method is obsolete')
    resdic = {}
    resdic['mass'] = self.mass
    resdic['result'] = self.result_dic.values()[0]
    resdic['conditions'] = {}
    conds = self.conditions_dic.values()
    for weight in conds[0].keys():
      resdic['conditions'].update({weight : []})
    for cond in conds:
      for weight in cond:
        resdic['conditions'][weight].append(cond[weight])

    return resdic
  
#  def getCondition ( self, order="NLL", sqrts=8 ):
#    """ get the condition that has been met, False if not available"""
#    key="%d TeV (%s)" % ( int(sqrts), order)

#    if not self.conditions_dic.has_key ( key ):
#      return False
#    return self.conditions_dic[key]
  

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
  
