#!/usr/bin/env python

"""
.. module:: TheoryPrediction
    :synopsis: Classes encapsulating the results of the computation of reference cross sections
    
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
    
"""

class ClusterOutput:
  """ a wrapper to store the theory predictions from evaluateCluster """

  def __init__ (self):
    self.result_dic = {} ## dictionary of the results, keys are the constraints, values are  TheoryPredicitions
    self.conditions_dic = {} ## dictionary of the conditions, keys are conditions, values are  TheoryPredicitions
    self.mass = None ## average mass inside the cluster
    self.explimit = None ## experimental limit on production cross section for the average mass
    
  def oldformat(self):
    """ Returns a dictionary with the old output format """
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

  def prediction ( self ):
    """ get the theory prediction FIXME (simplified version with new output. What is missing?)  FIXME m1, m2, etc not used???

      :returns: cross section prediction, with units
    """

    return self.result_dic.values()[0]


class TheoryPrediction:
  """ a wrapper for the result of EAnalysis.evaluteResult,
      make it easier to access the theoretical xsec prediction for 
      a particular EElement """

  def equal ( self, m1, m2 ):
    from Tools.PhysicsUnits import rmvunit
    """ is array of masses m1 equal to array m2? """
    if len(m1)!=len(m2): return False
    for (i,m) in enumerate(m1):
      d=abs(rmvunit(m,"GeV")-rmvunit(m2[i],"GeV"))
      if d>1e-1: return False
    return True

  def findIn ( self, m1, dmass ):
    """ search for array of masses m1 in array of array of masses dmass """
    for m2 in dmass:
      if self.equal ( m1, m2 ): return True
    return False

  def __init__ ( self, data ):
    self.data=data
  # make it behave much like a dictionary
  def __len__ ( self ): return len(self.data)
  def __getitem__ ( self, i ): return self.data[i]
  def items ( self ): return self.data.items()
  def __str__ ( self ): return str(self.data)

  def predictionFor ( self, m1=None, m2=None, sqrts=None, order=None, condition=None ):
    """ get the theory xsec prediction for specific conditions:
        m1: get it for this array of masses 
        m2: specify also second array of masses for other branch.
        sqrts: 7 or 8 
        order: LO or NLO
        condition: the condition as is in the database 

        :returns: reference cross section, with units, None if no result available.
    """
    runs=None
    if sqrts!=None and order!=None:
      runs = [ "%d TeV (%s)" % ( int(sqrts), order ) ]
    if sqrts!=None and order==None:
      runs= [ "%d TeV (NLL)" % ( int(sqrts) ), "%d TeV (LO)" % int (sqrts) ]
    if sqrts==None and order!=None:
      runs= [ "7 TeV (%s)" % ( order), "8 TeV (%s)" % order ]

    # print "runs=",runs
    ## print "[TheoryPrediction] runs=",runs
    ret=None
    count=0
    for d in self.data:
      ## check if condition matches, if given
      if condition!="None" and not condition in d['conditions']: continue
      ## check if masses match, if given
      if m1!=None and not self.findIn ( m1, d['mass'] ): continue
      if m2!=None and not self.findIn ( m2, d['mass'] ): continue
      #print "match!"
      res=d['result']
      for (key,value) in res.items():
        if runs==None or key in runs:
          count+=1
          ret=value
          # print "[TheoryPrediction.py] match:",res
    if count>1:
      print "[TheoryPrediction.py] error: more than one result matches description."
    if count==0:
      print "[TheoryPrediction.py] error: no result matches description m1=",m1,"m2=",m2,"sqrts=",sqrts,"order=",order,"condition=",condition
    return ret


