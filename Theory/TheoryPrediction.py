class TheoryPrediction:
  """ basically a wrapper for the result of EAnalysis.evaluteResult,
      make it easier to access the theoretical xsec prediction for 
      a particular EElement """

  def __init__ ( self, data ):
    self.data=data
  # make it behave much like a dictionary
  def __len__ ( self ): return len(self.data)
  def __getitem__ ( self, i ): return self.data[i]
  def items ( self ): return self.data.items()
  def __str__ ( self ): return str(self.data)

  def predictionFor ( self, m1=None, m2=None, sqrts=None, order=None, condition=None ):
    """ get the theory prediction for specific conditions:
        m1: get it for this array of masses 
        m2: specify also second array of masses for other branch.
        sqrts: 7 or 8 
        order: LO or NLO
        condition: the condition as is in the database """
    runs=None
    if sqrts!=None and order!=None:
      runs = [ "%d TeV (%s)" % ( int(sqrts), order ) ]
    if sqrts!=None and order==None:
      runs= [ "%d TeV (NLL)" % ( int(sqrts) ), "%d TeV (LO)" % int (sqrts) ]
    if sqrts==None and order!=None:
      runs= [ "7 TeV (%s)" % ( order), "8 TeV (%s)" % order ]

    ## print "[TheoryPrediction] runs=",runs
    ret=None
    count=0
    for d in self.data:
      ## check if condition matches, if given
      if condition!="None" and not condition in d['conditions']: continue
      ## check if masses match, if given
      if m1!=None and not m1 in d['mass']: continue
      if m2!=None and not m2 in d['mass']: continue
      res=d['result']
      for (key,value) in res.items():
        if runs==None or key in runs:
          count+=1
          ret=value
          # print "[TheoryPrediction.py] match:",res
    if count>1:
      print "[TheoryPrediction.py] error: more than one result matches description."
    if count==0:
      print "[TheoryPrediction.py] error: no result matches description."
    return ret


