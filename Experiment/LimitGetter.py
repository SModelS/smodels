#!/usr/bin/env python

"""
.. module:: LimitGetter
    :synopsis: missing

.. moduleauthor:: missing <email@example.com>

"""

def limit(analysis, addTheoryPrediction=True):
  """ the next generation limit retrieval function, should get all information
      from the analysis object """
  import SMSResults
  from Tools.PhysicsUnits import rmvunit
  run=analysis.run
  sqrts=rmvunit(analysis.sqrts,"TeV")
  ## lead=analysis.Top.leadingElement()
  ret=[]
  ##print "[LimitGetter.py] get limit for",analysis,"run=",run
  for (constraint,condition) in analysis.results.items():
    if addTheoryPrediction:
      if not analysis.evaluateResult() or len(analysis.ResultList) == 0: continue
      theoRes=analysis.ResultList[0]
    ##print"[LimitGetter.py] theoRes=",theoRes
    Tx=analysis.plots[constraint][0]
    for ana in analysis.plots[constraint][1]:
      for (index,element) in enumerate(analysis.Top.elements()):
        for (mi,masses1) in enumerate(element.B[0].masses):
        ## masses1=element.B[0].masses[0] ## ,lead.B[1].masses[0]
          masses2=element.B[1].masses[mi] ## ,lead.B[1].masses[0]
          ul=SMSResults.getSmartUpperLimit(ana,Tx,masses1,masses2)
          tmp={ "ul": ul, "analysis": ana, "Tx": Tx, "m1": masses1, "m2": masses2, "sqrts": sqrts }
          if addTheoryPrediction:
            ## theory=theoRes.predictionFor ( masses1, masses2, sqrts, "NLL", condition )
            theory=theoRes.prediction ( )
            ##print "[LimitGetter.py] %s %s ul=%s theory xsec=%s" % ( Tx, ana, ul, theory )
            #          excluded=None
            #print "theory=",theory,"ul=",ul,"type=",type(theory)
            #if True: # theory!=None and ul!=None:
            excluded=rmvunit(theory,"fb")>rmvunit(ul,"fb")
            tmp["theory"]=theory
          ret.append ( tmp )
  return ret 

    
def GetPlotLimit(inmass,Analysis,complain = False):
  """ Get upper limit on sigma*BR for a specific array of masses from plot
        inmass: array of masses in SModelS graph
        Analysis: SMSDataObjects.EAnalysis"""
  from Tools.PhysicsUnits import rmvunit
  import copy, SMSResults

  massarray = copy.deepcopy(inmass)

#Skip empty mass arrays:
  if len(massarray) < 1: 
    if complain: print "[LimitGetter.py] length of massarray < 1"
    return False

#Make sure the two branches have equal masses:
  if massarray[0] != massarray[1]:
    if complain: print "[LimitGetter.py] masses differ between branches"
    return False

  masslist = [rmvunit(mass,'GeV') for mass in massarray[0]]

#Run label:
  run = Analysis.run
  if run == "": run = None   #If run has not been defined, use latest run

  CMSlabel = Analysis.plots.values()[0][0]   #CMS-type label
  analysis = Analysis.plots.values()[0][1][0]  #analysis name
  
  limit = SMSResults.getSmartUpperLimit(analysis,CMSlabel,masslist,run,debug=complain)
 
  return limit

