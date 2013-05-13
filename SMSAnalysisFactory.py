#!/usr/bin/env python

def getRealTopo ( Tx ):
  """ T3w025 -> T3w, etc """
  ret=Tx
  ret.replace("050","").replace("x1C180","").replace("025","")
  if ret.find("x")>-1: ret=ret[:ret.find("x")]
  return ret

def getArray ( constraint ):
  """ a helper method that is intended to be used to make it possible
      to extract the number of vertices, branches, and insertions from
      the constraint string. It maps e.g.
      2*([[['L'],['L']],[['L'],['nu']]] + [[['L'],['L']],[['nu'],['L']]])
      to
      [[['L'],['L']],[['L'],['nu']]] """
  # print "get array for",constraint
  c=constraint.replace(" ","")
  c=c.strip("1234567890.*%()")
  if c.find("+[")>5:
    c=c[:c.find("+[")]
  if c.find("-[")>5:
    c=c[:c.find("-[")]
  ret=eval(c)
  return ret

def load():
  """ This method creates the analysis objects from the info given in the 
      SMS database.
      These EAnalysis objects are registered at SMSglobals.ListOfAnalyses """
  from SMSmethods import EAnalysis
  import SMSResults, SMSglobals

  debug=False
  if debug:
    SMSResults.verbosity ( "info" )
  else:
    SMSResults.verbosity ( "error" )
  #anas=[ "alphaT", "RA48TeV", "RA2b8TeV", "Weakinos8TeV", "alphaT8TeV", "ATLAS_CONF_2013-024" ]
  anas=SMSResults.getAllResults().keys()
  for analysis in anas:
    if debug:
      print 
      print "Building analysis",analysis
    for Tx in SMSResults.getTopologies(analysis):
      if debug: print Tx,
      Analysis = EAnalysis()
      Analysis.sqrts=SMSResults.getSqrts( analysis )
      stopo=getRealTopo ( Tx )
      Analysis.label = analysis+":"+stopo
      Analysis.masscomp = 0.2
      Analysis.run = SMSResults.getRun ( analysis ) ##  "2012"
      constraint=SMSResults.getConstraints ( analysis, topo=stopo )
      if not constraint: 
        if debug:
          print "dont have a constraint for",analysis,Tx,"(",stopo,")"
        continue
      #print "constraint >>%s<<" % constraint
      constrarray=getArray ( constraint )
      #print "array >>%s<<" % constrarray
      #Global Topology:
      for branch in [0,1]:
        Analysis.Top.B[branch].vertnumb = len(constrarray[branch])+1    #Number of vertices of branch
        vertparts1 = [ len(x) for x in constrarray[branch] ]
        vertparts1.append(0)
        Analysis.Top.B[branch].vertparts = vertparts1
        #print vertparts1
      
      cond=SMSResults.getConditions ( analysis, topo=stopo )
      if cond==None: cond=""
      ## andreconstraint=constraint.replace("'","").replace("[[[","[[").replace("]]]","]]").replace(" ","")
      andreconstraint=constraint.replace("'","").replace(" ","")
      #print "constraint=",constrarray,andreconstraint
      andrecond=cond.replace("'","").replace(" ","")
      #print "cond=",cond
      Analysis.results={ andreconstraint: andrecond }
      analyses=[ x for x in SMSResults.getAnalyses ( stopo ) if SMSResults.getConditions ( x ).has_key(stopo) and SMSResults.getConditions ( x )[stopo] == cond ]
      #print "analyses= ",analyses
      Analysis.plots = { andreconstraint: [ stopo, analyses ] }

      ul=SMSResults.getUpperLimit ( analysis, Tx, mx=400., my=100., interpolate=True )
      print "[SMSAnalysisFactory] test for %s:%s ul(mx=400,my=100)=%s" % ( analysis, Tx, ul )
      
  #Add analysis to list of analyses:
      SMSglobals.ListOfAnalyses.append(Analysis)
      #print "done."

# Build list of elements from constraints and conditions with zero weights
# to be computed later with theoretical weights
  print "Now generate elements"
  for Analy in SMSglobals.ListOfAnalyses:
    #print "Generate",Analy.label
    Analy.GenerateElements()
    Analy.GetPlots(verbose=debug)
    
if __name__ == "__main__":
  import SMSglobals
  load()
  print "List of analyses/results: "
  for (ct,ana) in enumerate(SMSglobals.ListOfAnalyses):
    print ct,ana.label # .label,ana.sqrts
