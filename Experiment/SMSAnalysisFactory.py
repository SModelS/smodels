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

def load( anas = None, topos=None ):
  """ This method creates the analysis objects from the info given in the 
      SMS database, returns ListOfAnalyses. 
      If anas is given as a list, then we create only objects for these analyses
      (the database naming convention is used). 
      If topos is given as a list, then only these topos are considered """
  from Theory.SMSmethods import EAnalysis
  from Experiment import SMSResults
  import types

  ## lets make sure the user can also supply a single topo/ana
  ## without having to code an array
  if type(topos)==types.StringType: topos = [ topos ]
  if type(anas)==types.StringType: anas = [ anas ]
  
  ListOfAnalyses = []

  debug=False
  if debug:
    SMSResults.verbosity ( "info" )
  else:
    SMSResults.verbosity ( "error" )
  #anas=[ "alphaT", "RA48TeV", "RA2b8TeV", "Weakinos8TeV", "alphaT8TeV", "ATLAS_CONF_2013-024" ]
  if anas==None: anas=SMSResults.getAllResults().keys()
  for analysis in anas:
    if debug:
      print 
      print "Building analysis",analysis
    for Tx in SMSResults.getTopologies(analysis):
      if topos!=None and Tx not in topos: continue
      if debug: print Tx,
      # if Tx[:2]!="T2": continue
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
      constrarray=getArray ( constraint )
      #Global Topology:
      Analysis.Top.vertnumb = [ len(constrarray[0])+1, len(constrarray[1])+1 ]
      Analysis.Top.vertparts = [ [ len(x) for x in constrarray[0] ], [ len(x) for x in constrarray[1] ] ]
      Analysis.Top.vertparts[0].append(0)
      Analysis.Top.vertparts[1].append(0)
      #for branch in [0,1]:
        # Analysis.Top.B[branch].vertnumb = len(constrarray[branch])+1    #Number of vertices of branch
        # vertparts1 = [ len(x) for x in constrarray[branch] ]
      #  vertparts1.append(0)
      #  Analysis.Top.B[branch].vertparts = vertparts1
      
      cond=SMSResults.getConditions ( analysis, topo=stopo )
      if cond==None: cond=""
      andreconstraint=constraint.replace("'","").replace(" ","")
      andrecond=cond.replace("'","").replace(" ","")
      Analysis.results={ andreconstraint: andrecond }
      analyses=[ analysis ]
      ## analyses=[ x for x in SMSResults.getAnalyses ( stopo ) if SMSResults.getConditions ( x ).has_key(stopo) and SMSResults.getConditions ( x )[stopo] == cond ]
      Analysis.plots = { andreconstraint: [ stopo, analyses ] }

      # ul=SMSResults.getUpperLimit ( analysis, Tx, mx=400., my=100., interpolate=True )
#      print "[SMSAnalysisFactory] test for %s:%s ul(mx=400,my=100)=%s" % ( analysis, Tx, ul )
      
  #Add analysis to list of analyses:
      ListOfAnalyses.append(Analysis)
      #print "done."

# Build list of elements from constraints and conditions with zero weights
# to be computed later with theoretical weights
#  print "[SMSAnalysisFactory.py] Now generate elements"
  for Analy in ListOfAnalyses:
#    print "[SMSAnalysisFactory.py] Generate element ``%s [%s]'' " % ( Analy.label,Analy.run )
    Analy.GenerateElements()
    Analy.GetPlots(verbose=debug)
    
  return ListOfAnalyses  
    
if __name__ == "__main__":
  load()
  print "List of analyses/results: "
  ListOfAnalyses=load()
  for (ct,ana) in enumerate(ListOfAnalyses):
    print ct,ana.label,ana.Top.vertnumb,ana.Top.vertparts # .label,ana.sqrts
