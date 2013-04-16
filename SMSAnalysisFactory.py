#!/usr/bin/env python

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

  for analysis in [ "alphaT", "RA48TeV", "RA2b8TeV", "Weakinos8TeV" ]:
  #for analysis in [ "Weakinos8TeV" ]:
    for Tx in SMSResults.getTopologies(analysis):
    #for Tx in [ "TChiChipmSlepL" ]:
      # if Tx=="TChiChipmSlepL": continue
      #print "Create analysis object for",analysis,Tx
      Analysis = EAnalysis()
      Analysis.sqrts=SMSResults.getSqrts( analysis )
      Analysis.label = analysis+":"+Tx
      Analysis.masscomp = 0.2
      Analysis.run = SMSResults.getRun ( analysis ) ##  "2012"
      constraint=SMSResults.getConstraints ( analysis, topo=Tx )
      if not constraint: continue
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
      
      cond=SMSResults.getConditions ( analysis, topo=Tx )
      ## andreconstraint=constraint.replace("'","").replace("[[[","[[").replace("]]]","]]").replace(" ","")
      andreconstraint=constraint.replace("'","").replace(" ","")
      #print "constraint=",constrarray,andreconstraint
      andrecond=cond.replace("'","").replace(" ","")
      #print "cond=",cond
      Analysis.results={ andreconstraint: andrecond }
      analyses=[ x for x in SMSResults.getAnalyses ( Tx ) if SMSResults.getConditions ( x ).has_key(Tx) and SMSResults.getConditions ( x )[Tx] == cond ]
      #print "analyses= ",analyses
      Analysis.plots = { andreconstraint: [ Tx, analyses ] }
      
  #Add analysis to list of analyses:
      SMSglobals.ListOfAnalyses.append(Analysis)
      #print "done."

# Build list of elements from constraints and conditions with zero weights
# to be computed later with theoretical weights
  print "Now generate elements"
  for Analy in SMSglobals.ListOfAnalyses:
    #print "Generate",Analy.label
    Analy.GenerateElements()
    Analy.GetPlots(True)
    
if __name__ == "__main__":
  import SMSglobals
  load()
  print "List of analyses/results: "
  for ana in SMSglobals.ListOfAnalyses:
    print ana.label # .label,ana.sqrts
