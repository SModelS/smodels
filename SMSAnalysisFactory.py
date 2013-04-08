#!/usr/bin/env python

def load():
    """ This method creates the analysis objects from the info given in the SMS database.
        These EAnalysis objects are registered at SMSglobals.ListOfAnalyses """
    from SMSmethods import EAnalysis
    import SMSResults, SMSglobals

    analysis="alphaT"
    Tx="T1tttt"

    Analysis = EAnalysis()
    Analysis.sqrts=SMSResults.getSqrts( analysis )
    Analysis.label = analysis
    Analysis.masscomp = 0.2
    Analysis.run = SMSResults.getRun ( analysis ) ##  "2012"
    constraint=SMSResults.getConstraints ( analysis )[Tx]
    constrarray=eval(constraint)
    #Global Topology:
    for branch in [0,1]:
      Analysis.Top.B[branch].vertnumb = len(constrarray[branch])+1    #Number of vertices of branch
      vertparts1 = [ len(x) for x in constrarray[branch] ]
      vertparts1.append(0)
      Analysis.Top.B[branch].vertparts = vertparts1
      print vertparts1
    
    cond=SMSResults.getConditions ( analysis)[Tx]
    andreconstraint=constraint.replace("'","").replace("[[[","[[").replace("]]]","]]").replace(" ","")
    print "constraint=",constrarray,andreconstraint
    Analysis.results={ andreconstraint: cond }
    analyses=[ x for x in SMSResults.getAnalyses ( Tx ) if SMSResults.getConditions ( x ).has_key(Tx) and SMSResults.getConditions ( x )[Tx] == cond ]
    print "analyses= ",analyses
    Analysis.plots = { andreconstraint: [ Tx, analyses ] }
    
#Add analysis to list of analyses:
    SMSglobals.ListOfAnalyses.append(Analysis)

# Build list of elements from constraints and conditions with zero weights
# to be computed later with theoretical weights
    for Analy in SMSglobals.ListOfAnalyses:
        Analy.GenerateElements()
        Analy.GetPlots()
    
if __name__ == "__main__":
  import SMSglobals
  load()
  print "List of analyses: "
  for ana in SMSglobals.ListOfAnalyses:
    print ana # .label,ana.sqrts
