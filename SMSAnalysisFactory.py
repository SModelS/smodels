#!/usr/bin/env python

def load():
    """ This method creates the analysis objects from the info given in the SMS database.
        These EAnalysis objects are registered at SMSglobals.ListOfAnalyses """
    from SMSmethods import EAnalysis
    import SMSResults, SMSglobals

    analysis="alphaT"
    topo="T1tttt"

    Analysis = EAnalysis()
    Analysis.sqrts=SMSResults.getSqrts( analysis )
    Analysis.histname
    Analysis.label = analysis
    Analysis.masscomp = 0.2
    Analysis.run = "2012"
    condition=SMSResults.getConditions ( analysis )[topo]
    condarray=eval(condition)
    print "condition=",condarray
    #Global Topology:
    for branch in [0,1]:
      Analysis.Top.B[branch].vertnumb = len(condarray[branch])+1    #Number of vertices of branch
      vertparts1 = [ len(x) for x in condarray[branch] ]
      vertparts1.append(0)
      Analysis.Top.B[branch].vertparts = vertparts1
      print vertparts1
    
    Analysis.results={ condition: SMSResults.getConstraints ( analysis)[topo] }
    Analysis.plots = { condition: [ topo, [ analysis ] ] }
    
#Add analysis to list of analyses:
    SMSglobals.ListOfAnalyses.append(Analysis)

# Build list of elements from constraints and conditions with zero weights
# to be computed later with theoretical weights
    for Analy in SMSglobals.ListOfAnalyses:
        Analy.GenerateElements() 
    
if __name__ == "__main__":
  import SMSglobals
  load()
  print "List of analyses: "
  for ana in SMSglobals.ListOfAnalyses:
    print ana # .label,ana.sqrts
