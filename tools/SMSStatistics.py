#!/usr/bin/python

""" A simple tool that prints out statistics """

import sys
sys.path.append ( "../" )

import SMSResults

SMSResults.ResultsForSqrts ( 0 )

All=SMSResults.getAllResults()

Count={"total":0, "public":0, "haveconstraint": 0, "2012": 0, "2011": 0, "ATLAS8TeV": 0, "8TeV":0, "checked": 0 }

for analysis in All:
  for topo in All[analysis]:
    Count["total"]+=1
    if not SMSResults.isPublic ( analysis ): continue
    Count["public"]+=1
    constr=SMSResults.getConstraints ( analysis, topo )
    if not constr: continue
    Count["haveconstraint"]+=1
    run=SMSResults.getRun ( analysis )
    Count[run]+=1
    checked=SMSResults.getCheckedBy ( analysis )
    if checked: Count["checked"]+=1


print Count["total"],"results total ..."
print "... of which",Count["public"],"are public ..."
print "... of which",Count["haveconstraint"],"have a SModelS description."
print "... of which:"
print Count["2011"]+Count["2012"]+Count["8TeV"],"are CMS results"
print Count["ATLAS8TeV"],"are ATLAS results"
print Count["2011"]+Count["2012"],"are 7 TeV results"
print Count["8TeV"]+Count["ATLAS8TeV"],"are 8 TeV results"
print Count["checked"], "out of these", Count["haveconstraint"],"have been verified so far"
