#!/usr/bin/python

""" This is a simple tool that I (WW) need to work out with the CMS susy group
what digitized results are not yet published """

import set_path
from experiment import SMSResults, SMSHelpers
from tools.physicsUnits import rmvunit

canas={ "cms": {}, "atlas": {} }


def pas  ( i ):
  tmp=SMSResults.getPAS ( i )
  if tmp: return tmp
  return i
  
All=SMSResults.getDatabaseResults ( run=None, category=None )
for (analysis,topos) in All.items():
  if not SMSResults.isPublic ( analysis ): continue
  if not analysis: continue
  panalysis=analysis.replace("_","-")
  for topo in topos:
    constr=SMSResults.getConstraints ( analysis, topo )
    if not constr: continue
    cats=SMSResults.getCategories ( analysis, topo )
    sqrts=SMSResults.getSqrts ( analysis )
    run=SMSResults.getRun ( analysis )
    # print analysis,topo,constr,cat
    experiment="cms"
    if analysis.find("ATLAS")>-1 or analysis.find("SUSY")>-1: 
      experiment="atlas"
    if cats==None: 
      print "category missing for",analysis
      continue
    for cat in cats.split(","):
      cat=cat.strip()
      if not canas[experiment].has_key ( cat ): canas[experiment][cat]=[]
      if not panalysis in canas[experiment][cat]:
        canas[experiment][cat].append ( panalysis )

for (experiment, tanas) in canas.items():
  for (cat, anas) in tanas.items():
    tmp=pas(anas[0])
    sanas="%s~\cite{%s}" %  ( tmp, tmp )
    for i in anas[1:]: 
      tmp=pas(i)
      sanas+=", %s~\cite{%s}" %  ( tmp, tmp )
    print "%s %s: %s" % ( experiment, cat, sanas )
