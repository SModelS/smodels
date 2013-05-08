#!/usr/bin/python

""" This is a simple tool that I (WW) need to work out with the CMS susy group what digitized results are not yet published """

import SMSResults

SMSResults.ResultsForSqrts ( 0 )

All=SMSResults.getAllResults()

nones={}
falses={}
trues={}

for analysis in All:
  for topo in All[analysis]:
    if not SMSResults.isPublic ( analysis ): continue
    constr=SMSResults.getConstraints ( analysis, topo )
    if not constr: continue
    run=SMSResults.getRun ( analysis )
    if run=="ATLAS8TeV": continue
    ## ok, now we have the subset of CMS results that are public and SModelS-described
    ## the fun begins here.
    published=SMSResults.hasDataPublished ( analysis )
    if published==None: nones[analysis]=True
    if published==True: trues[analysis]=True
    if published==False: falses[analysis]=True

print len(nones),"analyses are not yet checked:"
for a in nones.keys(): print a,
print "\n",len(trues),"analyses have published data:"
for a in trues.keys(): print a,
print "\n",len(falses),"analyses have not published data:"
for a in falses.keys(): print a,
print
