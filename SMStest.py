#!/usr/bin/env python

import SMSglobals, sys, SMSmethods, SMSxsec, SMSgetlimit
import SMSanalysesTest as SMSanalyses ## make it short

SMSglobals.initglob()
SMSanalyses.load()

#Generate events and compute cross-sections:
nevts = 1000
slhafile = "AndreSLHA/andrePT5.slha"
Wv = SMSxsec.pytinit(nevts,slhafile,rpythia = True, donlo = True)
W = Wv["Wdic"]
lhefile = Wv["lhefile"]
LHEfile = open(lhefile,"r")

SMSTopList = []     
#Read events and get topologies (fills SMSTopList)
for iev in range(nevts):
##Read event    
    PList = SMSmethods.getNextEvent(LHEfile) 
##Get mother PDGs:
    momPDG = tuple(SMSmethods.getMom(PList))
#Get event weight list:
    weight = {}
    
    for k in W.keys(): 
        if W[k].has_key(momPDG):
            weight.update({k : W[k][momPDG]})
        else:
            print "Error getting weight"
            sys.exit()    

#Get event topology    
    SMSTopListEv = SMSmethods.getEventTop(PList, weight)
#    for top in SMSTopListEv: print "xxx: ",top,"<br>"
    
#Add event topology to topology list:  
    for TopEv in SMSTopListEv:  
        SMSmethods.AddToList(TopEv,SMSTopList)
    
for top in SMSTopList: 
  print "yyy: ",top,"<br>"

for i in SMSTopList:
  for el in i.B[0].ElList:
    print el
  #for el in range ( len ( SMSTopList[i].B[0].ElList ) ):
  #  print SMSTopList[i].B[0].ElList[el].particles
