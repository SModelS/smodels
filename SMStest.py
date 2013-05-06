#!/usr/bin/env python

import SMSglobals, sys, SMSmethods, SMSxsec, SMSgetlimit
import SMSanalysesTest as SMSanalyses ## make it short
import SMSTranslation

SMSglobals.initglob()
SMSanalyses.load()

nevts = 10
LHEfile=open("regression/T1_1.lhe" )

SMSTopList = []     
for iev in range(nevts):
    PList = SMSmethods.getNextEvent(LHEfile) 

#Get event topology    
    SMSTopListEv = SMSmethods.getEventTop(PList, {} )
    if SMSTopListEv == None: break
    
#Add event topology to topology list:  
    for TopEv in SMSTopListEv:  
        SMSmethods.AddToList(TopEv,SMSTopList)
    
for i in SMSTopList:
  SMSTranslation.getTx ( i )
  for el in i.B[0].ElList:
    print el
  for el in i.B[1].ElList:
    print el
