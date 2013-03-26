#!/usr/bin/python

import SMSglobals, SMSanalyses, sys, SMSmethods

#PYTHIA must have MSTP(42)=0 ! no mass smearing (narrow width approximation)
#Initialize global variables:
SMSglobals.initglob()
#Creat analyses list:
SMSanalyses.load()

SMSTranslation = { 
  "T2": [ [ 1, [1], ['jet'] ], [ 1, [1], ['jet'] ] ],
  "T1": [ [ 1, [2], ['jet','jet'] ], [ 1, [2], ['jet','jet'] ] ],
  "T3w": [ [ 2, [2], ['jet','jet'] ], [ 1, [2], ['jet','jet'] ] ],
}


nevts = 10000
Sqrts = 7
pytres = {"xsecfb" : 1.}


for (sms,andrecode) in SMSTranslation.items():
  print sms,andrecode

  File = open("regression/%s_1.lhe" % sms,"r")
      
  PList = SMSmethods.getNextEvent(File)
  weight = [pytres["xsecfb"]/float(nevts),pytres["xsecfb"]/float(nevts)]
  SMSmethods.GTop()
  SMSTop = SMSmethods.getEventTop(PList, weight)
  b1=[ SMSTop.B[0].vertnumb,SMSTop.B[0].vertins,SMSTop.B[0].ElList[0].particles ]
  b2=[ SMSTop.B[0].vertnumb,SMSTop.B[0].vertins,SMSTop.B[0].ElList[0].particles ]
  check=( [ b1, b2 ] == andrecode)
  print check
  if not check:
    print [ b1, b2 ],andrecode
