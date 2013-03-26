#!/usr/bin/python

import SMSglobals, SMSanalyses, sys, SMSmethods, SMSdraw

#PYTHIA must have MSTP(42)=0 ! no mass smearing (narrow width approximation)
#Initialize global variables:
SMSglobals.initglob()
#Creat analyses list:
SMSanalyses.load()

SMSTranslation = { 
  "T1": [ [ 1, [2], ['jet','jet'] ], [ 1, [2], ['jet','jet'] ] ],
  "T2": [ [ 1, [1], ['jet'] ], [ 1, [1], ['jet'] ] ],
}

nevts = 10000
Sqrts = 7
pytres = {"xsecfb" : 1.}


for (sms,andrecode) in SMSTranslation.items():
  print sms,andrecode

  File = open("regression/%s_1.lhe" % sms,"r")
      
  PList = SMSmethods.getNextEvent(File)
  weight = {"7 TeV" : pytres["xsecfb"]/float(nevts), "8 TeV" : pytres["xsecfb"]/float(nevts)}
  SMSTop = SMSmethods.getEventTop(PList, weight)
  b1=[ SMSTop.B[0].vertnumb,SMSTop.B[0].vertparts,SMSTop.B[0].ElList[0].particles ]
  b2=[ SMSTop.B[0].vertnumb,SMSTop.B[0].vertparts,SMSTop.B[0].ElList[0].particles ]
  SMSdraw.asciidraw ( SMSTop )
  check=( [ b1, b2 ] == andrecode)
  #print check
  #if not check:
  #  print [ b1, b2 ],andrecode
