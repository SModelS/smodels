#!/usr/bin/python

import SMSglobals, SMSanalyses, sys, SMSmethods, SMSTranslation

#PYTHIA must have MSTP(42)=0 ! no mass smearing (narrow width approximation)
#Initialize global variables:
SMSglobals.initglob()
#Creat analyses list:
SMSanalyses.load()

SMSTranslation = SMSTranslation.SMSTranslation

nevts = 10000
Sqrts = 7
pytres = {"xsecfb" : 1.}


for (sms,andrecode) in SMSTranslation.items():
  print sms,andrecode

  File = open("regression/%s_1.lhe" % sms,"r")
      
  PList = SMSmethods.getNextEvent(File)
  weight = [pytres["xsecfb"]/float(nevts),pytres["xsecfb"]/float(nevts)]
  SMSmethods.GTop()
  SMSTop = SMSmethods.getEventTop(PList, {})
#Sort branches
  Branch1 = SMSTop[0].B[0]
  Branch2 = SMSTop[0].B[1]
  if (Branch1.vertnumb < Branch2.vertnumb) or (Branch1.vertnumb == Branch2.vertnumb and sum(Branch1.vertparts) < sum(Branch2.vertparts)):
        SMSTop[0].B = [Branch2,Branch1]
#read branch infortmation
  b1=[ SMSTop[0].B[0].vertnumb,SMSTop[0].B[0].vertparts,SMSTop[0].B[0].ElList[0].particles ]
  b2=[ SMSTop[0].B[1].vertnumb,SMSTop[0].B[1].vertparts,SMSTop[0].B[1].ElList[0].particles ]
  check=( [ b1, b2 ] == andrecode)
  print check
  if not check:
    print [ b1, b2 ],andrecode
