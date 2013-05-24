#!/usr/bin/python

""" simple tool that is meant to draw lessagraphs, both as ascii plot 
    and via pyfeyn """

import set_path
from Experiment import SMSanalyses
from Theory import LHEReader, SMSmethods
import SMSdraw

#PYTHIA must have MSTP(42)=0 ! no mass smearing (narrow width approximation)
#Initialize global variables:
#Creat analyses list:
SMSanalyses.load()

SMSTranslation = { 
#  "T1": [ [ 1, [2], ['jet','jet'] ], [ 1, [2], ['jet','jet'] ] ],
  "T2": [ [ 1, [1], ['jet'] ], [ 1, [1], ['jet'] ] ],
}

nevts = 10000
Sqrts = 7
pytres = {"xsecfb" : 1.}


for (sms,andrecode) in SMSTranslation.items():
  print sms,andrecode

  reader = LHEReader.LHEReader("../regression/%s_1.lhe" % sms,nevts)
  Event = reader.next()
  PList = Event.particles
  weight = {"7 TeV" : pytres["xsecfb"]/float(nevts), "8 TeV" : pytres["xsecfb"]/float(nevts)}
  SMSTop = SMSmethods.getEventTop(PList, weight)
  b1=[ SMSTop[0].vertnumb[0],SMSTop[0].vertparts[0],SMSTop[0].ElList[0].B[0].particles ]
  b2=[ SMSTop[0].vertnumb[1],SMSTop[0].vertparts[1],SMSTop[0].ElList[0].B[1].particles ]
  print b1, b2
  SMSdraw.asciidraw ( SMSTop[0] )
  SMSdraw.draw ( SMSTop[0] )
  check=( [ b1, b2 ] == andrecode)
