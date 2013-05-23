#!/usr/bin/python

import set_path
from Experiment import TxNames, SMSanalyses
from Theory import SMSmethods, LHEReader

#PYTHIA must have MSTP(42)=0 ! no mass smearing (narrow width approximation)
#Initialize global variables:
#Creat analyses list:
# SMSanalyses.load()

green='\033[0;32m'
red='\033[0;31m'
reset='\033[;0m'

topolist = ['T1','T2','T1tttt', 'T2tt','T3W', 'T5WW', 'TChiWZ', 'T1bbbb', 'T2bb', 'T5WZ', 'T3Wb', 'T3Z', 'T5ZZ', 'T6bbZZ', 'TChiWW', 'TSlepSlep']

def ok ( A,B ):
  if A==B: return "%sok.%s" % ( green, reset )
  return "%sfailed. [%s]%s" % ( red, B, reset )

for topo in topolist:
        File=open("%s_1.lhe" % topo)
        reader = LHEReader.LHEReader("%s_1.lhe" % topo)
        Event = reader.next()
        PList = Event.particles
        SMSmethods.GTop()
        SMSTop = SMSmethods.getEventTop(PList, {})
	res = TxNames.getTx(SMSTop[0].ElList[0])
	print "Checking %7s: %s" % ("["+res+"]",ok(res,topo))
