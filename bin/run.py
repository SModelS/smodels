#!/usr/bin/python

import sys, argparse
sys.path.append ( "../" )
import SMSglobals, SMSmethods, TxNames, LHEReader
import LHEReader
import SMSAnalysisFactory
import SMSDecomp

SMSglobals.initglob()
print "loading analyses ...."
SMSglobals.ListOfAnalyses = SMSAnalysisFactory.load()
print "done loading analyses",len(SMSglobals.ListOfAnalyses)

# n=1000
# slhafile="AndreSLHA/andrePT4.slha"
lhefile="../regression/T1_1.lhe"
reader=LHEReader.LHEReader ( lhefile )
Event=reader.event()
print Event
SMSTop=SMSmethods.getEventTop ( Event.particles )
tx=TxNames.getTx(SMSTop[0].ElList[0])
print "SMSTop=",SMSTop,"tx=",tx
