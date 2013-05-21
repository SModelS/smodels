#!/usr/bin/python

import sys, argparse
sys.path.append ( "../" )
sys.path.append ( "../regression" )
import SMSglobals, SMSmethods, TxNames

import SMSAnalysisFactory
import SMSDecomp

SMSglobals.initglob()
print "loading analyses ...."
SMSglobals.ListOfAnalyses = SMSAnalysisFactory.load()
print "done loading analyses",len(SMSglobals.ListOfAnalyses)

n=1000
#slhafile="AndreSLHA/andrePT4.slha"
lhefile="../regression/T1_1.lhe"
File=open(lhefile )
PList=SMSmethods.getNextEvent(File)
SMSTop=SMSmethods.getEventTop ( PList )
print "SMSTop=",SMSTop
#SMSTopList=SMSDecomp.LHEdecomp ( lhefile, {}, 1000 )
