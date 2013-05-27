#!/usr/bin/env python

import set_path
import ROOT, os
import Theory.SMSxsec as XSEC
import Theory.SLHATools as SLHA
import Experiment.SMSResults as SMS
import Tools.PhysicsUnits as UNIT

topo = 'T2bb'
nevts = 100
pidmom = 1000005 #sbottom
pidlsp = 1000022 #lsp

SMS.considerRuns(["8TeV", "2012", "2011"])
analyses = SMS.getAnalyses(topo)
print analyses

x_low = 100000.
y_low = 100000.
x_up = 0.
y_up = 0.

for ana in analyses:
#   print ana
   if x_low > UNIT.rmvunit(SMS.getLowX(ana, topo), 'GeV'):
      x_low = UNIT.rmvunit(SMS.getLowX(ana, topo), 'GeV')
#   print 'xl',x_low, UNIT.rmvunit(SMS.getLowX(ana, topo), 'GeV')
   if x_up < UNIT.rmvunit(SMS.getUpX(ana, topo), 'GeV'):
      x_up = UNIT.rmvunit(SMS.getUpX(ana, topo), 'GeV')
#   print 'xu',x_up, UNIT.rmvunit(SMS.getUpX(ana, topo), 'GeV')
#   if y_low > UNIT.rmvunit(SMS.getLowY(ana, topo), 'GeV'):
#      y_low = UNIT.rmvunit(SMS.getLowY(ana, topo), 'GeV')
#   print 'yl',y_low, UNIT.rmvunit(SMS.getLowY(ana, topo), 'GeV')
#   if y_up < UNIT.rmvunit(SMS.getUpY(ana, topo), 'GeV'):
#      y_up = UNIT.rmvunit(SMS.getUpY(ana, topo), 'GeV')
#   print 'yu',y_up, UNIT.rmvunit(SMS.getUpY(ana, topo), 'GeV')

#print x_low, x_up, y_low, y_up
print x_low, x_up
mlsp = 0.
i = 0.
step = 25.
m = 0.
#y_low, y_up unnoetig weil lsp masse fix angenommen (sinnvolle annahme?)
#while m < x_up-12.5:
while i < 2:
   m = (x_low+12.5)+step*i
   i += 1
#   print m
   masses = {pidlsp: mlsp, pidmom: m}
   slhafile = SLHA.createSLHAFile(topo, masses)
   wkdir = os.getcwd()
   os.chdir('../')
   print wkdir, os.getcwd()
   dic = XSEC.pytinit(nevts, slhafile)
   os.chdir(wkdir)
   print dic
   if dic['Xsecdic']['7 TeV (NLO)'].has_key((pdgmom,pdgmom)):
      xs7 = dic['Xsecdic']['7 TeV (NLO)'][(pdgmom,pdgmom)]
      print xs7
   xs8 = dic['Xsecdic']['8 TeV (NLO)'][(pdg1,pdg2)]
