#!/usr/bin/env python

## https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections8TeVstopsbottom

import set_path
import ROOT, os
import Theory.SMSXSecs as XSEC
import Theory.SLHATools as SLHA
import Experiment.SMSResults as SMS
import Tools.PhysicsUnits as UNIT

topo = 'T2bb'
nevts = 100
pidmom = 1000005 #sbottom
pidlsp = 1000022 #lsp

SMS.considerRuns(["8TeV", "2012", "2011"])
analyses = SMS.getAnalyses(topo)
#print analyses

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


rootfile = ROOT.TFile('../data/%s.root' %topo, 'recreate')
txtfile7 = open('../data/%s7TeV.txt' %topo,'w')
txtfile8 = open('../data/%s8TeV.txt' %topo,'w')
txtfile7.write('mass[%s]      RefXSec[fb]\n' %pidmom)
txtfile8.write('mass[%s]      RefXSec[fb]\n' %pidmom)

hist7 = ROOT.TH1F('hist7', 'refxsec7', int((x_up-x_low)/step) , x_low, x_up)
hist8 = ROOT.TH1F('hist8', 'refxsec8', int((x_up-x_low)/step) , x_low, x_up)

while m < x_up-12.5:
#while i < 2:
   m = (x_low+12.5)+step*i
   i += 1
#   print m
   masses = {pidlsp: mlsp, pidmom: m}
   slhafile = SLHA.createSLHAFile(topo, masses)
   dic = XSEC.computeXSecs(nevts, slhafile)
#   print dic
   if dic['Xsecdic']['7 TeV (NLO)'].has_key((-pidmom,pidmom)):
      xs7 = UNIT.rmvunit(dic['Xsecdic']['7 TeV (NLO)'][(-pidmom,pidmom)], 'fb')
#      print xs7
      hist7.Fill(masses[pidmom], xs7)
      txtfile7.write('%f      %f\n' %(masses[pidmom], xs7))
   if dic['Xsecdic']['8 TeV (NLO)'].has_key((-pidmom,pidmom)):
      xs8 = UNIT.rmvunit(dic['Xsecdic']['8 TeV (NLO)'][(-pidmom,pidmom)], 'fb')
#      print xs8
      hist8.Fill(masses[pidmom], xs8)
      txtfile8.write('%f      %f\n' %(masses[pidmom], xs8))

txtfile7.close()
txtfile8.close()

hist7.Write()
hist8.Write()
