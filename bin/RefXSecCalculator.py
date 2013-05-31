#!/usr/bin/env python

## https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections8TeVstopsbottom
## https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections

import set_path
import ROOT, os, tempfile
import Theory.XSecComputer as XSEC
import Theory.SLHATools as SLHA
import Experiment.SMSResults as SMS
import Tools.PhysicsUnits as UNIT

topo = 'T2'
nevts = 10000
pidmom = [1000001, 1000002, 1000003, 1000004, 2000001, 2000002, 2000003, 2000004] #lightsquarks
#pidmom = [1000005, 2000005] #sbottoms
#pidmom = [1000006, 2000006] #stops
#pidmom = [1000021] #gluino
#pidmom = [1000024, 1000023] #chi1+, chi20
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
   if not SMS.exists(ana, topo):
      print "No limits found for %s in %s!!!" %(topo, ana)
      continue
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

print x_low, x_up
mlsp = 0.
i = 0.
step = 25.
m = 0.


rootfile = ROOT.TFile('../data/%s_%devts.root' %(topo,nevts) , 'recreate')
txtfile7 = open('../data/%s7TeV_%devts.txt' %(topo,nevts),'w')
txtfile8 = open('../data/%s8TeV_%devts.txt' %(topo,nevts),'w')
txtfile7.write('mass%s      RefXSec[fb]\n' %str(pidmom))
txtfile8.write('mass%s      RefXSec[fb]\n' %str(pidmom))

hist7 = ROOT.TH1F('hist7', 'refxsec7', int((x_up-x_low)/step) , x_low, x_up)
hist8 = ROOT.TH1F('hist8', 'refxsec8', int((x_up-x_low)/step) , x_low, x_up)

Tmp = tempfile.mkdtemp()

while m < x_up-(step/2.):
   m = (x_low+(step/2.))+step*i
   i += 1
#   print m
#   masses = {pidlsp: mlsp, pidmom: m}
   masses = {pidlsp: mlsp}
   for pid in pidmom:
      masses[pid] = m
   slhafile = SLHA.createSLHAFile(topo, masses)
   print "[RefXSecCalculator] compute for",masses
   dic = XSEC.compute(nevts, slhafile, datadir = Tmp)
#   print dic.crossSections()
   xs7 = dic.getSumOfCrossSections(pidmom, sqrts=7)
   xs8 = dic.getSumOfCrossSections(pidmom, sqrts=8)
#   print xs7, xs8
   if xs7:
      hist7.Fill(m, xs7)
      txtfile7.write('%f      %f\n' %(m, xs7))
      #txtfile7.flush()
   if xs8:
      hist8.Fill(m, xs8)
      txtfile8.write('%f      %f\n' %(m, xs8))
      #txtfile8.flush()
   os.unlink(slhafile)

XSEC.clean(Tmp)

txtfile7.close()
txtfile8.close()

hist7.Write()
hist8.Write()
