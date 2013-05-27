#!/usr/bin/env python

import set_path
import ROOT
import Theory.SMSxsec as XSEC
import Theory.SLHATools as SLHA
import Experiment.SMSResults as SMS
import Tools.PhysicsUnits as UNIT

topo = 'T2bb'
nevts = 100

SMS.considerRuns(["8TeV", "2012", "2011"])
analyses = SMS.getAnalyses(topo)
print analyses

x_low = 100000.
y_low = 100000.
x_up = 0.
y_up = 0.

for ana in analyses:
   print ana
   if x_low > UNIT.rmvunit(SMS.getLowX(ana, topo), 'GeV'):
      x_low = UNIT.rmvunit(SMS.getLowX(ana, topo), 'GeV')
   print 'xl',x_low, UNIT.rmvunit(SMS.getLowX(ana, topo), 'GeV')
   if x_up < UNIT.rmvunit(SMS.getUpX(ana, topo), 'GeV'):
      x_up = UNIT.rmvunit(SMS.getUpX(ana, topo), 'GeV')
   print 'xu',x_up, UNIT.rmvunit(SMS.getUpX(ana, topo), 'GeV')
   if y_low > UNIT.rmvunit(SMS.getLowY(ana, topo), 'GeV'):
      y_low = UNIT.rmvunit(SMS.getLowY(ana, topo), 'GeV')
   print 'yl',y_low, UNIT.rmvunit(SMS.getLowY(ana, topo), 'GeV')
   if y_up < UNIT.rmvunit(SMS.getUpY(ana, topo), 'GeV'):
      y_up = UNIT.rmvunit(SMS.getUpY(ana, topo), 'GeV')
   print 'yu',y_up, UNIT.rmvunit(SMS.getUpY(ana, topo), 'GeV')

print x_low, x_up, y_low, y_up


#y_low, y_up unnoetig weil lsp masse fix angenommen (sinnvolle annahme?)
#slha fuer massen erstellen (createSLHAFile())
#slhafile = SLHA.createSLHAFile(topo, masses)
#dic = XSEC.pytinit(nevts, slhafile) # im loop fuer mothermass im bin abstand (25 GeV???)

#xs7 = dic['Xsecdic']['7 TeV (NLO)'][(pdg1,pdg2)] # wie komm ich zu den pdgs???????
#xs8 = dic['Xsecdic']['8 TeV (NLO)'][(pdg1,pdg2)]
