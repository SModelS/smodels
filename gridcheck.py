#!/usr/bin/python

import sys
sys.path.append('/home/lessa/Pioneer/plots/')
import AuxPlot
from Theory import NLLXSec
from ROOT import *

energy = 8

hk = TH2F("","",100,0.,10000.,100,0.,10000.)
grid = open('grid_gg.dat','w')

for squarkmass in range(50,10000,100):
  for gluinomass in range(50,3000,100):
#for squarkmass in range(5000,5100,100):
#  for gluinomass in range(1000,1100,100):
  
    res = NLLXSec.getNLLfast(process = "gg", pdf = 'cteq', squarkmass=float(squarkmass), gluinomass=float(gluinomass), Energy = energy, base="./nllfast/", interpolate=True )
    kfactor = res['K_NLL_8TeV']*res['K_NLO_8TeV']
    hk.Fill(float(squarkmass),float(gluinomass),kfactor)
    grid.write(str(float(squarkmass))+' '+str(float(gluinomass))+' '+str(kfactor)+'\n')
    
grid.close()
plane = TCanvas("c1", "c1",0,0,700,500)
AuxPlot.Default(plane,"TCanvas")
plane.cd()
plane.SetRightMargin(0.14700422)
plane.SetTopMargin(0.05296053)
plane.SetBottomMargin(0.16796053)
AuxPlot.set_palette(gStyle,name="blue")

hk.SetStats(kFALSE)
hk.Draw("COL Z")
AuxPlot.Default(hk,"TH2")
hk.GetXaxis().SetTitle("m_sq")
hk.GetYaxis().SetTitle("m_gl")
hk.GetZaxis().SetTitle("k factor" )
hk.GetZaxis().SetRangeUser(0.,10.)
hk.GetYaxis().SetRangeUser(0.,3000.)
gPad.RedrawAxis()    

ans = raw_input("Hit any key to close\n")
