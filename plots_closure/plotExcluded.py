#!/usr/bin/env python

import argparse
import os, sys

argparser = argparse.ArgumentParser(description='plots the exclusion region and envelope')
argparser.add_argument( 'input', help='input data file (mother_mass,LSP_mass,R)')
args=argparser.parse_args()

from ROOT import *
import AuxPlot

#Define metadata tags:
tags = ['title','Root file','Out file','Kfactor','Root tag']
#Get metadata:
if not os.path.isfile(args.input):
  print 'Input file not found'
  sys.exit()
metadata = AuxPlot.getMetadata(args.input,tags)
Rmax = 1.
if metadata['Kfactor']:
  Rmax = Rmax/eval(metadata['Kfactor'][0])
  metadata['title'][0] += '*'+metadata['Kfactor'][0]

#Get data:
results = AuxPlot.getData(args.input,Rmax)
exc = results['exc']
allow = results['allow']
not_tested = results['not_tested']
#Get exclusion envelope:
exc_curve = AuxPlot.getEnvelope(exc)

#Get root plots
rootPlots = AuxPlot.getRootPlots(metadata)

#Set options
exc.SetMarkerStyle(20)
exc.SetMarkerColor(kRed-2)
allow.SetMarkerStyle(20)
allow.SetMarkerColor(kAzure-9)
not_tested.SetMarkerStyle(20)
not_tested.SetMarkerColor(kRed-10)
exc_curve.SetLineColor(kRed)
exc_curve.SetLineStyle(9)
exc_curve.SetLineWidth(3)
for iplot,plot in enumerate(rootPlots.keys()):
  rootPlots[plot].SetLineStyle(iplot+1)
  rootPlots[plot].SetLineWidth(3)
  rootPlots[plot].SetLineColor(kBlack)  

    
base = TMultiGraph()
base.Add(exc,"P")
base.Add(allow,"P")
base.Add(not_tested,"P")    
base.Add(exc_curve,"L")
for plot in rootPlots.keys(): base.Add(rootPlots[plot],'L')

#Legend
leg = TLegend(0.6325287,0.7408994,0.9827586,1)
AuxPlot.Default(leg,"Legend")
leg.AddEntry(exc,"Excluded","P")
leg.AddEntry(allow,"Allowed","P")
leg.AddEntry(not_tested,"Not Tested","P")
for plot in rootPlots.keys(): leg.AddEntry(rootPlots[plot],plot,'L')

#Canvas
plane = TCanvas("c1", "c1",0,0,800,500)
AuxPlot.Default(plane,"TCanvas")
base.Draw("AP")
AuxPlot.Default(base,"TGraph")
base.GetXaxis().SetTitle("M (GeV)")
base.GetYaxis().SetTitle("m_{LSP} (GeV)")
base.GetYaxis().SetTitleOffset(0.75)
gPad.RedrawAxis()  
leg.Draw()

#Add title
if metadata['title']:
  title = metadata['title'][0]
  tit = TPaveLabel(0.054253,0.8308351,0.5948276,0.9486081,title,"NDC")
  tit.SetBorderSize(4)
  tit.SetFillColor(0)
  tit.SetTextFont(42)
  tit.SetTextSize(0.2727273)
  tit.Draw()


#if metadata['Out file']: c1.Print(metadata['Out file'][0])
ans = raw_input("Hit any key to close\n")

  
