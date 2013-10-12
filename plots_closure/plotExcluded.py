#!/usr/bin/env python

import argparse
import os, sys

argparser = argparse.ArgumentParser(description='plots the exclusion region and envelope')
argparser.add_argument( 'input', help='input data file (mother_mass,LSP_mass,R)')
args=argparser.parse_args()

from ROOT import *
import AuxPlot

exc = TGraph()
allow = TGraph()
not_tested = TGraph()
exp_limit = TH2F("","",100,0.,2500.,100,0.,2000.)
exc.SetMarkerStyle(20)
exc.SetMarkerColor(kRed-2)
allow.SetMarkerStyle(20)
allow.SetMarkerColor(kAzure-9)
not_tested.SetMarkerStyle(20)
not_tested.SetMarkerColor(kRed-10)


#Define metadata tags:
tags = ['title','Root file','Out file','Kfactor','Root tag']
#Get metadata:
if not os.path.isfile(args.input):
  print 'Input file not found'
  sys.exit()
infile = open(args.input)
data = infile.read()
pts = data[:data.find('#END')-1].split('\n')
info = data[data.find('#END'):].split('\n')
metadata = {}
if len(info) > 0:
  for line in info:
    for tag in tags:
      if tag in line: metadata[tag] = line.lstrip(tag+' :').rstrip()
Rmax = 1.
if 'Kfactor' in metadata.keys():
  Rmax = Rmax/eval(metadata['Kfactor'])
  metadata['title'] += '*'+metadata['Kfactor']

#Get data:
for pt in pts:
  x,y,res,lim,tot = pt.split()  
#  R = float(eval(tot))/float(eval(lim))
  R = float(eval(res))/float(eval(lim))
  x = eval(x)
  y = eval(y)
  lim = eval(lim)
  exp_limit.Fill(x,y,lim)
  if R < 0.:
    not_tested.SetPoint(not_tested.GetN(),x,y)
  elif R >= Rmax:
    exc.SetPoint(exc.GetN(),x,y)
  elif R < Rmax:
    allow.SetPoint(allow.GetN(),x,y)
  else:
    print 'Unknown R value',R
    sys.exit()
    


base = TMultiGraph()
base.Add(exc,"P")
base.Add(allow,"P")
base.Add(not_tested,"P")

   
#Get experimental curve
if 'Root file' in metadata.keys() and os.path.isfile(metadata['Root file']):
  rootfile = TFile(metadata['Root file'],"read")
  objs =  gDirectory.GetListOfKeys()  
  for ob in objs:
    Tob = ob.ReadObj()
    if type(Tob) == type(TGraph()):
      if not 'Root tag' in metadata.keys() or metadata['Root tag'] in ob.GetName():
        if 'expected' in ob.GetName(): continue
        exp = Tob
        base.Add(exp,"L")

plane = TCanvas("c1", "c1",0,0,800,500)
AuxPlot.Default(plane,"TCanvas")
base.Draw("AP")
AuxPlot.Default(base,"TGraph")
base.GetXaxis().SetTitle("M (GeV)")
base.GetYaxis().SetTitle("m_{LSP} (GeV)")
base.GetYaxis().SetTitleOffset(0.75)
gPad.RedrawAxis()  

leg = TLegend(0.6325287,0.7408994,0.9827586,1)
AuxPlot.Default(leg,"Legend")
leg.AddEntry(exc,"Excluded","P")
leg.AddEntry(allow,"Allowed","P")
leg.AddEntry(not_tested,"Not Tested","P")
#if 'Root file' in metadata.keys(): leg.AddEntry(exp,"Official Exclusion","L")
leg.Draw()


if 'title' in metadata.keys():
  title = metadata['title']  
  tit = TPaveLabel(0.054253,0.8308351,0.5948276,0.9486081,title,"NDC")
  tit.SetBorderSize(4)
  tit.SetFillColor(0)
  tit.SetTextFont(42)
  tit.SetTextSize(0.2727273)
  tit.Draw()

#exp_limit.Draw('COLZ')
#exc.Draw("AP")

if 'Out file' in metadata.keys(): c1.Print(metadata['Out file'])
ans = raw_input("Hit any key to close\n")

  
