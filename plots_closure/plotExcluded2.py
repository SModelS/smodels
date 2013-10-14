#!/usr/bin/env python

import argparse
import os, sys

argparser = argparse.ArgumentParser(description='plots the exclusion region and envelope')
argparser.add_argument( 'input', help='input data file (mother_mass,LSP_mass,R)')
args=argparser.parse_args()

from ROOT import *
import AuxPlot

exc = TGraph()
exc_curve = TGraph()
allow = TGraph()
not_tested = TGraph()
exp_limit = TH2F("","",100,0.,2500.,100,0.,2000.)
exc.SetMarkerStyle(20)
exc.SetMarkerColor(kRed-2)
allow.SetMarkerStyle(20)
allow.SetMarkerColor(kAzure-9)
not_tested.SetMarkerStyle(20)
not_tested.SetMarkerColor(kRed-10)
exc_curve.SetLineColor(kRed)
exc_curve.SetLineStyle(9)
exc_curve.SetLineWidth(3)


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
for tag in tags: metadata[tag] = None
if len(info) > 0:
  for line in info:
    for tag in tags:
      if tag in line:
        if not metadata[tag]: metadata[tag] = []
        entry = line.lstrip(tag+' :').rstrip()
        if ':' in entry: entry = entry.split(':')
        metadata[tag].append(entry)
Rmax = 1.
if metadata['Kfactor']:
  Rmax = Rmax/eval(metadata['Kfactor'][0])
  metadata['title'][0] += '*'+metadata['Kfactor'][0]

#Get data:
for pt in pts:
  x,y,res,lim,cond,tot = pt.split()    
#  R = float(eval(tot))/float(eval(lim))
  R = float(eval(res))/float(eval(lim))
  if eval(res) < 0.: continue
  if cond == 'None': cond = '0.'
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

exc.Sort()
x1,y1 = Double(), Double()
exc.GetPoint(0,x1,y1)
yline = []
for ipt in range(exc.GetN()+1): 
  x,y = Double(), Double()
  if ipt < exc.GetN(): exc.GetPoint(ipt,x,y)
  if ipt != exc.GetN() and x == x1: yline.append(y)
  else:
    if len(yline) <= 3 or exc_curve.GetN() == 0: newy = max(yline)
    else:
      newy = max(yline)
      yline = sorted(yline,reverse=True)
      dy = [abs(yline[i]-yline[i+1]) for i in range(len(yline)-1)]
      dmin = min(dy)
      for iD in range(len(dy)-1):
        if dy[iD] == dmin and dy[iD+1] == dmin:
          newy = yline[iD]
          break   
    exc_curve.SetPoint(exc_curve.GetN(),x1,newy)
    x1 = x
    yline = [y]


leg = TLegend(0.6325287,0.7408994,0.9827586,1)
AuxPlot.Default(leg,"Legend")
leg.AddEntry(exc,"Excluded","P")
leg.AddEntry(allow,"Allowed","P")
leg.AddEntry(not_tested,"Not Tested","P")

   
if metadata['Root file'] and os.path.isfile(metadata['Root file'][0]):
  rootfile = TFile(metadata['Root file'][0],"read")
  objs =  gDirectory.GetListOfKeys()
  for ob in objs:
    add = False
    Tob = ob.ReadObj()
    if type(Tob) != type(TGraph()): continue
    if metadata['Root tag']:
      for rootTag in metadata['Root tag']:
        Tag = rootTag
        if type(Tag) == type([]) and len(Tag) > 1: Tag = Tag[0]
        if Tag == ob.GetName():  add = rootTag
    else:
      add = 'Official Exclusion'
    if add:
      exp = Tob
#      exp.Sort()
      exp.SetLineStyle(len(base.GetListOfGraphs())-2)
      base.Add(exp,"L")
      if type(add) == type([]): leg.AddEntry(exp,add[1],"L")
      else: leg.AddEntry(exp,add,"L")


if exc_curve.GetN() > 0:  base.Add(exc_curve,"L")
        

plane = TCanvas("c1", "c1",0,0,800,500)
AuxPlot.Default(plane,"TCanvas")
base.Draw("AP")
AuxPlot.Default(base,"TGraph")
base.GetXaxis().SetTitle("M (GeV)")
base.GetYaxis().SetTitle("m_{LSP} (GeV)")
base.GetYaxis().SetTitleOffset(0.75)
gPad.RedrawAxis()  
leg.Draw()


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

  
