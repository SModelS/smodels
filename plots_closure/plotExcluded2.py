#!/usr/bin/env python

import argparse
import os, sys

argparser = argparse.ArgumentParser(description='plots the exclusion region and envelope')
argparser.add_argument( 'input', help='input data file (mother_mass,LSP_mass,R)')
args=argparser.parse_args()

from ROOT import *
import AuxPlot

exc = TGraph()
exc.SetMarkerStyle(20)
exc.SetMarkerColor(kRed-2)
exc_curve = TGraph()
exc_curve.SetLineColor(kRed)
exc_curve.SetLineStyle(9)
exc_curve.SetLineWidth(4)


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
xv = []
yv = []
zv = []
for pt in pts:
  x,y,res,lim,cond,tot = pt.split() 
  if eval(res) < 0.: continue
  x = eval(x)
  y = eval(y)
  res = eval(res)
  lim = eval(lim)
  xv.append(x)
  yv.append(y)
  zv.append(lim)
  R = res/lim
  if R >= Rmax: exc.SetPoint(exc.GetN(),x,y)
   
xv = sorted(set(xv))    
yv = sorted(set(yv))
zv = sorted(set(zv))
xmax = max(xv)
ymax = max(yv)  
zmax = max(zv)
xmin = min(xv)
ymin = min(yv)
zmin = min(zv)
dxmin = min([abs(xv[i]-xv[i+1]) for i in range(len(xv)-1)])*1.1
dymin = min([abs(yv[i]-yv[i+1]) for i in range(len(yv)-1)])*1.1
exp_limit = TH2F("","",int((xmax-xmin)/dxmin),xmin,xmax,int((ymax-ymin)/dymin),ymin,ymax)
for pt in pts:
  x,y,res,lim,cond,tot = pt.split()    
  x = eval(x)
  y = eval(y)
  lim = eval(lim)
  if eval(res) < 0.: continue
  if exp_limit.GetBinContent(exp_limit.FindBin(x,y)) < lim:
    exp_limit.SetBinContent(exp_limit.FindBin(x,y),lim)


AuxPlot.set_palette(gStyle)
plane = TCanvas("c1", "c1",0,0,800,500)
AuxPlot.Default(plane,"TCanvas")
plane.SetRightMargin(0.14700422)
plane.SetTopMargin(0.05296053)
plane.SetBottomMargin(0.16796053)
exp_limit.SetStats(kFALSE)
AuxPlot.Default(exp_limit,'TH2')
c1.SetLogz()
exp_limit.Draw('COLZ')
exp_limit.GetXaxis().SetTitle("M (GeV)")
exp_limit.GetYaxis().SetTitle("m_{LSP} (GeV)")
exp_limit.GetYaxis().SetTitleOffset(0.75)
exp_limit.GetZaxis().SetTitle("Upper Limit (fb)")
exp_limit.GetZaxis().SetTitleOffset(0.75)
exp_limit.GetZaxis().SetRangeUser(zmin,zmax)
gPad.RedrawAxis()  


exc.Sort()
x1,y1 = Double(), Double()
exc.GetPoint(0,x1,y1)
yline = []
for ipt in range(exc.GetN()+1): 
  x,y = Double(), Double()
  dmin = 0.
  if ipt < exc.GetN(): exc.GetPoint(ipt,x,y)
  if ipt != exc.GetN() and x == x1: yline.append(y)
  else:
    yline = sorted(yline,reverse=True)
    dy = [abs(yline[i]-yline[i+1]) for i in range(len(yline)-1)]
    if len(yline) <= 3 or exc_curve.GetN() == 0:
      newy = max(yline)
      if len(dy) > 2: dmin = min([abs(yline[i]-yline[i+1]) for i in range(len(yline)-1)])
    else:
      newy = max(yline)     
      dmin = min(dy)
      for iD in range(len(dy)-1):
        if dy[iD] == dmin and dy[iD+1] == dmin:
          newy = yline[iD]
          break
    exc_curve.SetPoint(exc_curve.GetN(),x1,newy+dmin/2.)
    x1 = x
    yline = [y]

x2,y2 = Double(), Double()
exc_curve.GetPoint(exc_curve.GetN()-1,x2,y2)
exc_curve.SetPoint(exc_curve.GetN(),x2,0.)

leg = TLegend(0.6934673,0.8051392,0.9886935,0.99)
AuxPlot.Default(leg,"Legend")

ngraphs = 0
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
      exp.SetLineStyle(ngraphs+1)
      exp.SetLineWidth(4)
      exp.Draw("SAMEL")
      ngraphs += 1
      if type(add) == type([]): leg.AddEntry(exp,add[1],"L")
      else: leg.AddEntry(exp,add,"L")
      
if exc_curve.GetN() > 0:
  exc_curve.Draw("SAMEL")
  leg.AddEntry(exc_curve,'SmodelS',"L")

leg.Draw()

if metadata['title']:
  title = metadata['title'][0]
  tit = TPaveLabel(0.01005025,0.8672377,0.5150754,0.9807281,title,"brNDC")
  tit.SetBorderSize(4)
  tit.SetFillColor(0)
  tit.SetTextFont(42)
  tit.SetTextSize(0.2727273)
  tit.Draw()

gPad.Update()
if metadata['Out file']: c1.Print('2D-'+metadata['Out file'][0])
ans = raw_input("Hit any key to close\n")

  
