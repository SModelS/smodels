#!/usr/bin/env python

import argparse
import os, sys

argparser = argparse.ArgumentParser(description='plots the exclusion region and envelope')
argparser.add_argument( 'input', help='input data file (mother_mass,LSP_mass,R)')
args=argparser.parse_args()

from ROOT import *
import AuxPlot

theo = TH2F("","",100,0.,2000.,100,0.,1000.)
exp_limit = TH2F("","",100,0.,2000.,100,0.,1000.)
exp_data = TH2F("","",100,0.,2000.,100,0.,1000.)

infile = open(args.input)
data = infile.read()
pts = data[:data.find('#END')-1].split('\n')
info = data[data.find('#END'):].split('\n')

#Get data:
for pt in pts:
  x,y,res,lim,tot = pt.split()    
  x = eval(x)
  y = eval(y)
  res = eval(res)
  lim = eval(lim)
  if res <= 0.: continue
  theo.Fill(x,y,res)
  if exp_limit.GetBinContent(exp_limit.FindBin(x,y)) < lim:
    exp_limit.SetBinContent(exp_limit.FindBin(x,y),lim)
  if theo.GetBinContent(theo.FindBin(x,y)) < res:
    theo.SetBinContent(theo.FindBin(x,y),res)
  
infile.close()
infile = open('T2data.dat','r')
pts = infile.read()
pts = pts.split('\n')
pts = pts[:-1]
for pt in pts:
  x,y,lim = pt.split()
  x = eval(x)
  y = eval(y)
  lim = eval(lim)
  if exp_data.GetBinContent(exp_data.FindBin(x,y)) < lim:
    exp_data.SetBinContent(exp_data.FindBin(x,y),lim)

AuxPlot.set_palette(gStyle)
c1 = TCanvas()
c1.SetLogz()
exp_limit.Draw('COLZ')
exp_limit.GetZaxis().SetRangeUser(0.5,10.**3)
#c2 = TCanvas()
#c2.SetLogz()
#exp_data.Draw('COLZ')
c3 = TCanvas()
c3.SetLogz()
theo.Draw('COLZ')
theo.GetZaxis().SetRangeUser(0.5,10.**3)



ans = raw_input("Hit any key to close\n")

  
