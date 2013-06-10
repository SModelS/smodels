#!/usr/bin/python

import set_path, ROOT
ROOT.gROOT.SetBatch()

from numpy import arange
from Tools import PhysicsUnits
from Experiment import SMSResults as S
S.useUnits(False)
a="Weakinos8TeV"
tslep="TChiChipmSlepStau"
tstau="TChiChipmStauStau"

xRange=S.getLowX( a, tslep ), S.getUpX( a, tslep )
xWidth=S.getBinWidthX ( a, tslep )
yRange=S.getLowY( a, tslep ), S.getUpY( a, tslep )
yWidth=S.getBinWidthY ( a, tslep )
print "x,y=",xRange,yRange
print "wx,wy=",xWidth,yWidth
xRange=S.getLowX( a, tstau ), S.getUpX( a, tstau )
xWidth=S.getBinWidthX ( a, tstau )
yRange=S.getLowY( a, tstau ), S.getUpY( a, tstau )
yWidth=S.getBinWidthY ( a, tstau )
print "x,y=",xRange,yRange
print "wx,wy=",xWidth,yWidth

hslepstau=S.getUpperLimit ( a, tslep )
print hslepstau
hstaustau=S.getUpperLimit ( a, "TChiChipmStauStau" )
print hstaustau

hk=hslepstau.Clone()
for x in arange ( xRange[0]-xWidth, xRange[1]+xWidth, xWidth ):
  for y in arange ( yRange[0]-yWidth, yRange[1]+yWidth, yWidth ):
    Bin=hk.FindBin(x,y)
    Binhstaustau=hk.FindBin(x,y)
    numerator=hk.GetBinContent(Bin )
    denom=hstaustau.GetBinContent(Binhstaustau)
    if denom>0.:
      hk.SetBinContent ( Bin, numerator / denom )

hk.GetZaxis().SetRangeUser(0.,2.)
hk.Draw("colz")

ROOT.c1.Print("k.png")

