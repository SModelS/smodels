#!/usr/bin/python

#produce png-plot of observed and reproduced exclusion lines
#give ana with -a, topo with -m
#for topologies with intermediate masses: give -mz, if needed -axes

import set_path
import argparse, testGetLimit, ROOT, sys
from Experiment import SMSResults, SMSInterpolation
from Tools.PhysicsUnits import rmvunit

argparser=argparse.ArgumentParser()
argparser.add_argument('-a','--ana',help='input analysis',default='alphaT8TeV')
argparser.add_argument('-t','--topo',help='input topology',default='T2bb')
argparser.add_argument('-mz','--mz',help='intermediate mass information')
argparser.add_argument('-axes','--axes',help='axes information')
argparser.add_argument('-n','--nevts',help='number of events per point in refXSec', type=int, default=10000)
args=argparser.parse_args()

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(0000)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetCanvasBorderMode(0)
#ROOT.gStyle.SetPadLeftMargin(0.16) #set these if needed
ROOT.gStyle.SetPadRightMargin(0.15)
#ROOT.gStyle.SetPadBottomMargin(0.11)
#ROOT.gStyle.SetPadTopMargin(0.07)

tevIn=rmvunit ( SMSResults.getSqrts(args.ana), "TeV" )
run1=SMSResults.getRun(args.ana)

rootname = testGetLimit.recreateHist(args.ana, args.topo, mz=args.mz, axes=args.axes, line=True, tev=tevIn, nevents=args.nevts)

if not rootname:
  print "Could not run recreateHist, png was not produced."
  sys.exit()

toponame = args.topo
if args.mz: toponame = SMSInterpolation.gethistname(args.topo, args.mz)

rootfile = ROOT.TFile(rootname)
ul = rootfile.Get('h')
orig_ul = SMSResults.getUpperLimit(args.ana, toponame)
ul.GetXaxis().SetTitle(SMSResults.particleName(args.topo)+' mass [GeV]')
ul.GetYaxis().SetTitle(orig_ul.GetYaxis().GetTitle())
ul.GetZaxis().SetTitle(orig_ul.GetZaxis().GetTitle().replace('pb','fb'))
exclusion = SMSResults.getExclusionLine(toponame, args.ana, False, plusminussigma=0)
exclusionm1= SMSResults.getExclusionLine(toponame, args.ana, False, plusminussigma=-1)
exclusionp1= SMSResults.getExclusionLine(toponame, args.ana, False, plusminussigma=1)
reproduced_exclusion = rootfile.Get('Graph')
reproduced_exclusion.SetLineWidth(3)
reproduced_exclusion.SetLineColor(ROOT.kRed)

c1=ROOT.TCanvas()
c1.SetLogz()
ul.Draw("COLZ")
exclusion.Draw("SAME")
exclusionm1.Draw("SAME")
exclusionp1.Draw("SAME")
reproduced_exclusion.Draw("SAME")

legend = ROOT.TLegend(0.7, 0.85, 0.15, 0.67)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.AddEntry(exclusion, 'exclusion reported by experiment','L')
legend.AddEntry(exclusionm1, 'exclusion #pm 1#sigma reported by experiment','L')
legend.AddEntry(reproduced_exclusion, "SModelS exclusion",'L')
legend.Draw()

title = ROOT.TLatex()
title.SetNDC()
title.DrawLatex(0.1,0.93, SMSResults.getPrettyName(args.ana)+', '+str(int(tevIn))+' TeV, NLONLL,')
#title = ROOT.TMathText()
#title.SetNDC()
#title.DrawMathText(0.1,0.93, SMSResults.getPrettyName(args.ana)+', '+str(int(tevIn))+' TeV, NLONLL,')

decay = ROOT.TLatex()
decay.SetNDC()
decay.DrawLatex(0.55, 0.93, SMSResults.SMSInfo("decay", args.topo, args.ana))

c1.Print("../plots/%s_%s_%devts.png" %(args.ana,toponame,args.nevts))

rootfile.Close()
