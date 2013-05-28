#!/usr/bin/python

#produce png-plot of observed and reproduced exclusion lines
#give ana with -a, topo with -m
#for topologies with intermediate masses: give -mz, if needed -axes
#default sqrt(s) is 8 TeV, set -tev to change to 7 TeV

import set_path
import argparse, testGetLimit, ROOT, sys
from Experiment import SMSResults, SMSInterpolation, SMSResultsCollector

argparser=argparse.ArgumentParser()
argparser.add_argument('-a','--ana',help='input analysis',default='alphaT8TeV')
argparser.add_argument('-t','--topo',help='input topology',default='T2bb')
argparser.add_argument('-mz','--mz',help='intermediate mass information')
argparser.add_argument('-axes','--axes',help='axes information')
argparser.add_argument('-tev','--tev',type=int,help='input sqrt(s)',default=8)
args=argparser.parse_args()

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(0000)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetCanvasBorderMode(0)
#ROOT.gStyle.SetPadLeftMargin(0.16) #set these if needed
ROOT.gStyle.SetPadRightMargin(0.15)
#ROOT.gStyle.SetPadBottomMargin(0.11)
#ROOT.gStyle.SetPadTopMargin(0.07)

tevIn=8
if args.tev: tevIn=args.tev
tevIn=SMSResultsCollector.getSqrts(args.ana)
run1=SMSResults.getRun(args.ana)

rootname = testGetLimit.recreateHist(args.ana, args.topo, mz=args.mz, axes=args.axes, line=True, tev=tevIn)

if not rootname:
  print "Could not run recreateHist, png was not produced."
  sys.exit()

toponame = args.topo
if args.mz: toponame = SMSInterpolation.gethistname(args.topo, args.mz)

rootfile = ROOT.TFile(rootname)
ul = rootfile.Get('h')
orig_ul = SMSResults.getUpperLimit(args.ana, toponame)
ul.GetXaxis().SetTitle(SMSResultsCollector.particleName(args.topo)+' mass [GeV]')
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

legend = ROOT.TLegend(0.5, 0.85, 0.15, 0.67)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.AddEntry(reproduced_exclusion, "reproduced exclusion",'L')
legend.AddEntry(exclusion, 'observed exclusion','L')
legend.AddEntry(exclusionm1, 'observed exclusion #pm 1#sigma','L')
legend.Draw()

title = ROOT.TLatex()
title.SetNDC()
title.DrawLatex(0.1,0.93, args.ana+', '+str(tevIn)+' TeV, NLONLL')

decay = ROOT.TLatex()
decay.SetNDC()
decay.DrawLatex(0.55, 0.93, SMSResultsCollector.SMSInfo("decay", args.topo, args.ana))

c1.Print("../plots/%s_%s.png" %(args.ana,toponame))

rootfile.Close()
