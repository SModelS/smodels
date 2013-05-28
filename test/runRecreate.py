#!/usr/bin/python

#produce png-plot of observed and reproduced exclusion lines
#give ana with -a, topo with -m
#for topologies with intermediate masses: give -mz, if needed -axes
#default sqrt(s) is 8 TeV, set -tev to change to 7 TeV

import set_path
import argparse, testGetLimit, ROOT
from Experiment import SMSResults, SMSInterpolation

argparser=argparse.ArgumentParser()
argparser.add_argument('-a','--ana',help='input analysis')
argparser.add_argument('-t','--topo',help='input topology')
argparser.add_argument('-mz','--mz',help='intermediate mass information')
argparser.add_argument('-axes','--axes',help='axes information')
argparser.add_argument('-tev','--tev',type=int,help='input sqrt(s)')
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
run1=SMSResults.getRun(args.ana)

rootname = testGetLimit.recreateHist(args.ana, args.topo, mz=args.mz, axes=args.axes, line=True, tev=tevIn)

toponame = args.topo
if args.mz: toponame = SMSInterpolation.gethistname(args.topo, args.mz)

rootfile = ROOT.TFile(rootname)
orig_rootfile = f2=ROOT.TFile("/afs/hephy.at/user/w/walten/public/sms/%s/%s/sms.root" % (run1, args.ana))
ul = rootfile.Get('h')
orig_ul = orig_rootfile.Get('limit_%s' %toponame)
ul.GetXaxis().SetTitle(orig_ul.GetXaxis().GetTitle())
ul.GetYaxis().SetTitle(orig_ul.GetYaxis().GetTitle())
ul.GetZaxis().SetTitle(orig_ul.GetZaxis().GetTitle())
reproduced_exclusion = rootfile.Get('Graph')
reproduced_exclusion.SetLineColor(ROOT.kRed)
exclusion = orig_rootfile.Get('exclusion_%s' % toponame)
exclusionm1= orig_rootfile.Get('exclusionm1_%s' % toponame)
exclusionp1= orig_rootfile.Get('exclusionp1_%s' % toponame)

c1=ROOT.TCanvas()
ul.Draw("COLZ")
exclusion.Draw("SAME")
exclusionm1.Draw("SAME")
exclusionp1.Draw("SAME")
reproduced_exclusion.Draw("SAME")

legend = ROOT.TLegend(0.45, 0.8, 0.2, 0.65)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.AddEntry(reproduced_exclusion, "reproduced exclusion",'L')
legend.AddEntry(exclusion, 'observed exclusion','L')
legend.AddEntry(exclusionm1, 'observed exclusion #pm 1#sigma','L')
legend.Draw()

c1.Print("../plots/%s_%s.png" %(args.ana,toponame))

rootfile.Close()
orig_rootfile.Close()
