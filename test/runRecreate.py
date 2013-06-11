#!/usr/bin/python

#produce png-plot of observed and reproduced exclusion lines
#give ana with -a, topo with -t
#for topologies with intermediate masses: give -mz, if needed -axes

import set_path
import argparse, testGetLimit, ROOT, sys, TestTools, os
from Experiment import SMSResults, SMSInterpolation
from Tools.PhysicsUnits import rmvunit
from Tools import ROOTTools

def zero(x,y):
  if (x+2)%4==0 and (y+2)%4==0: return False                                   
  if x%4!=0 or y%4!=0: return True                                              
  return False

argparser=argparse.ArgumentParser()
argparser.add_argument('-a','--ana',help='input analysis [alphaT8TeV]',default='alphaT8TeV')
argparser.add_argument('-t','--topo',help='input topology [T2bb]',default='T2bb')
argparser.add_argument('-mz','--mz',help='intermediate mass information')
argparser.add_argument('-axes','--axes',help='axes information')
argparser.add_argument('-n','--nevts',help='number of events per point in refXSec [10000]', type=int, default=10000)
argparser.add_argument('-b','--binsize',help='binsize in GeV', type=int)
argparser.add_argument('-text','--text',help='print upper limit each 100 GeV',action='store_true')

args=argparser.parse_args()

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(0000)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetCanvasBorderMode(0)
##ROOT.gStyle.SetPadLeftMargin(0.16) #set these if needed
ROOT.gStyle.SetPadRightMargin(0.18)
#ROOT.gStyle.SetPadBottomMargin(0.11)
#ROOT.gStyle.SetPadTopMargin(0.07)
ROOTTools.useNiceColorPalette ( "temperature", 0.5, 50 ) 

tevIn=rmvunit ( SMSResults.getSqrts(args.ana), "TeV" )
run1=SMSResults.getRun(args.ana)
binsize=args.binsize
rootname = testGetLimit.recreateHist(args.ana, args.topo, mz=args.mz, axes=args.axes, line=True, tev=tevIn, nevents=args.nevts, binsize=args.binsize)

if not rootname:
  print "Could not run recreateHist, png was not produced."
  sys.exit()

toponame = args.topo
if args.mz: toponame = SMSInterpolation.gethistname(args.topo, args.mz)

rootfile = ROOT.TFile(rootname)
ul = rootfile.Get('h')
orig_ul = SMSResults.getUpperLimit(args.ana, toponame)
ul.GetXaxis().SetTitle(SMSResults.particleName(args.topo)+' mass [GeV]')
if not SMSResults.hasDictionary(args.ana):
  ul.GetYaxis().SetTitle(orig_ul.GetYaxis().GetTitle())
  ul.GetZaxis().SetTitle('#splitline{%s}{reported by experiment}' % orig_ul.GetZaxis().GetTitle().replace('pb','fb'))
else:
  ul.GetYaxis().SetTitle('LSP mass [GeV]')
  if args.axes and SMSInterpoltaion.getaxis(y,args.axes)!="M0":
    ul.GetYaxis().SetTitle(SMSInterpolation.getaxis(y,args.axes)+" [GeV]")
  ul.GetZaxis().SetTitle('#splitline{95% CL upper limit on #sigma (fb)}{reported by experiment}')
ul.GetZaxis().SetTitleOffset(1.4)
ul.GetYaxis().SetTitleOffset(1.15)
exclusion = SMSResults.getExclusionLine(toponame, args.ana, False, plusminussigma=0)
exclusionm1= SMSResults.getExclusionLine(toponame, args.ana, False, plusminussigma=-1)
exclusionp1= SMSResults.getExclusionLine(toponame, args.ana, False, plusminussigma=1)
reproduced_exclusion = rootfile.Get('Graph')
reproduced_exclusion.SetLineWidth(3)
reproduced_exclusion.SetLineColor(ROOT.kRed)

c1=ROOT.TCanvas()
c1.SetLogz()
ul.Draw("COLZ")
if exclusion: exclusion.Draw("SAME")
if exclusionm1: exclusionm1.Draw("SAME")
if exclusionp1: exclusionp1.Draw("SAME")
reproduced_exclusion.Draw("SAME")

if args.text:
  ulText = ul.Clone("ulText")
  for x in range(1,ulText.GetXaxis().GetLast()+1):
    for y in range(1,ulText.GetYaxis().GetLast()+1):
      bincontent=ulText.GetBinContent(x,y)
      if zero(x,y): bincontent=0
      ulText.SetBinContent(x,y,bincontent)
  ulText.Draw("SAMETEXT") 

legend = ROOT.TLegend(0.7, 0.85, 0.15, 0.67)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
if exclusion: legend.AddEntry(exclusion, 'exclusion reported by experiment (#pm 1#sigma)','L')
#legend.AddEntry(exclusionm1, 'exclusion #pm 1#sigma reported by experiment','L')
#legend.AddEntry(reproduced_exclusion, "SModelS exclusion",'L')
legend.AddEntry(reproduced_exclusion, "                exclusion",'L')
legend.SetMargin(0.12)
legend.Draw()

if exclusionm1:
  line = ROOT.TLine()
  line.SetLineStyle(exclusionm1.GetLineStyle())
  line.SetLineWidth(exclusionm1.GetLineWidth())
  line.DrawLineNDC(0.21, 0.785, 0.16, 0.785)
  line.DrawLineNDC(0.21, 0.825, 0.16, 0.825)

#title = ROOT.TMathText()
title = ROOT.TLatex()
title.SetNDC()

#title_str='%s , %s (%s-%s) , %d TeV , NLONLL' %(SMSResults.SMSInfo("decay", args.topo, args.ana),SMSResults.getPrettyName(args.ana),SMSResults.getExperiment(args.ana),SMSResults.getPAS(args.ana),tevIn)
#title.DrawMathText(0.1,0.93,title_str.replace('#','\\'))

order=SMSResults.getOrder ( args.ana )
if order=="nlonll":
  order="NLO-NLL"
if run1!="ATLAS8TeV":
  title.DrawLatex(0.1,0.93, '%s, %s (%s-%s), %d TeV, %s' %(SMSResults.SMSInfo("decay", args.topo, args.ana),SMSResults.getPrettyName(args.ana),SMSResults.getExperiment(args.ana),SMSResults.getPAS(args.ana),tevIn, order))
else:
  title.DrawLatex(0.1,0.93, '%s (%s), %d TeV, %s' %(SMSResults.SMSInfo("decay", args.topo, args.ana),SMSResults.getPAS(args.ana),tevIn, order ))

img=ROOT.TImage.Create()
img.FromPad(c1)
logo=ROOT.TASImage("../plots/smodels75.png")
logo.Draw("SAME")
img.Merge(logo,"alphablend",147,123)
if not args.binsize:
  plotname= "../plots/%s_%s_%devts.png" %(args.ana,toponame,args.nevts)
else:
  plotname= "../plots/%s_%s_%devts_%dGeVbin.png" %(args.ana,toponame,args.nevts,args.binsize)
print "[runRecreate.py] now write",plotname
img.WriteImage(plotname)


#c1.Print("../plots/%s_%s_%devts.png" %(args.ana,toponame,args.nevts))
#c1.Print("../plots/%s_%s_%devtsROOT.pdf" %(args.ana,toponame,args.nevts))
#TestTools.convertROOTpdf("../plots/%s_%s_%devtsROOT.pdf" %(args.ana,toponame,args.nevts), "../plots/%s_%s_%devts.pdf" %(args.ana,toponame,args.nevts),"../plots/%s_%s_%devts.png" %(args.ana,toponame,args.nevts))#convertROOTpdf not working like it should be!
#os.remove("../plots/%s_%s_%devtsROOT.pdf" %(args.ana,toponame,args.nevts))

rootfile.Close()
