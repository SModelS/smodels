#!/usr/bin/env python

import ROOT

f = ROOT.TFile("CMS-SUS-16-049_Figure_007.root", "update")

h = ROOT.gDirectory.Get("hXsec_obs_corr")

hFix = h.Clone()
hFix.SetName("hXsec_obs_corrRemoved")
for i in range(1,hFix.GetNbinsX() + 1):
    for j in range(1,hFix.GetNbinsY() + 1):
        xval = hFix.GetXaxis().GetBinCenter(i)
        yval = hFix.GetYaxis().GetBinCenter(j)
        if -200. < yval-xval < -150. and  xval+yval < 370.:
            ibin = hFix.GetBin(i,j)
            hFix.SetBinContent(ibin,0.)
hFix.Draw("COLZ")

hFix.Write("hXsec_obs_corrRemoved")
f.Close()

