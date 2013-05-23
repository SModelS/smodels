""" collection of tools needed for use and manipulation of ROOT files """

import ROOT


def getTGraphfromContour(exclhisto):
   """ returns the contour of an exculsion histogram as TGraph"""
   ROOT.gROOT.SetBatch()
   c1 = ROOT.TCanvas()
   c1.cd()
   exclhisto.Draw("CONT LIST")
   ROOT.gPad.Update()
   ROOT.gROOT.GetListOfSpecials().ls()
   gr = ROOT.gROOT.GetListOfSpecials().FindObject('contours')('TList').At(0)
   return gr
 
