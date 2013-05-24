""" collection of tools needed for use and manipulation of ROOT files """


def getTGraphfromContour(exclhisto):
  import ROOT
  """ returns the contour of an exclusion histogram as TGraph"""
  ROOT.gROOT.SetBatch()
  c1 = ROOT.TCanvas()
  c1.cd()
  exclhisto.Draw("CONT LIST")
  ROOT.gPad.Update()
  ROOT.gROOT.GetListOfSpecials().ls()
  gr = ROOT.gROOT.GetListOfSpecials().FindObject('contours')('TList').At(0)
  return gr

