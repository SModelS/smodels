""" collection of tools needed for use and manipulation of ROOT files """


def getTGraphFromContour(exclhisto):
  """ returns the contour of an exclusion histogram as TGraph"""
  import ROOT
  ROOT.gROOT.SetBatch()
  c1 = ROOT.TCanvas()
  c1.cd()
  exclhisto.Draw("CONT LIST")
  ROOT.gPad.Update()
  ROOT.gROOT.GetListOfSpecials().ls()
  gr = ROOT.gROOT.GetListOfSpecials().FindObject('contours')('TList').At(0)
  return gr

def useNiceColorPalette( palette="temperature", f=0., ngradientcolors=20 ):
  """ create a fine-grained temperature color palette,
      :param palette: which palette. values are temperature, blackwhite, darkbody, deepsea, blueyellow, rainbow, inverteddarkbody, yellowpurple, greenpurple, bluepurple
      :type palette: str
      :param ngradientcolors: how many colors
      :type ngradientcolors: int
      :param f: 0 means full color, 0.5 lighter palette, 1.0 all white
      :type f: float
  """
  from array import array
  import ROOT
  foundpalette=False
  stops,red,green,blue=[],[],[],[]
  if palette=="brown":
    foundpalette=True
    stops = [0.0, 0.5, 0.75, 1.0 ]
    red   = [106./256, 0.38, 1.00, 1.0 ]
    green = [44. /256, 0.21, 0.81, 0.9 ]
    blue  = [0.05 , 0.07, 0.05, 0.05 ]
  if palette=="temperature":
    foundpalette=True
    stops = [0.00, 0.34, 0.61, 0.84, 1.00]
    red   = [0.00, 0.00, 0.87, 1.00, 0.51]
    green = [0.00, 0.91, 1.00, 0.20, 0.00]
    blue  = [0.51, 1.00, 0.12, 0.00, 0.00]
  if palette=="yellowpurple":
    foundpalette=True
    red = [ 0., 0.0, 1.0, 1.0, 1.0 ]
    green = [ 0., 0.0, 0.0, 1.0, 1.0 ]
    blue  = [ 0., 1.0, 0.0, 0.0, 1.0 ]
    stops = [ 0., .25, .50, .75, 1.0 ]
  if palette=="greenpurple":
    foundpalette=True
    red = [ 1.00, 0.50, 0.00 ]
    green = [ 0.50, 0.00, 1.00 ]
    blue  = [ 1.00, 0.00, 0.50 ]
    stops = [ 0.00, 0.50, 1.00 ]
  if palette=="bluepurple":
    foundpalette=True
    red   = [ 1.00, 0.00, 0.00 ]
    green = [ 0.00, 1.00, 0.00 ]
    blue  = [ 1.00, 0.00, 1.00 ]
    stops = [ 0.00, 0.50, 1.00 ]

  if foundpalette:
    if f!=None and f>0.0:
      for (i,r) in enumerate(red):
        r=1-(1.0-r)*f
        red[i]=r
      for (i,r) in enumerate(green):
        r=1-(1.0-r)*f
        green[i]=r
      for (i,r) in enumerate(blue):
        r=1-(1.0-r)*f
        blue[i]=r
    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)
    nstops = len(s)
    ROOT.TColor.CreateGradientColorTable(nstops, s, r, g, b, ngradientcolors )
    ROOT.gStyle.SetNumberContours( ngradientcolors )
    ROOT.gStyle.cd()
  if palette=="deepsea":
    foundpalette=True
    ROOT.gStyle.SetPalette(51) ## black-and-white
  if palette=="blackwhite":
    foundpalette=True
    ROOT.gStyle.SetPalette(52) ## black-and-white
  if palette=="darkbody":
    foundpalette=True
    ROOT.gStyle.SetPalette(53) 
  if palette=="blueyellow":
    foundpalette=True
    ROOT.gStyle.SetPalette(54)
  if palette=="rainbow":
    foundpalette=True
    ROOT.gStyle.SetPalette(55)
  if palette=="inverteddarkbody":
    foundpalette=True
    ROOT.gStyle.SetPalette(56) 
  if not foundpalette:
    print "[ROOTTools.py] error: did not find palette %s. Existing palettes are: temperature, blackwhite, darkbody, deepsea, blueyellow, rainbow, inverteddarkbody "
