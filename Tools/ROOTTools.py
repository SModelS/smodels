#!/usr/bin/env python

"""
.. module:: ROOTTools
    :synopsis: Collection of methods used in the context of ROOT

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

def getRootVersionFromConfig_ ( astuple=False ):
  """ get the ROOT version, via root-config --version

    :param astuple: false returns string, true returns tuple of integers.
    :returns: ROOT version
  """
  from VariousHelpers import logging
  log = logging.getLogger(__name__)
  try:
    import commands
    S=commands.getoutput("root-config --prefix")
    if S.find("not found")>-1:
      log.error ( S )
      return None
    if not astuple: return S
    T,C=S.split("/")
    A,B=T.split(".")
    return (int(A),int(B),int(C))
  except Exception,e:
    log.error ( e )
    return None


def getRootVersion ( astuple=False, useconfig=False ):
  """ get the ROOT version.

    :param astuple: false returns string, true returns tuple of integers.
    :returns: ROOT version
  """
  if useconfig: return getRootVersionFromConfig_ ( astuple )
  from VariousHelpers import logging
  log = logging.getLogger(__name__)
  try:
    import ROOT
    S=ROOT.gROOT.GetVersion()
    if not astuple: return S
    T,C=S.split("/")
    A,B=T.split(".")
    return (int(A),int(B),int(C))
  except Exception,e:
    log.error ( e )
    return None
    
def getRootPath ( ):
  """ get the ROOT path, first try via root-config, then query ROOTSYS
    :returns: ROOT path
  """
  import logging
  log = logging.getLogger(__name__)
  try:
    import commands
    out=commands.getoutput("root-config --prefix")
    if out.find("not found")>-1:
      log.info ( out )
      import os
      ret=os.getenv("ROOTSYS")
      if not ret:
        log.error ( "ROOTSYS not set, either" )
      return None
    else:
      return out
  except Exception,e:
    log.error ( e )
    return None
      
def getRootLibraryPath ( ):
  """ get the ROOT library path, first try via root-config, then systematically
      try candidate paths.

    :returns: ROOT library path
  """
  import logging
  log = logging.getLogger(__name__)
  try:
    import commands
    out=commands.getoutput("root-config --libdir")
    if out.find("not found")>-1:
      log.info ( out )
      import os
      ret=os.getenv("ROOTSYS")
      for Dir in [ "lib/x86_64-linux-gnu", "lib64/root", "lib/root" ]:
        F="%s/%s/libRint.so" % ( ret, Dir )
        if os.path.exists ( F ):
          return F
      log.error ( "no suitable libdir found." )
      return None
    else:
      return out
  except Exception,e:
    log.error ( e )
    return None

def getRootPythonPath ( ):
  """ get the ROOT python path, via .getRootLibraryPath() and .getRootVersion().

    :returns: ROOT python path
  """
  import os
  version = getRootVersion(True)
  libpath = getRootLibraryPath()
  if not version or not libpath:
    return None
  V=str(version[0])+"."+str(version[1])
  for SubDir in [ V, "root"+V ]:
    Dir=libpath+"/"+SubDir
    if os.path.exists ( Dir+"/ROOT.py" ):
      return Dir 
  return None
    
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
