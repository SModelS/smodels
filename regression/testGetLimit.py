#!/usr/bin/python

import sys
sys.path.append ( "../" )

import ROOT, SMSmethods, SMSgetlimit
from Experiment import SMSResults, SMSInterpolation
# from Experiment.SMSUnits import addunit, rmvunit

def recreateHist(ana,topo,mz=None,axes=None, run=''):
  """ recreate ROOT TH2F histogram of a given analysis and topology (only topologies with no intermediate masses) """
  if SMSResults.hasDictionary(ana):
    print "Cannot recreate histogram from dictionary"
    return None
  print mz, axes
  toponame=topo
  if mz: toponame=SMSInterpolation.gethistname(topo,mz)
  xmin=rmvunit(SMSResults.getLowX(ana,toponame),'GeV')
  ymin=rmvunit(SMSResults.getLowY(ana,toponame),'GeV')
  xmax=rmvunit(SMSResults.getUpX(ana,toponame),'GeV')
  ymax=rmvunit(SMSResults.getUpY(ana,toponame),'GeV')
  bwx=rmvunit(SMSResults.getBinWidthX(ana,toponame),'GeV')
  bwy=rmvunit(SMSResults.getBinWidthY(ana,toponame),'GeV')
  bx=int((xmax-xmin)/bwx)
  by=int((ymax-ymin)/bwy)
#  print bx,by
#  xmin=0.
#  xmax=1000.
#  bwx=20.
#  ymin=0.
#  ymax=1000.
#  bwy=20.
#  bx=int((xmax-xmin)/bwx)
 #       by=int((ymax-ymin)/bwy)
  h = ROOT.TH2F('h','h',bx,xmin,xmax,by,ymin,ymax)
  x=xmin+bwx/2
  y=ymin+bwy/2
  a=SMSmethods.EAnalysis()
  D=None
  L=None
  if mz.find('D')>-1: D=float(mz.split('=')[1])
  if mz.find('LSP')>-1: L=float(mz[mz.find('P')+1:mz.find('P')+4])
  while x<xmax:
    while y<ymax:
      if D:
        massv=[0.,0.,0.]
        massv[SMSInterpolation.getaxis('x',axes)]=x
        massv[SMSInterpolation.getaxis('y',axes)]=y
        if massv[SMSInterpolation.getaxis('x',mz)]==0.: massv[SMSInterpolation.getaxis('x',mz)]=massv[SMSInterpolation.getaxis('y',mz)]+D
        if massv[SMSInterpolation.getaxis('y',mz)]==0.: massv[SMSInterpolation.getaxis('y',mz)]=massv[SMSInterpolation.getaxis('x',mz)]-D
      elif L:
        massv=[0.,0.,L]
        massv[SMSInterpolation.getaxis('x',axes)]=x
        massv[SMSInterpolation.getaxis('y',axes)]=y
      elif mz: massv=[x,SMSInterpolation.getxval(x,y,mz,mass=True),y]
      else: massv=[x,y]
      v=rmvunit(SMSgetlimit.GetPlotLimit([massv,massv],[topo,[ana]],a)[0][1],'pb')
      if v: h.Fill(x,y,v)
      y+=bwy
    y=ymin+bwy/2
    x+=bwx
  return h
