import set_path
import ROOT
from Experiment import SMSResults, SMSInterpolation, LimitGetter, ROOTTools
from Theory import SMSAnalysis
from Tools.PhysicsUnits import addunit, rmvunit

def recreateHist(ana,topo,mz=None,axes=None, run='',line=False,tev=8,nevents=10000,binsize=None):
  """ recreate ROOT TH2F histogram of a given analysis and topology
        needs mz for a topology with intermediate mass,
        axes, for histograms with axes different from M1-M0
        if line=True is selected, produce rootfile with produced histograms
        and line, return filename"""

  toponame=topo
  run1=SMSResults.getRun(ana)
  if mz: toponame=SMSInterpolation.gethistname(topo,mz)
  print toponame
  xmin=rmvunit(SMSResults.getLowX(ana,toponame),'GeV')
  ymin=rmvunit(SMSResults.getLowY(ana,toponame),'GeV')
  xmax=rmvunit(SMSResults.getUpX(ana,toponame),'GeV')
  ymax=rmvunit(SMSResults.getUpY(ana,toponame),'GeV')
  bwx=rmvunit(SMSResults.getBinWidthX(ana,toponame),'GeV')
  bwy=rmvunit(SMSResults.getBinWidthY(ana,toponame),'GeV')
  if binsize:
    bwx=float(binsize)
    bwy=float(binsize)
  print bwx, bwy, xmin, xmax
  bx=int((xmax-xmin)/bwx)
  by=int((ymax-ymin)/bwy)

  h = ROOT.TH2F('h','h',bx,xmin,xmax,by,ymin,ymax)

  if line:
    hL = ROOT.TH2F('hL','hL',bx,xmin,xmax,by,ymin,ymax)
    prod_mode = topo
    if topo=="T1tttt" or topo=="T1bbbb": prod_mode = "T1"
    rootname = "../data/%s_%devts.root" %(prod_mode,nevents)
    if binsize: rootname =  "../data/%s_%devts_%dGeVbin.root" %(prod_mode,nevents,binsize)
    f=ROOT.TFile(rootname)
    if tev==8: rXsName="hist8"
    if tev==7:rXsName="hist7"
    rXs=f.Get(rXsName)

  if line and not rXs:
    print "No refXSec histogram found"
    return None

  x=xmin+bwx/2
  y=ymin+bwy/2

  a=SMSAnalysis.EAnalysis()
  D=None
  L=None

  if mz and mz.find('D')>-1: D=float(mz.split('=')[1])
  if mz and mz.find('LSP')>-1: L=float(mz[mz.find('P')+1:mz.find('P')+4])

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

      v=rmvunit(LimitGetter.GetPlotLimit([massv,massv],[topo,[ana]],a)[0][1],'fb')

      if v: h.Fill(x,y,v)

      if line:
        if v and v<rXs.GetBinContent(rXs.FindBin(x)):
          hL.Fill(x,y,0)
        else: hL.Fill(x,y)

      y+=bwy
    y=ymin+bwy/2
    x+=bwx

  if line:
    f.Close()
    rootname_out="%s_%devts.root" % (topo,nevents)
    if binsize:rootname_out="%s_%devts_%dGeVbin.root" % (topo,nevents,binsize)
    f1=ROOT.TFile(rootname_out,"recreate")
    h.Write()
    hL.Write()
    exclusion=ROOTTools.getTGraphfromContour(hL)
    exclusion.Write()
    f1.Close()
    print rootname_out
    return rootname_out
  return h
