import set_path
import ROOT
from Experiment import SMSResults, SMSInterpolation, SMSgetlimit, ROOTTools
from Theory import SMSDataObjects
from Tools.PhysicsUnits import addunit, rmvunit

def recreateHist(ana,topo,mz=None,axes=None, run='',line=False,tev=8):
  """ recreate ROOT TH2F histogram of a given analysis and topology
        needs mz for a topology with intermediate mass,
        axes, for histograms with axes different from M1-M0
        if line=True is selected, returns histogram and exclusion line """

  toponame=topo
  run1=SMSResults.getRun(ana)
  if mz: toponame=SMSInterpolation.gethistname(topo,mz)

  xmin=rmvunit(SMSResults.getLowX(ana,toponame),'GeV')
  ymin=rmvunit(SMSResults.getLowY(ana,toponame),'GeV')
  xmax=rmvunit(SMSResults.getUpX(ana,toponame),'GeV')
  ymax=rmvunit(SMSResults.getUpY(ana,toponame),'GeV')
  bwx=rmvunit(SMSResults.getBinWidthX(ana,toponame),'GeV')
  bwy=rmvunit(SMSResults.getBinWidthY(ana,toponame),'GeV')

  bx=int((xmax-xmin)/bwx)
  by=int((ymax-ymin)/bwy)

  h = ROOT.TH2F('h','h',bx,xmin,xmax,by,ymin,ymax)

  if line:
    hL = ROOT.TH2F('hL','hL',bx,xmin,xmax,by,ymin,ymax)
    f=ROOT.TFile("../data/%s.root" %topo)
    if tev==8: rXsName="hist8"
    if tev==7:rXsName="hist7"
    rXs=f.Get(rXsName)

  if line and not rXs:
    print "No refXSec histogram found"
    return None

  x=xmin+bwx/2
  y=ymin+bwy/2

  a=SMSDataObjects.EAnalysis()
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

      v=rmvunit(SMSgetlimit.GetPlotLimit([massv,massv],[topo,[ana]],a)[0][1],'pb')

      if v: h.Fill(x,y,v)

      if line:
        if v and v<rXs.GetBinContent(rXs.FindBin(x)):
          hL.Fill(x,y)

      y+=bwy
    y=ymin+bwy/2
    x+=bwx

  f.Close()
  if line:
    f1=ROOT.TFile("%s.root" % topo,"recreate")
    h.Write()
    exclusion=ROOTTools.getTGraphfromContour(hL)
    exclusion.Write()
    f1.Close()
    f2=ROOT.TFile("/afs/hephy.at/user/w/walten/public/sms/%s/%s/sms.root" % (run1, ana))
    orig_exclusion=f2.Get("exclusion_%s" % toponame)
    c1=ROOT.TCanvas()
    h.Draw("COLZ")
    exclusion.Draw("SAME")
    orig_exclusion.Draw("SAME")
    c1.Print("../plots/%s_%s.png" %(ana,topo))
    f2.Close()
    return "%s.root" %topo
  return h
