import set_path
import ROOT,sys
from Experiment import SMSResults, SMSInterpolation, LimitGetter, SMSAnalysisFactory
from Theory import SMSAnalysis, LHEDecomposer
from Tools import ROOTTools
from Tools.PhysicsUnits import addunit, rmvunit

def recreateHist(ana,topo,mz=None,axes=None, run='',line=False,tev=8,nevents=10000,binsize=None, fromSlha=True):
  """ recreate ROOT TH2F histogram of a given analysis and topology
        needs mz for a topology with intermediate mass,
        axes, for histograms with axes different from M1-M0
        if line=True is selected, produce rootfile with produced histograms
        and line, return filename"""

  lhefile = "../lhe/%s_1.lhe" %topo
  if fromSlha:
    import os
    os.system("cp ../slha/%s.slha fort.61" %topo)
    os.system("../pythia_lhe < pyIn.dat") #run with 20 events
    lhefile = "fort.68" 
  topoList=LHEDecomposer.decompose ( lhefile, {})#create default topologylist with only topo


  toponame=topo
  run1=SMSResults.getRun(ana)
  if mz: toponame=SMSInterpolation.gethistname(topo,mz)
  xmin=rmvunit(SMSResults.getLowX(ana,toponame),'GeV')
  ymin=rmvunit(SMSResults.getLowY(ana,toponame),'GeV')
  xmax=rmvunit(SMSResults.getUpX(ana,toponame),'GeV')
  ymax=rmvunit(SMSResults.getUpY(ana,toponame),'GeV')
  bwx=float(rmvunit(SMSResults.getBinWidthX(ana,toponame),'GeV'))
  bwy=float(rmvunit(SMSResults.getBinWidthY(ana,toponame),'GeV'))
  if binsize:
    bwx=float(binsize)
    bwy=float(binsize)
  print bwx, bwy, xmin, xmax, ymin, ymax
  bx=int((xmax-xmin)/bwx)
  by=int((ymax-ymin)/bwy)

  h = ROOT.TH2F('h','h',bx,xmin,xmin+bx*bwx,by,ymin,ymin+by*bwy)

  if line:
    hL = ROOT.TH2F('hL','hL',bx,xmin,xmin+bx*bwx,by,ymin,ymin+by*bwy)
    prod_mode = topo
    if topo=="T1tttt" or topo=="T1bbbb": prod_mode = "T1"
    if topo=="TChiChipmSlepL" or topo=="TChiChipmSlepStau" or topo=="TChiChipmStauStau": prod_mode = "TChiWZ_Wino"
    rootname = "../data/%s_%devts.root" %(prod_mode,nevents)
    if binsize: rootname =  "../data/%s_%devts_%sGeVbin.root" %(prod_mode,nevents,binsize)
    f=ROOT.TFile(rootname)
    if tev==8: rXsName="hist8"
    if tev==7:rXsName="hist7"
    rXs=f.Get(rXsName)
    print rootname

  if line and not rXs:
    print "No refXSec histogram found"
    return None

  x=xmin+bwx/2
  y=ymin+bwy/2

#  a=SMSAnalysis.EAnalysis()
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

      ana_obj = SMSAnalysisFactory.load( anas=ana, topos=topo )
      for idx in range(len(massv)):
        massv[idx]=addunit(massv[idx],'GeV')
      for eachtopo in topoList:
        for eachEl in eachtopo.ElList:
          if len(eachEl.B[0].masses)==len(eachEl.B[1].masses)==len(massv):
            eachEl.B[0].masses=massv
            eachEl.B[1].masses=massv #for now only equal branches
      ana_obj[0].add( topoList )
#      print topos, topos[0].ElList[0].B[0].masses, ana_obj[0].Top
#      ana_obj[0].Top.ElList[0].B[0].masses=massv
#      ana_obj[0].Top.ElList[0].B[1].masses=massv
      lims=LimitGetter.limit(ana_obj[0], addTheoryPrediction=False)
      v=None
      if lims: v=rmvunit(lims[0]['ul'],'fb')

#      v=rmvunit(LimitGetter.GetPlotLimit([massv,massv],[topo,[ana]],a)[0][1],'fb')

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
    if binsize:rootname_out="%s_%devts_%sGeVbin.root" % (topo,nevents,binsize)
    f1=ROOT.TFile(rootname_out,"recreate")
    h.Write()
    hL.Write()
    exclusion=ROOTTools.getTGraphFromContour(hL)
    exclusion.Write()
    f1.Close()
    return rootname_out
  return h
