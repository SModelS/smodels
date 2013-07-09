import set_path
import ROOT,sys, os, tempfile
from Experiment import SMSResults, SMSInterpolation, LimitGetter, SMSAnalysisFactory
from Theory import SMSAnalysis, LHEDecomposer, SLHATools, SLHADecomposer, XSecComputer
from Tools import ROOTTools
from Tools.PhysicsUnits import addunit, rmvunit

pid_dic = {"T1": [[1000021], [1000022]],
           "T1tttt": [[1000021], [1000022]],
           "T2": [[1000001, 1000002, 1000003, 1000004, 1000005, 2000001, 2000002, 2000003, 2000004, 2000005], [1000022]],
           "T2bb": [[1000005],[1000022]],
           "T2tt": [[1000006],[1000022]],
           "TChiWZ": [[1000024, 1000023], [1000022]],
           "TChiChipmSlepL": [[1000024, 1000023], [1000011, 1000012, 1000013, 1000014, 1000015, 1000016], [1000022]],
           "TChiChipmSlepStau": [[1000024, 1000023], [1000011, 1000012, 1000013, 1000014, 1000015, 1000016], [1000022]],#??set all slepton masses here
           "TChiChipmStauStau": [[1000024, 1000023], [1000011, 1000012, 1000013, 1000014, 1000015, 1000016], [1000022]] }#??set all slepton masses here

diff_dic = {"T1tttt": 360., "T2bb": 5.}


def recreateHist(ana,topo,mz=None,axes=None, run='',line=False,tev=8,nevents=10000,binsize=None, fromSlha=True):
  """ recreate ROOT TH2F histogram of a given analysis and topology
        needs mz for a topology with intermediate mass,
        axes, for histograms with axes different from M1-M0
        if line=True is selected, produce rootfile with produced histograms
        and line, return filename"""

#  lhefile = "../lhe/%s_1.lhe" %topo
#  if fromSlha:
#    import os
#    os.system("cp ../slha/%s.slha fort.61" %topo)
#    os.system("../pythia_lhe < pyIn.dat") #run with 20 events
#    lhefile = "fort.68" 
#  topoList=LHEDecomposer.decompose ( lhefile, {})#create default topologylist with only topo (is list of GTop objects)


  Tmp = tempfile.mkdtemp()


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
    hUL = ROOT.TH2F('hUL','hUL',bx,xmin,xmin+bx*bwx,by,ymin,ymin+by*bwy)
    hTheory = ROOT.TH2F('hTheory','hTheory',bx,xmin,xmin+bx*bwx,by,ymin,ymin+by*bwy)
    hLine = ROOT.TH2F('hLine','hLine',bx,xmin,xmin+bx*bwx,by,ymin,ymin+by*bwy)
    prod_mode = topo
    if topo=="T1tttt" or topo=="T1bbbb": prod_mode = "T1"
    if topo=="TChiChipmSlepL" or topo=="TChiChipmSlepStau" or topo=="TChiChipmStauStau" or topo=="TChiChipmStauL" or topo=="TChiWZ": prod_mode = "TChiWZ_Wino"
    if topo=="T2": prod_mode = "T24sq"
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
    while y<ymax and y<x-diff_dic[topo]:

      print x, y

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

      ana_obj = SMSAnalysisFactory.load( anas=ana, topos=topo ) #create list of analysis objects, but give one analysis, one topo, so list has one entry only!

      slha_dic = {}
      for ii in range(len(pid_dic[topo])):
        for ent in pid_dic[topo][ii]:
          slha_dic[ent]=massv[ii]

      sigmacut = addunit(0.1,'fb')
      slha_name = SLHATools.createSLHAFile(topo, slha_dic)
      xsec_dic = XSecComputer.compute(nevents, slha_name, datadir = Tmp)
      XSec = xsec_dic.crossSections()
      topoList = SLHADecomposer.decompose(slha_name, XSec, sigmacut)
      os.unlink(slha_name)

      ana_obj[0].add( topoList ) # the list ana_obj has only one entry here!, all EElements consisten with the constraints for this topology are added to the EAnalysis object
      lims=LimitGetter.limit(ana_obj[0], addTheoryPrediction=False)
      v=None
      if lims: v=rmvunit(lims[0]['ul'],'fb')

#      v=rmvunit(LimitGetter.GetPlotLimit([massv,massv],[topo,[ana]],a)[0][1],'fb')

      if v: h.Fill(x,y,v)

      if line:
        sumTheory=0.
        if v and v<rXs.GetBinContent(rXs.FindBin(x)):
          hL.Fill(x,y,0)
        else: hL.Fill(x,y)
        if ana_obj[0].computeTheoryPredictions():
          hUL.Fill(x,y,rmvunit(ana_obj[0].ResultList[0].explimit,'fb'))
          for values in ana_obj[0].ResultList[0].result_dic.values():
            hTheory.Fill(x,y,rmvunit(values['8 TeV (NLL)'],'fb'))
            sumTheory+=rmvunit(values['8 TeV (NLL)'],'fb')
          if rmvunit(ana_obj[0].ResultList[0].explimit,'fb') and rmvunit(ana_obj[0].ResultList[0].explimit,'fb')<sumTheory:
            hLine.Fill(x,y,0)
          else: hLine.Fill(x,y)
        else: hLine.Fill(x,y)

      y+=bwy

    while y < ymax:
      hL.Fill(x,y)
      hLine.Fill(x,y)
      y+=bwy

    y=ymin+bwy/2
    x+=bwx

  XSecComputer.clean(Tmp)

  if line:
    f.Close()
    rootname_out="%s_%devts.root" % (topo,nevents)
    if binsize:rootname_out="%s_%devts_%sGeVbin.root" % (topo,nevents,binsize)
    f1=ROOT.TFile(rootname_out,"recreate")
    h.Write()
    hL.Write()
    hUL.Write()
    hTheory.Write()
    hLine.Write()
    exclusion=ROOTTools.getTGraphFromContour(hL)
    exclusion.Write()
    exclusionER=ROOTTools.getTGraphFromContour(hLine)
    exclusionER.Write()
    f1.Close()
    return rootname_out
  return h
