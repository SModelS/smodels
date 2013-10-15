#!/usr/bin/env python
import os,sys
from ROOT import TTree,TColor,TCanvas,TF1,TGraph


Var_dic = {'So4_mass': "m_{#tilde{#chi}^{0}_{4}} (GeV)", 'pSI': "#sigma(p,SI) (pb)", 'H3_mass': "m_{A} (GeV)", 'SmL_mass': "m_{#tilde{e}_{L}} (GeV)", 'Mr1': "#tilde{m}_{#tilde{e}_{R}} (GeV)", 'pSD': "#sigma(p,SD) (pb)", 'Mr3': "#tilde{m}_{#tilde{#tau}_{R}} (GeV)", 'Ml2': "#tilde{m}_{#tilde{#mu}_{L}} (GeV)", 'Ml3': "#tilde{m}_{#tilde{#tau}_{L}} (GeV)", 'MG1': "M_{1} (GeV)", 'Ml1':  "#tilde{m}_{#tilde{e}_{L}} (GeV)", 'SmR_mass': "m_{#tilde{e}_{LR} (GeV)", 'bsmumu': "BR(B_{s} #rightarrow #mu #mu)", 'Bsmumu': "BR(B_{s} #rightarrow #mu #mu)", 'gmuon': "a_{#mu}",'Gmuon': "#Deltaa_{#mu}", 'h_mass': "m_{h} (GeV)", 'Mq1': "#tilde{m}_{#tilde{Q}}(1)} (GeV)", 'Mq3':  "#tilde{m}_{#tilde{Q}}(3)} (GeV)", 'Mq2':  "#tilde{m}_{#tilde{Q}}(2)} (GeV)", 'Am': 'A_{#tilde{#mu}} (GeV)', 'Md2':  "#tilde{m}_{#tilde{d}}(2)} (GeV)", 'Md3': "#tilde{m}_{#tilde{d}}(3)} (GeV)", 'Md1': "#tilde{m}_{#tilde{d}}(1)} (GeV)", 'Omega': "#Omega_{#tilde{#chi}_{1}^{0}} h^{2}", 'MH3': "M_{A} (GeV)", 'So2_mass': "m_{#tilde{#chi}^{0}_{2}} (GeV)", 'Mu1': "#tilde{m}_{#tilde{u}}(1)} (GeV)", 'Mu3': "#tilde{m}_{#tilde{u}}(3)} (GeV)", 'Mu2': "#tilde{m}_{#tilde{u}}(2)} (GeV)", 'Mtp': "m_{top} (GeV)", 'Ab': 'A_{b} (GeV)', 'So1_mass': "m_{#tilde{#chi}^{0}_{1}} (GeV)", 'MG3': "M_{3} (GeV)", 'Snm_mass': "m_{#tilde{#nu}_{#mu}} (GeV)", 'MG2': "M_{2} (GeV)", 'tb': "tan#beta", 'So3_mass': "m_{#tilde{#chi}^{0}_{3}} (GeV)", 'Al': 'A_{l} (GeV)', 'mu': '#mu (GeV)', 'C1p_mass': "m_{#tilde{#chi}^{#pm}_{1}} (GeV)", 'At': 'A_{t} (GeV)', 'Mr2': "#tilde{m}_{#tilde{#mu}_{R}} (GeV)", 'Bsgm': "BR(b #rightarrow s #gamma)", 'bsgmm': "BR(b #rightarrow s #gamma)", 'C2p_mass': "m_{#tilde{#chi}^{#pm}_{2}} (GeV)","epRat" : "#epsilon^{MSSM}/#epsilon^{SM}", "epMSSM" : "#epsilon^{MSSM}", "epSM" : "#epsilon^{SM}", "Smu1_mass": "m_{#tilde{#mu}_{1}} (GeV)", "pSI_eff" : "#sigma(p,SI)^{eff} (pb)", "mC1-mN1" : "m_{#tilde{#chi}^{#pm}_{1}}-m_{#tilde{#chi}^{0}_{1}} (GeV)","mSmu1-mN1" : "m_{#tilde{#mu}_{1}}-m_{#tilde{#chi}^{0}_{1}} (GeV)","mSlep-mN1" : "m_{#tilde{l}}-m_{#tilde{#chi}^{0}_{1}} (GeV)", "mSl1-mN1" : "m_{#tilde{#tau}_{1}}-m_{#tilde{#chi}^{0}_{1}} (GeV)"}  

Exp_dic = {"Smu1_mass" : "min(TREENAME.SmL_mass,TREENAME.SmR_mass)","Smu2_mass" : "max(TREENAME.SmL_mass,TREENAME.SmR_mass)","pSI_eff" : "TREENAME.pSI*TREENAME.Omega/0.11", "mC1-mN1" : "TREENAME.C1p_mass-TREENAME.So1_mass","mSmu1-mN1" : "min(TREENAME.SmL_mass,TREENAME.SmR_mass)-abs(TREENAME.So1_mass)","mSlep-mN1" : "min(TREENAME.SmL_mass,TREENAME.SmR_mass,TREENAME.Sl1_mass)-abs(TREENAME.So1_mass)","mSl1-mN1" : "TREENAME.Sl1_mass-abs(TREENAME.So1_mass)"}



def infiles(argv,SMS="smsComp",tryConst=True):
  if len(argv) > 1:
    filename = argv[1]
  else:  
    filename = "gmu_scanCT.root"
    
  if SMS and os.path.isfile(filename.replace(".root","_"+SMS+".root")):
    hasSMS = filename.replace(".root","_"+SMS+".root")
  else:
    hasSMS = False
  if tryConst and os.path.isfile(filename.replace(".root","_const.root")):
    hasConst = filename.replace(".root","_const.root")
  else:
    hasConst = False  

  return filename,hasSMS,hasConst

def getTree(gDirectory,hasSMS,hasConst,verbose=False):
    
#Get Trees:
  objs =  gDirectory.GetListOfKeys()
  trees = []
  for ob in objs:
    Tob = ob.ReadObj()
    if type(Tob) == type(TTree()):
      trees.append(Tob)

  if len(trees) <= 0:
    print "No trees found in file"
    sys.exit()
  elif len(trees) > 1:
    pstr = "More than one tree found, specify which one to use ("
    for tree in trees:
      pstr += tree.GetName()+","
    pstr += ")"
    usetree = raw_input(pstr)
    tree = gDirectory.Get(usetree)
  else:
    tree = trees[0]
    
  if hasSMS:
    print "Adding friend",hasSMS
    tree.AddFriend("tree",hasSMS)
  if hasConst:
    print "Adding friend",hasConst
    tree.AddFriend("tree",hasConst)
    
#Get Branches and Leaves (variables):
  branches = tree.GetListOfBranches()
  leaves = tree.GetListOfLeaves()
  lnames = [leaf.GetName() for leaf in leaves]
  
#Get Input parameters:
  inputpars = {}
  for br in branches:  
    if br.GetName().lower() != "input": continue
    for leaf in br.GetListOfLeaves():
      lname = leaf.GetName()
      inputpars.update({lname : [tree.GetMinimum(lname),tree.GetMaximum(lname)]})
      
  if verbose:      
    print "Nevts = ",tree.GetEntries()
    print "Input parameters and range:"
    for key in inputpars.keys():
      print key,inputpars[key]      

  return tree,lnames,inputpars

#Use Var_dic to create latex-compliant names for the leaves.
#If Var_dic does not contain the leaf name, use the own name:
def GetVarNames(allvars):
  varnames = {}
  for vname in allvars:
    if Var_dic.has_key(vname):
      varnames.update({vname : Var_dic[vname]})
    else:
      varnames.update({vname : vname})
    
  return varnames

#Check if variables in AllowedRange.keys() have their values in the allowed range
def GetExcluded(AllowedRange,tree,hasConst):
  import bounds

  excluded = {}

  if not hasConst and "Xe100" in AllowedRange.keys():
    excluded["Xe100"] = bounds.xexcl(tree.So1_mass,tree.pSI,tree.Omega)
    if excluded["Xe100"] <= AllowedRange["Xe100"][1] and excluded["Xe100"] >= AllowedRange["Xe100"][0]:
      excluded["Xe100"] = False
    else:
      excluded["Xe100"] = True

  for exc in AllowedRange.keys():    
    try:
      if getattr(tree,exc) <= AllowedRange[exc][1] and getattr(tree,exc) >= AllowedRange[exc][0]:
        excluded[exc] = False
      else:
        excluded[exc] = True
    except:
      pass


  
  return excluded  
    
    
def GetValue(tree,x):

  if hasattr(tree,x): return getattr(tree,x)
  if x in Exp_dic.keys():
    exp = Exp_dic[x]
    exp = exp.replace("TREENAME","tree")
    return eval(exp)
  else:
    print "[GetValue]: Unknown variable",x
    sys.exit()
  
  
    
def Default(obj,Type):
  
  if Type == "TCanvas":
    obj.SetLeftMargin(0.1097891)
    obj.SetRightMargin(0.02700422)
    obj.SetTopMargin(0.02796053)
    obj.SetBottomMargin(0.14796053)
    obj.SetFillColor(0)
    obj.SetBorderSize(0)
    obj.SetFrameBorderMode(0)
  elif "TGraph" in Type or "TH" in Type:
    obj.GetYaxis().SetTitleFont(132)
    obj.GetYaxis().SetTitleSize(0.065)
    obj.GetYaxis().CenterTitle(True)
    obj.GetYaxis().SetTitleOffset(0.9)
    obj.GetXaxis().SetTitleFont(52)
    obj.GetXaxis().SetTitleSize(0.065)
    obj.GetXaxis().CenterTitle(True)
    obj.GetXaxis().SetTitleOffset(1.0)
    obj.GetYaxis().SetLabelFont(132)
    obj.GetXaxis().SetLabelFont(132)
    obj.GetYaxis().SetLabelSize(0.05)
    obj.GetXaxis().SetLabelSize(0.05)
    if "TGraph2D" in Type or "TH2" in Type:
      obj.GetZaxis().SetTitleFont(132)
      obj.GetZaxis().SetTitleSize(0.06)
      obj.GetZaxis().CenterTitle(True)
      obj.GetZaxis().SetTitleOffset(0.7)
      obj.GetZaxis().SetLabelFont(132)
      obj.GetZaxis().SetLabelSize(0.05)
  elif "Leg" in Type:
   obj.SetBorderSize(1)
   obj.SetTextFont(132)
   obj.SetTextSize(0.05)
   obj.SetLineColor(1)
   obj.SetLineStyle(1)
   obj.SetLineWidth(1)
   obj.SetFillColor(0)
   obj.SetFillStyle(1001)


def Print(canvas,prefix,hasSMS):
  filename = prefix
  if hasSMS:
    addname = hasSMS.rstrip(".root")
    while addname.count("_") > 0:
      addname = addname[addname.index("_")+1:]
    filename += "_"+addname
  filename += ".eps"
  canvas.Print(filename)

def printInput(tree,inputpars):
  
  for key in inputpars.keys():
    print key,getattr(tree,key)    
    
    
def set_palette(gStyle,name="none", ncontours=999):
    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""
    
    from array import array

    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)

    npoints = len(s)
    TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    gStyle.SetNumberContours(ncontours)

def getContours(gROOT,vals,hist,fit=None):
    """Get contour graphs from a histogram with contour values vals"""
    from array import array

    res = {}
    canv = TCanvas("getContour_canvas","getContour_canvas",0,0,500,500)
    canv.cd()
    h2 = hist.Clone()
    contours = array('d',vals)
    h2.SetContour(len(vals),contours)
    h2.Draw("CONT Z LIST")
    canv.Update()
    conts = gROOT.GetListOfSpecials().FindObject("contours")
    for ival,val in enumerate(vals):
        contLevel = conts.At(ival)
        res[val] = contLevel.First()
        for cont in contLevel:
            if cont.GetN() > res[val].GetN(): res[val] = cont   #Get the contour wiht highest number of points 

    if fit:
        f1 = TF1("f1",fit,0.,1000.)
        for val in vals:
            curv = res[val].Clone()
            curv.Sort()
            xmin = curv.GetX()[0]
            xmax = curv.GetX()[curv.GetN()-1]
            f1.SetMinimum(xmin)
            f1.SetMaximum(xmax)
            curv.Fit(f1)
            fit = TGraph(curv.FindObject("f1"))
            res[val] = fit.Clone()

    canv.Clear()

    return res

