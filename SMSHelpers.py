#!/usr/bin/python

""" some helper functions that are not to be used by the end user """

import os

if 'lessa' in os.getcwd():
    Base = "/home/lessa/SMS_database/"
else:
    Base = "/afs/hephy.at/user/w/walten/public/sms/"

runs=[ "8TeV", "2012", "2011", "ATLAS8TeV" ]
## runs=[ "2012" ]

verbose=True

## track the open root files
openFiles={}

useUnits=False

def addunit ( value, unitstring ):
  """ a function that can add units to values, but also
      makes it easy to turn this functionality off, in case "units" isnt installed
      """
  if not useUnits: return value
  if useUnits: 
    import SMSUnits as u
    import types as t
    if type(value)!=t.FloatType and type(value)!=t.IntType:
      return value
    if unitstring=="GeV":
      return value * u.GeV
    if unitstring=="TeV":
      return value * u.TeV
    if unitstring=="fb":
      return value * u.fb
    if unitstring=="pb":
      return value * u.pb
    if unitstring=="fb-1":
      return value / u.fb
    print "[SMSHelpers.py] Warning: dont know what to do with unit",unitstring
  return value


def rmvunit ( value, unitstring ):
  """ a function that can remove units from values, but also
      makes it easy to turn this functionality off, in case "units" isnt installed
      """
  if not useUnits: return value
  if useUnits: 
    import SMSUnits as u
    import types as t
    if type(value) != type(1.*u.GeV):
      return value
    if unitstring=="GeV":
      return value.asNumber(u.GeV)
    if unitstring=="TeV":
      return value.asNumber(u.TeV)
    if unitstring=="fb":
      return value.asNumber(u.fb)
    if unitstring=="pb":
      return value.asNumber(u.pb)
    if unitstring=="fb-1":
      return value.asNumber(1/u.fb)

    print "[SMSHelpers.py] Warning: dont know what to do with unit",unitstring
    return value

    

def close():
  """ close all open files """
  for (name,tfile) in openFiles.items():
    tfile.Close()
  openFiles={}

logLevel= [ "error" ]

def log ( text, level="error" ):
  if level in logLevel:
    print "[SMSHelpers.py] %s: %s" % ( level, text )
  return

def getRun ( analysis, run=None ):
  """ search for an analysis, return the run,
      or return None if not found anywhere.
      if a specific run is given, then just check
      if we really have results. """
  if run:
    if os.path.exists ( "%s/%s/%s" % ( Base, run, analysis ) ):
      return run
    else:
      log ( "dont know about %s/%s" % (run,analysis), "warning" )
      return None
  for trun in runs:
    if os.path.exists ( "%s/%s/%s" % ( Base, trun, analysis ) ):
      return trun
  return None

def parseMetaInfo ( analysis, run ):
  """ get all the meta information for a given analysis/run pair """
  info="%s/%s/%s/info.txt" % ( Base, run, analysis )
  ret={}
  if not os.path.exists ( info ):
    log ( "cannot find %s" % info, "warning" )
    return ret
  f=open(info)
  lines=f.readlines()
  f.close()
  for line in lines:
    tokens=line.split(":",1)
    if not len(tokens)==2:
      log ( "cannot parse this line (1): %s in %s" % (line, info) )
      continue
    if tokens[0]=="exclusions":
      # we treat these separately
      continue
    ret[tokens[0]]=tokens[1]
  return ret

def motherParticleExclusions ( analysis, run ):
  """ get all the exclusion numbers for a given analysis/run pair """
  info="%s/%s/%s/info.txt" % ( Base, run, analysis )
  ret={}
  if not os.path.exists ( info ):
    log ("cannot find %s" % info )
    return ret
  f=open(info)
  lines=f.readlines()
  f.close()
  for line in lines:
    tokens=line.split(":",1)
    if not len(tokens)==2:
      log ( "cannot parse this line (2): %s in %s" % (line, info) )
      continue
    if tokens[0]!="exclusions":
      # we're only interested in the exclusions
      continue
    excl=tokens[1]
    while excl[0]==" ": excl=excl[1:]
    if excl[-1]=='\n': excl=excl[:-1]
    excls=excl.split(" ")
    if len(excls)<2:
      log ( "cannot parse this line (3): %s in %s" % (line, info) )
      continue
    values=map(int,excls[1:])
    if len(values)==2:
      values.insert(0,0)
    ret[excls[0]]=values
  return ret

def getLines ( analysis, run, label="condition" ):
  """ get all the conditions for a analysis/run pair """
  info="%s/%s/%s/info.txt" % ( Base, run, analysis )
  ret={}
  if not os.path.exists ( info ):
    log ("cannot find %s" % info )
    return ret
  f=open(info)
  lines=f.readlines()
  f.close()
  for line in lines:
    tokens=line.split(":",1)
    if not len(tokens)==2:
      log ( "cannot parse this line (2): %s in %s" % (line, info) )
      continue
    if tokens[0]!=label:
      # we're only interested in the conditions
      continue
    excl=tokens[1]
    while excl[0]==" ": excl=excl[1:]
    if excl[-1]=='\n': excl=excl[:-1]
    keyvalue=excl.split(" -> ")
    if len(keyvalue)!=2:
      log ( "cannot parse the following line: %s" % keyvalue )
    ret[ keyvalue[0] ] = keyvalue[1]
    # ret.append(excl)
  return ret

def conditions ( analysis, run ):
  """ get all the conditions for a analysis/run pair """
  return getLines ( analysis, run, "condition" )

def constraints ( analysis, run ):
  """ get all the conditions for a analysis/run pair """
  return getLines ( analysis, run, "constraint" )


def getUpperLimitHisto ( analysis, topo, run ):
  import ROOT
  rootfile="%s/%s/%s/sms.root" % ( Base, run, analysis )
  if not os.path.exists(rootfile):
    log("root file %s doesnt exist" % rootfile )
    return None
  f=None
  if openFiles.has_key ( rootfile ):
    f=openFiles[rootfile]
  else:
    f=ROOT.TFile(rootfile)
    if not f or not f.IsOpen():
      log("root file %s cannot be opened" % rootfile )
      return None
    openFiles[rootfile]=f
  histo=f.Get("limit_%s" % topo )
  if not histo:
    log("histogram %s not found in %s" % ( topo, rootfile ))
    return None
  return histo

def getUpperLimitAtPoint ( histo, mx, my ):
  """ given already a histogram, we get the upper limit for mx and my,
      we return None if we're out of bounds """

  if rmvunit(my,'GeV')==None:
    log ( "inconsistent mx/my, mx=None, my=%s" % my )
    return None
  if not histo: return 'no histogram'
  xmax=histo.GetXaxis().GetNbins()
  xbin=histo.GetXaxis().FindBin(rmvunit(mx,'GeV'))
  if xbin==0 or xbin>xmax:
    return None
  ymax=histo.GetYaxis().GetNbins()
  ybin=histo.GetYaxis().FindBin(rmvunit(my,'GeV'))
  if ybin==0 or ybin>ymax:
    return None
  c=histo.GetBinContent(xbin,ybin)

  if c == 0.: c = None
  
  return c

def getUpperLimitPng(analysis,topo,run):
  pngfile="%s/%s/%s/results/h_limit_%s_%s.png" % ( Base, run, analysis, topo, analysis )
  if not os.path.exists(pngfile):
    log("png file %s doesnt exist" % rootfile )
    return None
  return pngfile

def getEfficiencyHisto ( analysis, topo, run ):
  import ROOT
  rootfile="%s/%s/%s/sms.root" % ( Base, run, analysis )
  if not os.path.exists(rootfile):
    log("root file %s doesnt exist" % rootfile )
    return None
  f=None
  if openFiles.has_key ( rootfile ):
    f=openFiles[rootfile]
  else:
    f=ROOT.TFile(rootfile)
    if not f or not f.IsOpen():
      log("root file %s cannot be opened" % rootfile )
      return None
    openFiles[rootfile]=f
  histo=f.Get("efficiency_%s" % topo )
  if not histo:
    log("histogram %s not found in %s" % ( topo, rootfile ))
    return None
  return histo

def getEfficiencyAtPoint ( histo, mx, my ):
  """ given already a histogram, we get the efficiency for mx and my,
      we return None if we're out of bounds """
  if my==None:
    log ( "inconsistent mx/my, mx=None, my=%s" % my )
    return None
  if not histo: return None
  xmax=histo.GetXaxis().GetNbins()
  xbin=histo.GetXaxis().FindBin(mx)
  if xbin==0 or xbin>xmax:
    return None
  ymax=histo.GetYaxis().GetNbins()
  ybin=histo.GetYaxis().FindBin(my)
  if ybin==0 or ybin>ymax:
    return None
  c=histo.GetBinContent(xbin,ybin)
  if not c>=0. or not c<=1.:
    print "[SMSHelpers.py] warning efficiency not between 0 and 1? c=",c
  return c

def getErrorMessage ( histo, mx, my ):
  """ return an explication for why we have no result at mx, my """
  if not histo:
    return [ 1, "No 2d histogram available" ]
  xbin=histo.GetXaxis().FindBin(mx)
  titlex=histo.GetXaxis().GetTitle().replace("NLO","").replace( "NLL", "").replace("[GeV]","").replace("_","").replace("7TeV","").replace("8TeV","").replace(" mass ","")
  titley=histo.GetYaxis().GetTitle().replace("NLO","").replace( "NLL", "").replace("[GeV]","").replace("_","").replace("7TeV","").replace("8TeV","").replace(" mass ","")
  if xbin==0:
    return [ 2, "m("+titlex+") too low" ]
  xmax=histo.GetXaxis().GetNbins()
  if xbin>xmax:
    return [ 3, "m("+titlex+") too high" ]
  ybin=histo.GetYaxis().FindBin(my)
  if ybin==0:
    return [ 4, "m("+titley+") too low" ]
  ymax=histo.GetYaxis().GetNbins()
  if ybin>ymax:
    return [ 5, "m("+titley+") too high" ]
  c=histo.GetBinContent(xbin,ybin)
  if c==0.:
    return [ 6, "out of histogram bounds" ]
  return [ 0, str(c) ]

def hasMetaInfoField ( analysis, field, run=None ):
  run=getRun ( analysis, run )
  metainfo=parseMetaInfo ( analysis, run )
  return metainfo.has_key ( field )

def getMetaInfoField ( analysis, field, run=None, complain=True ):
  """ get one specific entry of the meta info """
  run=getRun ( analysis, run )
  metainfo=parseMetaInfo ( analysis, run )
  if not metainfo.has_key ( field ):
    if complain:
      log ( "%s/%s doesnt have a ``%s'' field." % ( run, analysis, field ) )
    return None
  f=metainfo[field]
  while f[0]==' ': f=f[1:]
  if f[-1]=='\n': f=f[:-1]
  return f

dicparticle = { 
   "T1": "#tilde{g}",
   "T1bbbb": "#tilde{g}",
   "T1tttt": "#tilde{g}",
   "T1lg": "#tilde{g}",
   "T1gg": "#tilde{g}",
   "T5gg": "#tilde{g}",
   "T5Wg": "#tilde{g}",
   "T2": "#tilde{q}",
   "T2bb": "#tilde{b}",
   "T2tt": "#tilde{t}",
   "T2ttww": "#tilde{b}",
   "T3w": "#tilde{g}",
   "T3lh": "#tilde{g}",
   "T5zz": "#tilde{g}",
   "T5zzInc": "#tilde{g}", 
   #"T5zzh": "#tilde{g}", 
   #"T5zzl": "#tilde{g}", 
   "T5lnu": "#tilde{g}",
   "T5zzgmsb": "#tilde{g}", 
   "TChiwz": "#tilde{#chi}^{0}_{2}#tilde{#chi}^{#pm}_{1}",
   #"TChiwz": "#tilde{#chi}^{0}",
   "TChizz": "#tilde{#chi}^{0}_{2}",
   "TChiSlep": "#tilde{#chi}^{0}_{2}",
   "TChiSlepSlep": "#tilde{#chi}^{0}_{2}#tilde{#chi}^{#pm}_{1}",
   "TChiNuSlep": "#tilde{#chi}^{0}_{2}",
}

def hasDictionary ( analysis, run=None ):
  """ are the upper limits available in dictionary format? """
  hd=getMetaInfoField ( analysis, "dictionary", run, complain=False )
  if hd:
    if hd in [ '1', 'True' ]: return True
  return False

def getUpperLimitDictionary ( analysis, topo, run ):
  dictfile="%s/%s/%s/sms.py" % ( Base, run, analysis )
  if not os.path.exists(dictfile):
    log("dictionary file %s doesnt exist" % dictfile )
    return None
  Globals={}
  execfile(dictfile,Globals)
  Dict=Globals["Dict"]
  if not Dict.has_key ( topo ):
    log("dictionary doesnt have topology"+topo )
    return None
  return Dict[topo]
