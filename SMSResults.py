#!/usr/bin/python

""" centralized facility to access the SMS results """

import SMSHelpers
from SMSHelpers import addunit, rmvunit

def setLogLevel ( l=[ "error" ] ):
  """ defines what is written out, and what isnt """
  SMSHelpers.logLevel=l


def useUnits ( b=True ):
  SMSHelpers.useUnits = b

def ResultsForSqrts ( sqrts ):
  import SMSResultsCollector
  """ If this is called, only results for a
      center of mass energy of sqrts in TeV will be retrieved.
      if sqrts is None or 0, then all results will be made available """
  if sqrts==7:
    SMSHelpers.runs=[ "2012", "2011" ]
    SMSResultsCollector.alldirectories=[ "2011", "2012" ]
    return
  if sqrts==8:
    SMSHelpers.runs=[ "8TeV", "ATLAS8TeV" ]
    SMSResultsCollector.alldirectories=[ "8TeV", "ATLAS8TeV" ]
    return
  if sqrts==0 or sqrts==None:
    SMSHelpers.runs=[ "8TeV", "2012", "2011", "ATLAS8TeV" ]
    SMSResultsCollector.alldirectories=[ "8TeV", "2012", "2011", "ATLAS8TeV" ]
    return
  print "[SMSResults.py:ResultsForSqrts] error: dont have any results for",sqrts

def verbosity ( level="error" ):
  if level=="error":
    SMSHelpers.verbose=False
  else:
    SMSHelpers.verbose=True

def getExclusion ( analysis, topo, run=None ):
  """ get the exclusions on the mother particles,
      for the summary plots """
  run=SMSHelpers.getRun ( analysis, run )
  ex=SMSHelpers.motherParticleExclusions ( analysis, run )
  for t in SMSHelpers.getPotentialNames ( topo ):
    if ex.has_key ( t ): return ex[t]
  return None

def exists ( analysis, topo, run=None ):
  """ does a result for analysis/topo/run exist? """
  run=SMSHelpers.getRun ( analysis, run )
  histo=SMSHelpers.getUpperLimitFromHisto ( analysis, topo, run )
  return histo!=None

def getBinWidthX ( analysis, topo, run=None ):
  """ get the bin width in X """
  run=SMSHelpers.getRun ( analysis, run )
  if SMSHelpers.hasDictionary ( analysis, run ):
    return None
  histo=SMSHelpers.getUpperLimitFromHisto ( analysis, topo, run )
  if not histo: return None
  w=histo.GetXaxis().GetBinWidth(1)
  return SMSHelpers.addunit ( w, "GeV" )

def getLowX ( analysis, topo, run=None ):
  """ get the lower edge of the x axis """
  run=SMSHelpers.getRun ( analysis, run )
  if SMSHelpers.hasDictionary ( analysis, run ):
    print "[SMSResults.py] implement this function"
    return None
  histo=SMSHelpers.getUpperLimitFromHisto ( analysis, topo, run )
  if not histo: return None
  w=histo.GetXaxis().GetXmin()
  return SMSHelpers.addunit ( w, "GeV" )

def getLowY ( analysis, topo, run=None ):
  """ get the lower edge of the y axis """
  run=SMSHelpers.getRun ( analysis, run )
  if SMSHelpers.hasDictionary ( analysis, run ):
    print "[SMSResults.py] implement this function"
    return None
  histo=SMSHelpers.getUpperLimitFromHisto ( analysis, topo, run )
  if not histo: return None
  w=histo.GetYaxis().GetXmin()
  return SMSHelpers.addunit ( w, "GeV" )

def getBinWidthY ( analysis, topo, run=None ):
  """ get the bin width in Y """
  run=SMSHelpers.getRun ( analysis, run )
  if SMSHelpers.hasDictionary ( analysis, run ):
    return None
  histo=SMSHelpers.getUpperLimitFromHisto ( analysis, topo, run )
  if not histo: return None
  w=histo.GetYaxis().GetBinWidth(1)
  return SMSHelpers.addunit ( w, "GeV" )

def getExclusionLine(topo,ana,expected=False,plusminussigma=0,extendedinfo=True,xvalue=None,factor=1.0):
  """ get the exclusion line, as a TGraph """
  if xvalue==None: xvalue=''
  import SMSResultsCollector
  ex=SMSResultsCollector.exclusionline(topo,ana,xvalue,expected,plusminussigma)
  return ex

def getTopologies ( analysis, run=None ):
  """ return all topologies that this analysis has results for """
  run=SMSHelpers.getRun ( analysis, run )
  # we used the exclusion info to get the list
  x=SMSHelpers.motherParticleExclusions ( analysis, run )
  return x.keys()

def getRun ( analysis, run=None ):
  """ tell us, which run the results will be fetched for.
      None if the request cannot be met """
  run=SMSHelpers.getRun ( analysis, run )
  return run

def getAnalyses ( topo, run=None ):
  import os
  """ return all analyses that have results for topo """
  runs=SMSHelpers.runs
  if run: runs= [ run ]
  analyses={}
  for r in runs:
    ## so thats the runs I really have to think about
    dirs=os.listdir ( "%s/%s/" % ( SMSHelpers.Base, r ) )
    for ana in dirs:
      if os.path.exists ( "%s/%s/%s/info.txt" % ( SMSHelpers.Base, r, ana ) ):
        e=getExclusion ( ana, topo, r )
        if e: analyses[ana]=True

  return analyses.keys()

def getAllResults ( run=None ):
  import os
  """ returns all analyses and the topologies they have results for """
  runs=SMSHelpers.runs
  if run: runs= [ run ]
  ret={}
  for r in runs:
    ## so thats the runs I really have to think about
    dirs=os.listdir ( "%s/%s/" % ( SMSHelpers.Base, r ) )
    for ana in dirs:
      if os.path.exists ( "%s/%s/%s/info.txt" % ( SMSHelpers.Base, r, ana ) ):
        topos=getTopologies ( ana, run )
        ret[ana]=topos
  return ret

def getClosestValue ( Dict, mx, my ):
  """ assuming that Dict is a dictionary of mx,my,ul, get the upper limit
      of the point in Dict that is closest to mx and my. """
  import math
  closest=9999999
  retul=None
  for (dmx,dmv) in Dict.items():
    for (dmy,ul) in dmv.items():
      dist= (mx-dmx)**2 + (my-dmy)**2
      if dist<closest:
        closest=dist
        retul=ul
  return retul

def getInterpolatedUpperLimit ( Dict, mx, my ):
  if Dict==None: return None
  if len(Dict)==0: return None
  """ get the upper limit, interpolate from (at most) four neighbours in Dict """
  #print "[SMSResults.py] debug closest value is",getClosestValue ( Dict,mx,my)
  # Large=99999.
  # Neighbors=[ [Large,0.],[Large,0.],[Large,0.],[Large,0.] ]
  Neighbors=[]
  for (dmx,dmv) in Dict.items():
    for (dmy,ul) in dmv.items():
      dist= (dmx-mx)**2 + (dmy-my)**2 
      if dist==0.: return ul
      Neighbors.append ( [ dist, ul ] )
  Neighbors.sort()
  Neighbors=Neighbors[:4]
  ret=sum( [ y/x for (x,y) in Neighbors ] ) / sum( [ 1./x for (x,y) in Neighbors ] )
  #print "interpolating",Neighbors,ret
  return ret

def getUpperLimitFromDictionary ( analysis, topo, mx=None, my=None, run=None, png=None, interpolate=False ):
  """ shouldnt have to call this directly. It's obtaining an upper limit from the python dictionary """
  if interpolate:
    print "[SMSResults.py] error: need to implement interpolation function for getUpperLimitFromDictionary"
    import sys
    sys.exit(0)
  Dict=SMSHelpers.getUpperLimitDictionary ( analysis, topo, run )
  if Dict==None: return Dict
  ## Dict=addunit ( Dict, "pb" )
  # print "[SMSResults.py] mx=",mx
  if mx==None: return Dict
  ## return getClosestValue ( Dict, mx, my )
  return addunit ( getInterpolatedUpperLimit ( Dict, mx, my ), "pb" )
 
def getSmartUpperLimit ( analysis, topo, masses, massesbranch2=None ):
  """ returns the upper limit for analysis/topo, given an ordered sequence of
      the mass (mother, intermediate, LSP) """
  return getUpperLimit ( analysis, topo, mx=masses[0], my=masses[-1], interpolate=True )

def getUpperLimit ( analysis, topo, mx=None, my=None, run=None, png=None, interpolate=False ):
  """ get the upper limit for run/analysis/topo.
      return none if it doesnt exist.
      if mx and my are none, return the entire histogram,
      if mx and my are floats, return the upper limit at this
      point
      if png==True, return path of pngfile containing the histogram"""
  run=SMSHelpers.getRun ( analysis, run )
  if SMSHelpers.hasDictionary ( analysis, run ):
    return getUpperLimitFromDictionary ( analysis, topo, mx, my, run, interpolate )
  histo=SMSHelpers.getUpperLimitFromHisto ( analysis, topo, run )
  if png==True:
    pngfile=SMSHelpers.getUpperLimitPng(analysis,topo,run)
    return pngfile
  if rmvunit(mx,'GeV') == None:
    return histo
  value=SMSHelpers.getUpperLimitAtPoint ( histo, mx, my, interpolate=interpolate )
  if value==0.0: value=None # 0.0 usually means out of bounds
  return SMSHelpers.addunit ( value, "pb" )

def getEfficiency ( analysis, topo, mx=None, my=None, run=None ):
  """ get the efficiency for run/analysis/topo.
      return none if it doesnt exist.
      if mx and my are none, return the entire histogram,
      if mx and my are floats, return the upper limit at this
      point """
  run=SMSHelpers.getRun ( analysis, run )
  histo=SMSHelpers.getEfficiencyHisto ( analysis, topo, run )
  if mx==None:
    return histo
  value=SMSHelpers.getEfficiencyAtPoint ( histo, mx, my )
  return value

def getExplanationForLackOfUpperLimit ( analysis, topo, mx=None, my=None, \
                                        run=None, number=False ):
  """ if there's no upper limit, we want to know what's wrong.
      If number is false, return a text, if number is true
      return the error code """
  value=getUpperLimit ( analysis, topo, run=run )
  msg=SMSHelpers.getErrorMessage ( value, mx, my )
  if number: return msg[0]
  return msg[1]

def isPublic ( analysis, run=None ):
  """ is the result from this analysis public? """
  value=SMSHelpers.getMetaInfoField ( analysis, "private", run )
  if not value:
    return None
  try:
    d=not bool ( int ( value ) )
    return d
  except Exception:
    print "[SMSResults.py] couldnt parse ``private'' field"
  return None

def hasDataPublished ( analysis, run=None ):
  """ has the analysis published their data in digitized form? """
  value=SMSHelpers.getMetaInfoField ( analysis, "publisheddata", run )
  if type(value)==type("string"):
    if value=="None": return None
    if value=="True" or value=="1": return True
    if value=="False" or value=="0": return False
  if not value:
    return None
  try:
    return bool(value)
  except Exception,e:
    print "[SMSResults.py] couldnt parse ``publisheddata'' field"
  return None

def getLumi ( analysis, run=None ):
  """ get the integrated luminosity for this analysis """
  lumifb=float(SMSHelpers.getMetaInfoField ( analysis, "lumi", run ))
  return SMSHelpers.addunit ( lumifb, "fb-1" )

def getSqrts ( analysis, run=None ):
  """ get s_hat for this analysis """
  sqrts=SMSHelpers.getMetaInfoField ( analysis, "sqrts", run )
  try:
    return SMSHelpers.addunit ( float(sqrts), "TeV" )
  except:
    pass
  return sqrts

def getPAS ( analysis, run=None ):
  """ get the PAS for this analysis """
  return SMSHelpers.getMetaInfoField ( analysis, "pas", run )

def hasDictionary ( analysis, run=None ):
  """ are the upper limits available in dictionary format? """
  return SMSHelpers.hasDictionary ( analysis, run )

def getx ( analysis, topo=None, run=None ):
  """ get the description of the x-values for this analysis, if you supply a
      topo, then the return value is the x-values only for this topo """
    
  st = SMSHelpers.getMetaInfoField ( analysis, "x", run )
  if not st:
     return None
  st = st.split(',')
  d = {}
  for i in range(len(st)):
     l = st[i].split(':')
     x = l[1].split()
     d[l[0].replace(" ","")] = x
     
  if topo:
    topo = topo.replace(" ","")
    if not d or not d.has_key ( topo ): return None
    else: return d[topo]
       
  return d

def getaxes (analysis, topo=None, run=None):
  """ get axes for this analysis: for each topo list of dictionary, each dictionary corresponds
	to one histogram, the key axes gives string (mx-my), the key mz gives information on other masses,
	if you supply a topo, returns list for this topo only"""
  st = SMSHelpers.getMetaInfoField (analysis, "axes", run)
  if not st: return None
  st = st.split(',')
  d = {}
  for i in range(len(st)):
    l = st[i].split(':')
    nm = l[0].replace(" ","")
    d[nm] = []
    m = l[1].split('-')
    for j in range(len(m)):
      n = m[j].split()
      if len(n)==2:
	d[nm].append({'axes': n[0]+'-'+n[1], 'mz': None})
      else: d[nm].append({'axes': n.pop(0)+'-'+n.pop(0), 'mz': n})
      
  if topo:
    topo = topo.replace(" ","")
    if not d or not d.has_key ( topo ): return None
    else: return d[topo]

  return d 
	
      
	
def getFigures ( analysis, run=None ):
  """ get the figure number for this analysis """
  return SMSHelpers.getMetaInfoField ( analysis, "figures", run )

def getComment ( analysis, run=None ):
  """ an option comment? """
  return SMSHelpers.getMetaInfoField ( analysis, "comment", run )

def getConditions ( analysis, topo="all", run=None ):
  """ get the conditions. if topo is "all",
      returns a dictionary, else it returns the condition
      only for the given topo, None if non-existent. """
  run=SMSHelpers.getRun ( analysis, run )
  ret = SMSHelpers.conditions ( analysis, run )
  if topo=="all": return ret
  if not ret.has_key ( topo ): return None
  return ret[topo]

def getConstraints ( analysis, topo="all", run=None ):
  """ get the constraints. if topo is "all", 
      returns a dictionary, else it returns the constraint
      only for the given topo, None if non-existent. """
  run=SMSHelpers.getRun ( analysis, run )
  ret = SMSHelpers.constraints ( analysis, run )
  if topo=="all": return ret
  if not ret.has_key ( topo ): return None
  return ret[topo]

def getRequirement ( analysis, run=None ):
  """ any requirements that come with this analysis? (e.g. onshellness) """
  return SMSHelpers.getMetaInfoField ( analysis, "requires", run )

def getCheckedBy ( analysis, run=None ):
  """ has the result been validated? """
  return SMSHelpers.getMetaInfoField ( analysis, "checked", run )

def getJournal ( analysis, run=None ):
  """ get the journal of this analysis """
  return SMSHelpers.getMetaInfoField ( analysis, "journal", run )

def getBibtex ( analysis, run=None ):
  """ get the inspire page with the bibtex entry for this analysis """
  return SMSHelpers.getMetaInfoField ( analysis, "bibtex", run )

def getURL ( analysis, run=None ):
  """ get the URL for this analysis """
  return SMSHelpers.getMetaInfoField ( analysis, "url", run )

def hasURL ( analysis, run=None ):
  """ see if an URL is known """
  return SMSHelpers.hasMetaInfoField ( analysis, "url", run )

def getContact ( analysis, run=None ):
  """ get the contact for this analysis """
  return SMSHelpers.getMetaInfoField ( analysis, "contact", run )

def getPerturbationOrder ( analysis, run=None ):
  """ get the contact for this analysis """
  return SMSHelpers.getMetaInfoField ( analysis, "order", run )

def particles(topo,plot = 'ROOT'):
  """return the production mode for a given topology:
     latex code either compatible to ROOT or Python"""
  if SMSHelpers.dicparticle.has_key(topo):
    part = SMSHelpers.dicparticle[topo]
    if plot == 'ROOT':
      #print "[SMSResults.py] debug",part
      return part
    if plot == 'python':
      return part.replace('#','\\')
  else:
    return None

def particleName(topo):
  """return the production mode for a given topology:
     write out the name in plain letters, no latex """
  if topo[:2]=="TGQ": return "associate"
  if topo=="TChiSlep" or topo=="TChiNuSlep": return "chargino neutralino"
  if topo=="TChiSlepSlep": return "chargino neutralino"
  if topo=="TChiwz": return "chargino neutralino"
  if not SMSHelpers.dicparticle.has_key(topo):
    return "???"
  part = SMSHelpers.dicparticle[topo].replace("#tilde","").replace("{","").replace("}","")
  if part=="g": part="gluino"
  if part=="b": part="sbottom"
  if part=="t": part="stop"
  if part=="q": part="squark"
  return part

def massDecoupling_ ( topo ):
  if topo=="T2tt":
    return "m(#tilde{g},#tilde{q})>>m(#tilde{t})"
  if topo=="T2FVttcc":
    return "m(#tilde{g},#tilde{q})>>m(#tilde{t})"
  if topo=="T2bb":
    return "m(#tilde{g},#tilde{q})>>m(#tilde{b})"
  if topo=="TChiSlep":
    #return "m(#tilde{g}),m(#tilde{q})>>m(#tilde{#chi}^{0}_{2}),m(#tilde{#chi}^{#pm})"
    return "m(#tilde{g},m(#tilde{q})>>m(#tilde{#chi}^{0}_{2},#tilde{#chi}^{#pm})"
  if topo=="TChiSlepSlep":
#return "m(#tilde{g}),m(#tilde{q})>>m(#tilde{#chi}^{0}_{2})"
    return "m(#tilde{g},#tilde{q})>>m(#tilde{#chi}^{0}_{2},#tilde{#chi}^{#pm})"
  if topo=="TChiNuSlep":
    return "m(#tilde{g},#tilde{q})>>m(#tilde{#chi}^{0}_{2},#tilde{#chi}^{#pm})"
    #return "m(#tilde{g}),m(#tilde{q})>>m(#tilde{#chi}^{0}_{2}),m(#tilde{#chi}^{#pm})"
  if topo=="TChiwz":
    return "m(#tilde{g},#tilde{q})>>m(#tilde{#chi}^{0}_{2},#tilde{#chi}^{#pm})"
  T2=topo[:2]
  if T2=="T1" or T2=="T3" or T2=="T5":
    return "m(#tilde{q})>>m(#tilde{g})"
  if T2=="T2" or T2=="T4" or T2=="T6":
    return "m(#tilde{g})>>m(#tilde{q})"
  return ""

def massDecoupling ( topo, plot='ROOT',kerning=True ):
  """ describe the assumed mass decoupling """
  md=massDecoupling_ ( topo )
  if kerning:
    md=md.replace(">>",">#kern[-.2]{>}")
  if plot!='ROOT':
    md = '$' + md.replace('#','\\') + '$'
  return md



def exists(analysis, topo, run = None):
  """ check if the histogram run/analysis/sms.root(limit_topo) exists."""
  """ or sms.py, if it's stored in dictionaries """
  run2=SMSHelpers.getRun( analysis, run )
  hasDict=SMSHelpers.hasDictionary ( analysis, run2 )
  if hasDict:
    histo=SMSHelpers.getUpperLimitDictionary ( analysis, topo, run2 )
    if not histo or len(histo)==0: return False
    return True
  histo=SMSHelpers.getUpperLimitFromHisto(analysis, topo, run2 )
  if not histo: return False
      
  return True
