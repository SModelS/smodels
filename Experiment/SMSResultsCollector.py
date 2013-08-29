#!/usr/bin/env python

"""
.. module:: SMSResultsCollector
    :synopsis: Getting information from the results database, should not be used directly. Use SMSResults.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>, Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>, Doris Proschofsky <Doris.Proschofsky@assoc.oeaw.ac.at>

"""


import os, sys, types

dicpath, lumis, dicexclusions = {},{},{}

Base = "/afs/hephy.at/user/w/walten/public/sms/"

verbose=True

alldirectories = ['8TeV', 'ATLAS8TeV', '2012', '2011']


def SMSInfo(obj,topo = None, ana = None, xvalue = '',plot = 'ROOT', kerning = True,year=None, name=None ):
  """ calls several functions to retrieve informations and objects from "/afs/hephy.at/user/w/walten/public/sms/"
      (acts as a distributor for inquiry)
      possible objects: alltopos, analyses, limit, efficiency, exclusionline, expectedexclusion,
      lumi, pas, contact, analysisname, particle, decay, exclusions, isPublic"""
  if topo != None:
    if topo[0:4] != 'T5zz':
      if topo[-3:] == '025' or topo[-3:] == '075' or topo[-3:] == '050':
        xvalue = topo[-3:]
        topo = topo[:-3]
    if topo[0:4] == 'T5zz':
      if topo[-1:] == 'l' or topo[-3:] == '075' or xvalue == '075':
        xvalue = 'l'
        topo = 'T5zz'
      if topo[-1:] == 'h' or topo[-3:] == '025'or xvalue == '025':
        xvalue = 'h'
        topo = 'T5zz'
      if topo[-3:] == '050' or xvalue == '050':
        xvalue = topo[-3:]
        topo = 'T5zz'

  if obj == 'alltopos':
    return topos(ana)

  if obj == 'analyses':
    return analyses(topo, year=year)

  if obj == 'limit' or obj == 'efficiency' or obj[:13] == 'exclusionline' or \
      obj[:17]== "expectedexclusion" or obj.find("exclusion")>-1:
    if ana == None or topo == None and verbose:
      print 'Object ',obj,'requires analysis and topology.'
      print 'given: analysis is ',ana,'topology is ', topo
    return SMSObjects(obj,topo,ana,xvalue = xvalue)

  if obj == 'lumi' or obj == 'pas' or obj == 'contact' or obj=='order' or obj=='perturbationorder':
    return AnalysisInfo(obj,ana)

  if obj == 'analysisname':
    return analysisname(ana,plot = plot)

  if obj == 'particle':
    return particles(topo,plot = plot)

  if obj == 'decay':
    return decays(topo,plot = plot, kerning = kerning)

  if obj == 'exclusions':
    return exclusions(topo,ana,xvalue=xvalue, plusminussigma=0)

  if obj == 'isPublic':
    return isPublic(topo, ana,xvalue=xvalue)

  else:
    print '[SMSInfo]', obj,' not supported'

def analyses(topo=None, year=None):
  """ retrieve all available analyses from "/afs/hephy.at/user/w/walten/public/sms/"
      using info.txt when called with topo, returns NONE for invalid topo """

  allyears = alldirectories
  if year != None:
    allyears = [year]

  allanalyses = os.listdir(Base+allyears[0])

  if topo != None:
    if not topo in SMSInfo('alltopos'):
      if verbose:
        print '[SMSResultsCollector] invalid topology: ', topo
      return None
    anal  = []
    for i in allanalyses:
      if os.path.exists(Base+allyears[0]+ '/'+ i + '/info.txt'):
        f = open(Base+allyears[0]+ '/' + i + '/info.txt')
        content = f.readlines()
        f.close()
        for j in content:
          store = j.split(' ')
          if store[0] == 'exclusions:' and store[1] == topo:
            anal.append(i)
            break
      allanalyses = anal

  for year in allyears:
    nextanalyses = os.listdir(Base+year)
    for i in nextanalyses:
      compare = True

      if topo != None and os.path.exists(Base+year+ '/'+ i + '/info.txt'):
        f = open(Base+year+ '/' + i + '/info.txt')
        content = f.readlines()
        f.close()
        for j in content:
          store = j.split(' ')
          if store[0] == 'exclusions:' and store[1] == topo:
            compare = False
            break

      if i in allanalyses:
        compare = True
      if compare == False:
        allanalyses.append(i)

  anal = []
  for ana in allanalyses:
    if not '.' in ana and ana != 'TODO' and ana != 'old':
      anal.append(ana)
  allanalyses = anal

  if len(allyears) == 1:
    allanalyses = [allyears[0]+ana for ana in allanalyses]


  return allanalyses


def topos(ana = None,Masssplitting = False):
  """retrieve all available topologies from "/afs/hephy.at/user/w/walten/public/sms/"
     using info.txt , returns NONE for invalid analysis"""

  allyears = alldirectories
  if ana != None:
    if ana[0:4] in alldirectories:
      allyears = [ana[0:4]]
      allanalyses = [ana[4:]]
      if not ana in SMSInfo('analyses',year=ana[0:4]):
        if verbose:
          print '[SMSResultsCollector] invalid analysis: ', ana
        return None
    else:
      allanalyses = [ana]
      if not ana in SMSInfo('analyses'):
        if verbose:
          print '[SMSResultsCollector] invalid analysis: ', ana
        return None
  alltopos = []

  for year in allyears:
    if ana == None:
      allanalyses = analyses(year=year)

    for anas in allanalyses:
      if ana == None:
        anas = anas[4:]

      if not os.path.exists(Base+year+ '/'+ anas + '/info.txt'):
        continue
      f = open(Base+year+ '/' + anas + '/info.txt')
      content = f.readlines()
      f.close()
      for j in content:
        store = j.split(' ')
        if store[0] == 'exclusions:':
	  if not Masssplitting:
            if store[1][-3:] == '025' or store[1][-3:] == '050' or store[1][-3:] == '075':
              store[1] = store[1][:-3]
            if store[1] == 'T5zzl' or store[1] =='T5zzh':
              store[1] = store[1][:-1]
          if not store[1] in alltopos:
            alltopos.append(store[1])
  alltopos.sort()
  return alltopos



def SMSObjects(obj,topo,ana,xvalue='',name=None):
  """ return exclusionlines, expectedexclusionlines, limit-histogramms and efficiency-histogramms"""

  Object={}
  if type(xvalue)==types.FloatType:
    xvalue=str(xvalue)
  if xvalue==None:
    xvalue=""
  if xvalue==".5" or xvalue=="0.5":
    xvalue='050'

  if topo != 'T5zz':
#    if xvalue == '':
#      x = ['025','','075']
#    elif xvalue == '050':
#      x = ['']
    x = [xvalue]
  if topo == 'T5zz':
    if xvalue == '':
      x = ['l','','h']
    elif xvalue == '050':
      x = ['']
    else:
      x = [xvalue]

  if obj == 'exclusionline' or obj == 'exclusionline3' or obj =='exclusionline13' or \
     obj.find('expectedexclusionline')==0:
    obj = obj.replace('line','')
  obj = obj + '_' + topo

  allyears = alldirectories

  if ana[0:4] in alldirectories:
    allyears = [ana[0:4]]
    ana = ana[4:]

  for year in allyears:
    filename=Base + year + "/" + ana + '/sms.root'
    if not os.path.exists(filename):
      continue
    from ROOT import TFile
    f = TFile(filename)
    if not f or not f.IsOpen():
      if verbose:
        print "[SMSResultsCollector.py] error: cannot open",filename
      continue
    for j in x:
      if j == "":
        key = '050'
      elif j == 'l':
        key = '075'
      elif j == 'h':
        key = '025'
      else:
        key = j
      Ob = f.Get(obj+j)
      if Ob:
        if obj[:9] != 'exclusion' and obj[:17] != 'expectedexclusion':
          Ob.SetDirectory(0)
        Object.update({key:Ob})
        if verbose:
          print '[SMSInfo]', obj+j,' for ', ana, ' found in ', year,'file',filename,"obj=",obj

    if Object !={}:
      return [year,Object]
    f.Close()

  if verbose:
    print '[SMSInfo] no %s %s for %s and %s in %s-%s\n' % ( obj,x,ana,topo,allyears[0],allyears[-1] )
  return None


def AnalysisInfo(obj,ana):
  """ retrieve the integrated luminosity or PAS-number
      or contact-information for a specific analysis """
  contend=[]
  allyears = alldirectories

  if ana[0:4] in alldirectories:
    allyears = [ana[0:4]]
    ana = ana[4:]

  for year in allyears:
      filename=Base + year + "/" + ana + '/info.txt'
      if not os.path.exists(filename):
        continue
      f = open(filename)
      content = f.readlines()
      content = [i.replace('\n','') for i in content]
      f.close()
      for i in content:
        store = i.split(':')

        if store[0] == obj:
          if verbose:
            print '[SMSInfo]', obj, ' found in', year,store
          if obj=="private" and len(store)>0:
            b=bool(int(store[1]))
            return [year,[b]]
          return [year,store[1][1:]]

  if obj=='order': return "nlo"
  if verbose:
    ##print '[SMSInfo]', 'no ',obj, ' for',ana,'in', allyears[0],'-',allyears[-1],"\n"
    print '[SMSInfo] no %s for %s in %s-%s.\n' % ( obj, ana, allyears[0], allyears[-1] )
  return None

def pas(ana):
  """ retrieve the PAS number for the analysis"""
  return AnalysisInfo('pas',ana)

def lumi(ana):
  """ retrieve the lumi for the analysis"""
  return AnalysisInfo('lumi',ana)

def order(ana):
  """ return what perturbation order the exclusion contours are.
      Values are nlo or nll """
  return AnalysisInfo('order',ana)

def perturbationOrder(ana):
  return order(ana)

def contact(ana):
  """ retrieve the contact person for the analysis"""
  return AnalysisInfo('contact',ana)


def isPublic(topo,ana,xvalue=''):
  """retrieve publicinfo for a Plot
     return None if plot does not exist
     return the latest public Plot if no "year" is given"""

  if topo != 'T5zz':
    if xvalue == '':
      x = ['025','','075']
    elif xvalue == '050':
      x = ['']
    else:
      x = [xvalue]
  if topo == 'T5zz':
    if xvalue == '':
      x = ['l','','h']
    elif xvalue == '050':
      x = ['']
    else:
      x = [xvalue]

  store ={}
  alltopos = topos(ana = ana,Masssplitting = True)
  for j in x:
    if  topo+j in alltopos:
      store[j]=True
  if store == {}:
    return None

  allyears = alldirectories
  if ana[0:4] in alldirectories:
    allyears = [ana[0:4]]
    ana = ana[4:]

  for year in allyears:
    PublicInfo = AnalysisInfo('private',year+ana)
    print "publicinfo=",PublicInfo
    hit = False
    if PublicInfo == None:
      break
    print "store keys=",store.keys()
    for k in store.keys():
      if topo+k in PublicInfo[1]:
        store[k] = PublicInfo[1][k]
        hit = True
    if hit == False:
      break

  Public ={}
  for (j,value) in store.items():
    if j == "":
      key = '050'
    elif j == 'l':
      key = '075'
    elif j == 'h':
      key = '025'
    else:
      key = j
    Public[key]=value

  if hit == True and year == allyears[-1] and len(allyears) != 1:
    year = ''

  return [year,Public]


def analysisname(ana, plot):
  """ return pretty name of a given analysis: latex code either compatible to ROOT or Python """

  if dicanalysesnames.has_key(ana):
    name = dicanalysesnames[ana]
    if plot == 'ROOT':
      return name.replace('$',"")
    if plot == 'python':
      return name.replace('#','\\') ## .replace('geq','lt')
  else:
    return ana

dicanalysesnames= {
   "RA4LS": "$e/#mu$ LS",
   "SS": "SS $e/#mu$",
   "LeptonsMT": "3l + $M_{T}$",
   "SSnoHT": "SS no $H_{T}$",
   "RAZOR": "razor",
   "RAZORb": "razor+b",
   "RAZORjets": "razor multijets",
   "MLSHAPE": "2l 2j",
   "LeptonsMET": "multilepton $(#geq 3) + #slash{E}_{T}$",
   "MULTILEPTON": "multilepton $(#geq 3)$",
   "COMBLEPTONS": "comb. leptons",
   "OSedge": "OS $e/#mu$ edge",
   "OSANN": "OS $e/#mu$ ANN",
   "OScounting": "OS $e/#mu + #slash{E}_{T}$",
   "OSshape": "OS shape",
   "RA4b": "$e/#mu #geq 3 b, Y_{MET}$",
   #"RA4b": "$e/#mu+b$",
   "RA4LP": "$e/#mu$ LP",
   "RA4MET": "$e/#mu #geq 2 b + #slash{E}_{T}$",
   #"RA4MET": "$e/#mu$ templates",
   "RA4ANN": "$e/#mu$ ANN",
   "PhotonsMET": "$#gamma#gamma j+ #slash{E}_{T}$",
   "PhotonJETS": "$#gamma jj +#slash{E}_{T}$",
#   "PhotonsMET": "$#gamma#gamma+#slash{E}_{T}$",
#   "PhotonJETS": "$#gamma$+jets",
   "SSb": "SS $e/#mu$ + b",
   "RA2": "$#slash{H}_{T}$ + jets",
   "RA2b": "$#slash{E}_{T}$+b",
   "alphaT": "$#alpha_{T}$",
   "alphaTb": "$#alpha_{T}$+b",
   "ZMET":"Z + $#slash{E}_{T}$",
   "JZB":"JZB",
   "VZMET": "?",
   "VZMT": "?",
   "SSEWK": "?",
   "MT2": "$M_{T2}$",
   "MT2b": "$M_{T2}$b",
   "LQ3": "LQ3",
   "HadronicStop": "hadronic stop",
   "LeptonicStop": "leptonic stop (8 TeV)",
   "MultiLepton8TeV": "multileptons (8 TeV)",
   "MonoJet8TeV": "monojet (8 TeV)",
   "RazorMono8TeV": "razor-mono (8 TeV)",
   "DileptonicStop8TeV": "stop + 2 leptons (8 TeV)",
   "SSb8TeV": "SS+b (8 TeV)",
}

def particleName(topo):
  """return the production mode for a given topology:
     write out the name in plain letters, no l
     atex """
  if topo[:2]=="TGQ": return "associate"
  if topo=="TChiSlep" or topo=="TChiNuSlep": return "chargino neutralino"
  if topo=="TChiSlepSlep": return "chargino neutralino"
  if not dicparticle.has_key(topo):
    return "???"
  part = dicparticle[topo].replace("#tilde","").replace("{","").replace("}","")
  if part=="g": part="gluino"
  if part=="b": part="sbottom"
  if part=="t": part="stop"
  if part=="q": part="squark"
  return part

def particles(topo,plot = 'ROOT'):
  """return the production mode for a given topology:
     latex code either compatible to ROOT or Python"""
  if dicparticle.has_key(topo):
    part = dicparticle[topo]
    if plot == 'ROOT':
      return part
    if plot == 'python':
      return part.replace('#','\\')
  else:
    return None

def production(topo,plot='ROOT'):
  a=particles(topo)
  pair=str(a)+" "+str(a)
  if topo=="TGQ":
    pair=particles("T1")+" "+particles("T2")
  if topo=="TChiSlep":
    pair="#tilde{#chi}^{0}_{2}#tilde{#chi}^{#pm}_{1}"
  if topo=="TChiSlep":
    pair="#tilde{#chi}^{0}_{2}#tilde{#chi}^{#pm}_{1}"
  if topo=="TChiwz":
    pair="#tilde{#chi}^{0}_{2}#tilde{#chi}^{#pm}_{1}"
  if topo=="TChizz":
    pair="#tilde{#chi}^{0}_{2}#tilde{#chi}^{0}_{3}"
  if topo=="TChiNuSlep":
    pair="#tilde{#chi}^{0}_{2}#tilde{#chi}^{#pm}_{1}"
  if topo=="TChiSlepSlep":
    pair="#tilde{#chi}^{0}_{2}#tilde{#chi}^{#pm}_{1}"
#pair="#tilde{#chi}^{0}_{2}#tilde{#chi}^{0}_{2}"

  if plot == 'ROOT':
    return pair
  if plot == 'python':
    return pair.replace('#','\\')

dicparticle = {
   "T1": "#tilde{g}",
   "T1bbbb": "#tilde{g}",
   "T1tttt": "#tilde{g}",
   "T1lg": "#tilde{g}",
   "T1gg": "#tilde{g}",
   "T5gg": "#tilde{g}",
   "T5wg": "#tilde{g}",
   "T2": "#tilde{q}",
   "T2bb": "#tilde{b}",
   "T2bw": "#tilde{t}",
   "T2tt": "#tilde{t}",
   "T2FVttcc": "#tilde{t}",
   "T6ttww": "#tilde{b}",
   "T6ttWW": "#tilde{b}",
   "T6bbWW": "#tilde{t}",
   "T6bbZZ": "#tilde{b}",
   "T2ttww": "#tilde{b}",
   "T3w": "#tilde{g}",
   "T3lh": "#tilde{g}",
   "T5zz": "#tilde{g}",
   "T5zzInc": "#tilde{g}",
   #"T5zzh": "#tilde{g}",
   #"T5zzl": "#tilde{g}",
   "T5lnu": "#tilde{g}",
   "T5zzgmsb": "#tilde{g}",
   "TChizz": "#tilde{#chi}^{0}_{1}",
   "TChiSlep": "#tilde{#chi}^{0}_{1}",
   "TChiSlepSlep": "#tilde{#chi}^{0}_{2}#tilde{#chi}^{#pm}_{1}",
   "TChiwz": "#tilde{#chi}^{0}_{2}#tilde{#chi}^{#pm}_{1}",
   "TChiNuSlep": "#tilde{#chi}^{0}_{1}",
}

def description(topo,plot='ROOT',kerning=True,short=False):
  """ give a long description of the topology,
      with production, decay, and mass decoupling """
  if short:
    return "pp #rightarrow %s %s; %s" % ( production(topo,plot),  decays ( topo,plot,kerning, omitleft=True), massDecoupling ( topo,plot,kerning) )
  return "pp#rightarrow %s, %s; %s" % ( production(topo,plot), decays ( topo,plot,kerning, omitleft=False ), massDecoupling ( topo,plot,kerning) )


def decays(topo,plot = 'ROOT', kerning=True, omitleft=False ):
  """ give the pretty decay string for a given topo.
      E.g. T1 -> ~g -> q q ~ch10.
      kerning: means smaller space between > and >
      omitleft means omit everything up to #rightarrow """

  if dicdecay.has_key(topo):
    part = dicdecay[topo]
    if omitleft and part.find("#rightarrow")>-1:
      part=part[part.find("#rightarrow"):]
    if plot=="ROOT":
      return part
    part = '$' + part.replace('#','\\') + '$'
    if kerning:
      part=part.replace(">>",">#kern[-.2]{>}")
    return part
  else:
    return None

def massDecoupling_ ( topo ):
  if topo=="T2tt":
    return "m(#tilde{g},#tilde{q}) >> m(#tilde{t})"
  if topo=="T2bw":
    return "m(#tilde{g},#tilde{q}) >> m(#tilde{t})"
  if topo=="T2FVttcc":
    return "m(#tilde{g},#tilde{q}) >> m(#tilde{t})"
  if topo=="T6ttww":
    return "m(#tilde{g},#tilde{q}) >> m(#tilde{b})"
  if topo=="T2ttww":
    return "m(#tilde{g},#tilde{q}) >> m(#tilde{b})"
  if topo=="T2bb":
    return "m(#tilde{g},#tilde{q}) >> m(#tilde{b})"
  if topo=="TChiSlep":
    return "m(#tilde{g}),m(#tilde{q}) >> m(#tilde{#chi}^{0}_{2}),m(#tilde{#chi}^{#pm}_{1})"
  if topo=="TChiSlepSlep":
#return "m(#tilde{g}),m(#tilde{q})>>m(#tilde{#chi}^{0}_{2})"
    return "m(#tilde{g}),m(#tilde{q}) >> m(#tilde{#chi}^{0}_{2}),m(#tilde{#chi}^{#pm}_{1})"
  if topo=="TChiNuSlep":
    return "m(#tilde{g}),m(#tilde{q}) >> m(#tilde{#chi}^{0}_{2}),m(#tilde{#chi}^{#pm}_{1})"
  if topo=="TChizz":
    return "m(#tilde{g}),m(#tilde{q}) >> m(#tilde{#chi}^{0}_{2}),m(#tilde{#chi}^{0}_{3})"
  if topo=="TChiwz":
    return "m(#tilde{g}),m(#tilde{q}) >> m(#tilde{#chi}^{0}_{2}),m(#tilde{#chi}^{#pm}_{1})"
  T2=topo[:2]
  if T2=="T1" or T2=="T3" or T2=="T5":
    return "m(#tilde{q}) >> m(#tilde{g})"
  if T2=="T2" or T2=="T4" or T2=="T6":
    return "m(#tilde{g}) >> m(#tilde{q})"
  return ""

def massDecoupling ( topo, plot='ROOT',kerning=True ):
  md=massDecoupling_ ( topo )
  if kerning:
    md=md.replace(">>",">#kern[-.2]{>}")
  if plot!='ROOT':
    md = '$' + md.replace('#','\\') + '$'
  return md

#dicdecayold = {

#  "": "",
#  "T1": "#tilde{g} #rightarrow qq #tilde{#chi}^{0}",
#  "T1bbbb": "#tilde{g} #rightarrow bb #tilde{#chi}^{0}",
#  "T1tttt": "#tilde{g} #rightarrow tt #tilde{#chi}^{0}_{1}",
#  "T2": "#tilde{q} #rightarrow q #tilde{#chi}^{0}",
#  "T2bb": "#tilde{b} #rightarrow b #tilde{#chi}^{0}",
#  "T2tt": "#tilde{t} #rightarrow t #tilde{#chi}^{0}_{1}",
#  "T6ttww": "unknown",
#  "T3w": "#tilde{g} #rightarrow qq [#tilde{#chi}^{0}|#tilde{#chi}^{#pm}]",
#  "T3lh": "#tilde{g} #rightarrow qq#tilde{#chi}^{0}_{2}|#tilde{#chi}^{0}",
#  "T5zz": "#tilde{g} #rightarrow qq #tilde{#chi}^{0}_{2}",
#  "T5zzlnc": "#tilde{g} #rightarrow qq #tilde{#chi}^{0}_{2}", #?
#  "T5zzh": "#tilde{g} #rightarrow qq #tilde{#chi}^{0}_{2}", #?
#  "T5zzl": "#tilde{g} #rightarrow qq #tilde{#chi}^{0}_{2}", #?
#  "T5lnu": "#tilde{g} #rightarrow qq #tilde{#chi}^{#pm}",
#  "T5zzgmsb": "unknown",
#  "TChiwz": "unknown",
#  "TChizz": "unknown",
#  "TChiSlep": "#tilde{#chi}^{0}_{2}|#tilde{#chi}^{#pm} #rightarrow 3l",
#}

dicdecay = { "T1": "#tilde{g} #rightarrow q#bar{q} #tilde{#chi}^{0}_{1}",
    "T2FVttcc": "#tilde{t} #rightarrow c #tilde{#chi}^{0}_{1}",
    "T2llnunubb": "#tilde{t} #rightarrow l #nu b #tilde{#chi}^{0}_{1}",
    "T1bbbb": "#tilde{g} #rightarrow b#bar{b} #tilde{#chi}^{0}_{1}",
    "T1tttt": "#tilde{g} #rightarrow t#bar{t} #tilde{#chi}^{0}_{1}",
    "T2":"#tilde{q} #rightarrow q #tilde{#chi}^{0}_{1}",
    "T2bb":"#tilde{b} #rightarrow b #tilde{#chi}^{0}_{1}",
    "T2bw":"#tilde{t} #rightarrow b (#tilde{#chi}^{#pm}_{1} #rightarrow W #tilde{#chi}^{0}_{1})",
    "T2bW":"#tilde{t} #rightarrow b (#tilde{#chi}^{#pm}_{1} #rightarrow W #tilde{#chi}^{0}_{1})",
    "T2tt": "#tilde{t} #rightarrow t #tilde{#chi}^{0}_{1}",
    "T2FVttcc": "#tilde{t} #rightarrow c #tilde{#chi}^{0}_{1}",
    "T6ttww": "#tilde{b} #rightarrow tW #tilde{#chi}^{0}_{1}",
    "T2ttww": "#tilde{b} #rightarrow tW #tilde{#chi}^{0}_{1}",
    "TChiwz":"#tilde{#chi}^{#pm} #tilde{#chi}^{0}_{2} #rightarrow W Z #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1}",
    "TChizz":"#tilde{#chi}^{0}_{3} #tilde{#chi}^{0}_{2} #rightarrow Z Z #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1}",
    "T6ttWW": "#tilde{b} #rightarrow tW #tilde{#chi}^{0}_{1}",
    "T2ttWW": "#tilde{b} #rightarrow tW #tilde{#chi}^{0}_{1}",
    "TChiWZ":"#tilde{#chi}^{#pm} #tilde{#chi}^{0}_{2} #rightarrow W Z #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1}",
    "TChiZZ":"#tilde{#chi}^{0}_{3} #tilde{#chi}^{0}_{2} #rightarrow Z Z #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1}",
    "TChiSlep":"#tilde{#chi}^{0}_{2} #tilde{#chi}^{#pm}_{1} #rightarrow l l l #nu #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} ",
    "TChiNuSlep":"#tilde{#chi}^{0}_{2} #tilde{#chi}^{#pm}_{1} #rightarrow l l l #nu #tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} ",
    "TChiSlepSlep":"#tilde{#chi}^{0}_{2} #tilde{#chi}^{#pm}_{1} #rightarrow lll #nu#tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} ",
#"TChiSlepSlep":"#tilde{#chi}^{0}_{2} #tilde{#chi}^{0}_{2} #rightarrow llll#tilde{#chi}^{0}#tilde{#chi}^{0} ",
    "T1gg":"#tilde{g} #rightarrow q#bar{q}(#tilde{#chi}^{0}_{2}#rightarrow #gamma#tilde{#chi}^{0}_{1})",
    "T5gg":"#tilde{g} #rightarrow q#bar{q}(#tilde{#chi}^{0}_{2}#rightarrow #gamma#tilde{#chi}^{0}_{1})",
    ## "T1gg":"#tilde{g} #rightarrow qq#tilde{#chi}^{0}_{2},#tilde{#chi}^{0}_{2} #rightarrow #gamma#tilde{#chi}^{0}",
    "T1lg":"#tilde{g} #rightarrow q#bar{q}(#tilde{#chi}^{0}_{2}#rightarrow#gamma#tilde{#chi}^{0}_{1}|#tilde{#chi}^{#pm}_{1}#rightarrow W#tilde{#chi}^{0}_{1})", ## ,#tilde{#chi}^{0}_{2} #rightarrow #gamma#tilde{#chi}^{0}, #tilde{#chi}^{#pm} #rightarrow W#tilde{#chi}^{0}",
    "T5wg":"#tilde{g} #rightarrow q#bar{q}(#tilde{#chi}^{0}_{2}#rightarrow#gamma#tilde{#chi}^{0}_{1}|#tilde{#chi}^{#pm}_{1}#rightarrow W#tilde{#chi}^{0}_{1})", ## ,#tilde{#chi}^{0}_{2} #rightarrow #gamma#tilde{#chi}^{0}, #tilde{#chi}^{#pm} #rightarrow W#tilde{#chi}^{0}",
    "T5Wg":"#tilde{g} #rightarrow q#bar{q}(#tilde{#chi}^{0}_{2}#rightarrow#gamma#tilde{#chi}^{0}_{1}|#tilde{#chi}^{#pm}_{1}#rightarrow W#tilde{#chi}^{0}_{1})",
    #"T1lg":"#tilde{g} #rightarrow qq(#tilde{#chi}^{0}_{2}|#tilde{#chi}^{#pm}),#tilde{#chi}^{0}_{2} #rightarrow #gamma#tilde{#chi}^{0}, #tilde{#chi}^{#pm} #rightarrow W#tilde{#chi}^{0}",
    #"T3w":"#tilde{#chi}^{#pm} #rightarrow W#tilde{#chi}^{0}",
    "T3w": "#tilde{g} #rightarrow q#bar{q}(#tilde{#chi}^{#pm}_{1}#rightarrow W#tilde{#chi}^{0}_{1} |#tilde{#chi}^{0}_{1})", ## #tilde{#chi}^{#pm} #rightarrow W#tilde{#chi}^{0}",
    #"T3w": "#tilde{g} #rightarrow qq(#tilde{#chi}^{#pm}|#tilde{#chi}^{0}), #tilde{#chi}^{#pm} #rightarrow W#tilde{#chi}^{0}",
    "T3wb":"#tilde{g} #rightarrow b#bar{b}(W)#tilde{#chi}^{0}_{1}",
    #"T3w":"#tilde{g} #rightarrow qq(W)#tilde{#chi}^{0}",
    "T3W": "#tilde{g} #rightarrow q#bar{q}(#tilde{#chi}^{#pm}_{1}#rightarrow W#tilde{#chi}^{0}_{1} |#tilde{#chi}^{0}_{1})",
    "T3Wb":"#tilde{g} #rightarrow b#bar{b}(W)#tilde{#chi}^{0}_{1}",
    "T1bbbb":"#tilde{g} #rightarrow b#bar{b} #tilde{#chi}^{0}_{1}",
   # "T1lnu":"#tilde{#chi}^{#pm} #rightarrow l^{#pm}#nu #tilde{#chi}^{0}",
    "T5lnu":"#tilde{g} #rightarrow q#bar{q}(#tilde{#chi}^{#pm} #rightarrow l^{#pm}#nu #tilde{#chi}^{0}_{1})",
    "T1lnu":"#tilde{g} #rightarrow q#bar{q}(#tilde{#chi}^{#pm} #rightarrow l^{#pm}#nu #tilde{#chi}^{0}_{1})",
    "T1lh":"#tilde{g} #rightarrow q#bar{q} #tilde{#chi}^{0}_{2},#tilde{#chi}^{0}_{2} #rightarrow l^{+}l^{-}#tilde{#chi}^{0}_{1}",
    "T3lh":"#tilde{g} #rightarrow q#bar{q} (#tilde{#chi}^{0}_{2}#rightarrow l^{+}l^{-}#tilde{#chi}^{0}_{1})",
    #"T3lh":"#tilde{g} #rightarrow q q #tilde{#chi}^{0}_{2},#tilde{#chi}^{0}_{2} #rightarrow l^{+}l^{-}#tilde{#chi}^{0}",
    "T5zz":"#tilde{g} #rightarrow q#bar{q} (#tilde{#chi}^{0}_{2}#rightarrow Z #tilde{#chi}^{0}_{1})",
    #"T5zz":"#tilde{g} #rightarrow qq #tilde{#chi}^{0}_{2}, #tilde{#chi}^{0}_{2} #rightarrow Z #tilde{#chi}^{0}",
    "T5zzInc":"#tilde{#chi}^{0}_{2} #rightarrow Z #tilde{#chi}^{0}_{1}",
    "T5ZZ":"#tilde{g} #rightarrow q#bar{q} (#tilde{#chi}^{0}_{2}#rightarrow Z #tilde{#chi}^{0}_{1})",
    "T5ZZInc":"#tilde{#chi}^{0}_{2} #rightarrow Z #tilde{#chi}^{0}_{1}",
    #"T5zzl":"#tilde{#chi}^{0}_{2} #rightarrow Z #tilde{#chi}^{0}",
    #"T5zzh":"#tilde{#chi}^{0}_{2} #rightarrow Z #tilde{#chi}^{0}",
    "T2tt":"#tilde{t} #rightarrow t #tilde{#chi}^{0}_{1}",
    "T2FVttcc":"#tilde{t} #rightarrow c #tilde{#chi}^{0}_{1}",
    "T6ttww":"#tilde{b} #rightarrow tW #tilde{#chi}^{0}_{1}",
    "T2ttww":"#tilde{b} #rightarrow tW #tilde{#chi}^{0}_{1}",
    "T6ttWW":"#tilde{b} #rightarrow tW #tilde{#chi}^{0}_{1}",
    "T2ttWW":"#tilde{b} #rightarrow tW #tilde{#chi}^{0}_{1}",
    "T5zzgmsb":"#tilde{g} #rightarrow q#bar{q} (#tilde{#chi}^{0}_{2}#rightarrow Z #tilde{#chi}^{0}_{1})",
    "TChiSlepSlep25":"#tilde{#chi}^{0}_{2} #tilde{#chi}^{#pm} #rightarrow lll #nu#tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} ",
    "TChiSlepSlep75":"#tilde{#chi}^{0}_{2} #tilde{#chi}^{#pm} #rightarrow lll #nu#tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} ",
    "T6bbWW":"#tilde{t} #rightarrow b(#tilde{#chi}^{+} #rightarrow W#tilde{#chi}^{0}_{1})",
    "TChiChipmSlepL":"#tilde{#chi}^{0}_{2} #tilde{#chi}^{#pm} #rightarrow lll #nu#tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} ",
    "TChiChipmSlepSlep":"#tilde{#chi}^{0}_{2} #tilde{#chi}^{#pm} #rightarrow lll #nu#tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} ",
    "TChiChipmSlepStau":"#tilde{#chi}^{0}_{2} #tilde{#chi}^{#pm} #rightarrow ll#tau #nu#tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} ",
    "TChiChipmStauStau":"#tilde{#chi}^{0}_{2} #tilde{#chi}^{#pm} #rightarrow #tau#tau#tau #nu#tilde{#chi}^{0}_{1} #tilde{#chi}^{0}_{1} ",
    "TChipChimSlepSnu":"#tilde{#chi}^{+}#tilde{#chi}^{-} #rightarrow l^{+}l^{-}#nu#nu#tilde{#chi}^{0}_{1}#tilde{#chi}^{0}_{1}",
    "TSlepSlep":"#tilde{l} #rightarrow l #tilde{#chi}^{0}_{1}",
    "Hadronic112q":"#tilde{q}_{R} #rightarrow qqqq  #lambda''_{112}",
    "Leptonic233q":"#tilde{q} #rightarrow qll#nu  #lambda_{233}",
    "Leptonic122g":"#tilde{g} #rightarrow qll#nu  #lambda_{122}",
    "SemiLeptonic233g":"#tilde{g} #rightarrow qbt#mu  #lambda'_{233}",
    "Leptonic233g":"#tilde{g} #rightarrow qll#nu  #lambda_{233}",
    "Hadronic122q":"#tilde{q}_{R} #rightarrow qqqq  #lambda_{122}",
    "SemiLeptonic233q":"#tilde{q} #rightarrow qbt#mu  #lambda'_{233}",
    "bprime":"b' #rightarrow bZ",
    "Stop233":"#tilde{t}_{R} #rightarrow #mu#tau#nut  #lambda_{233}",
    "Leptonic122q":"#tilde{q} #rightarrow qll#nu  #lambda_{122}",
    "Hadronic112g":"#tilde{g} #rightarrow qqqq  #lambda''_{112}",
    "Hadronic122g":"#tilde{g} #rightarrow qqqq  #lambda_{122}",
    "Stop122":"#tilde{t}_{R} #rightarrow #mue#nut  #lambda_{122}",
    "Stop123":"#tilde{t}_{R} #rightarrow #mu#tau#nut  #lambda_{123}",
    "Leptonic123g":"#tilde{g} #rightarrow qll#nu  #lambda_{123}",
    "Leptonic123q":"#tilde{q} #rightarrow qll#nu  #lambda_{123}",
    "StopLLE122":"#tilde{t}_{R} #rightarrow t#nu_{#mu}e#mu  #lambda_{122}",
    "StopLLE233":"#tilde{t}_{R} #rightarrow t#nu_{#tau}#mu#tau  #lambda_{233}",
    "StopLQD233":"#tilde{t}_{R} #rightarrow tbt#mu  #lambda'_{233}",
    "T7btW":"#tilde{g} #rightarrow btW#tilde{#chi}^{0}_{1}",
    "T5tttt":"#tilde{g} #rightarrow t(#tilde{t} #rightarrow t#tilde{#chi}^{0}_{1})",
    "T3tauh":"#tilde{g} #rightarrow q#bar{q}(#tilde{#chi}^{2}#rightarrow #tau #tau #tilde{#chi}^{0}_{1} |#tilde{#chi}^{0}_{1})",
    "T6bbZZ":"#tilde{b} #rightarrow bZ #tilde{#chi}^{0}_{1}",
    "Rstop":"#tilde{t} #rightarrow b#tau  #lambda'_{333}",
    "Rg3j":"#tilde{g} #rightarrow qqq  #lambda''_{112}",
    "SemiLeptonic231g":"#tilde{g} #rightarrow qbt#mu  #lambda'_{231}",
    "SemiLeptonic231q":"#tilde{q} #rightarrow qbt#mu  #lambda'_{231}",
    "T5WW":"#tilde{g} #rightarrow q#bar{q}(#tilde{#chi}^{#pm}_{1}#rightarrow W#tilde{#chi}^{0}_{1})",
    "T7btbtWW":"#tilde{g} #rightarrow b(#tilde{b} #rightarrow t(#tilde{#chi}^{#pm} #rightarrow W#tilde{#chi}^{0}_{1}))",
    "Gluino113/223":"#tilde{g} #rightarrow qqb  #lambda''_{113/223}",
    "Gluino323":"#tilde{g} #rightarrow tbs  #lambda''_{323}"
}


def exclusions(topo,ana,xvalue='',expected = False, plusminussigma=0):
  """ retrieve the minimum and maximum exclusion values, as
    they are needed for the summary plot """

  contend=[]
  glmin = {}
  glmax = {}
  minx = {}
  allyears = alldirectories

  if topo != 'T5zz':
    if xvalue == '':
      x = ['025','','075']
    elif xvalue == '050':
      x = ['']
    else:
      x = [xvalue]
  if topo == 'T5zz':
    if xvalue == '':
      x = ['l','','h']
    elif xvalue == '050':
      x = ['']
    else:
      x = [xvalue]

  if ana[0:4] in allyears:
    allyears = [ana[0:4]]
    ana = ana[4:]
  excltype = 'exclusions'
  if plusminussigma == 1: excltype = excltype+'p1'
  if plusminussigma == -1: excltype = excltype+'m1'
  if expected == True: excltype = 'expected' + excltype

  for year in allyears:
      if os.path.exists(Base + year + "/" + ana + '/info.txt'):
        f = open(Base + year + "/" + ana + '/info.txt')
        content = f.readlines()
        content = [i.replace('\n','') for i in content]
        f.close()
        for j in x:
          for i in content:
            store = i.split(" ")
            if store[0] == excltype + ':' and store[1] == topo + j:
              if j == '':
                key = '050'
              elif j == 'l':
                key = '075'
              elif j == 'h':
                key = '025'
              else:
                key = j
              if len(store) == 5:
                minx.update({key:float(store[2])})
                glmin.update({key:float(store[3])})
                glmax.update({key:float(store[4])})
              if len(store) == 4:
                minx.update({key:0.})
                glmin.update({key:float(store[2])})
                glmax.update({key:float(store[3])})

              if verbose:
                print '[SMSInfo]', excltype, key, ' found in', year

        if glmin != {} and glmax != {} and minx !={}:
          return [year,minx,glmin,glmax]

  if verbose:
    print '[SMSInfo]', 'no',excltype, 'for',ana," and ", topo ##, 'in', startyear,'-',endyear
  return None

def ul(topo,ana,xvalue=''):
  return exclusions(topo,ana,xvalue)

def upperlimit(topo,ana,xvalue=''):
  return exclusions(topo,ana,xvalue)

def exclusionline(topo,ana,xvalue='',factor=1.0,extendedinfo=True,expected=False,plusminussigma=0):
  """ get the exclusion line (TGraph)
     extendedinfo = True returns a list [run, {xvalue: exclusionline}] (if there is no xvalue e.g. T1 the key is 050),
     extendedinfo = False returns just the exclusionline """
  Type="exclusion"
  if abs(factor-3.0)<0.1:
    Type="exclusion3"
  if abs(factor-1./3.)<0.1:
    Type="exclusion13"
  if expected:
    Type="expected"+Type
  if plusminussigma==1:
    Type+="p1"
  if plusminussigma==-1:
    Type+="m1"
  ret=SMSObjects(Type,topo,ana,xvalue)
  if not extendedinfo:
    if not ret: return None
    if xvalue == '':
      key = '050'
    else:
      key = xvalue
    return ret[1][key]
  return ret

def efficiency ( topo,analysis,xfrac,histo='efficiency' ):
  return SMSObjects ( histo, topo, analysis, xfrac )

def getUpperLimit ( topo, analysis, xfrac, mx=None, my=None ):
  """ get the upper limit from an upper limit TH2F, query once
      specific bin given by mx and my .
      Return the entire histogram if mx,my==None,None
  """
  histo=SMSObjects("limit",topo,analysis,xfrac)
  if verbose:
    print "[SMSResultsCollector.getUpperLimit] histo=",histo


  if histo==None: return None
  tmp=histo[1]
  keys=tmp.keys()
  h,c=None,None
  if keys.count(xfrac)==1:
    h=tmp[xfrac]
  else:
    if len(keys)==1:
      h=tmp[keys[0]]
  if h:
    if mx==None and my==None: return h
    xbin=h.GetXaxis().FindBin(mx)
    ybin=h.GetYaxis().FindBin(my)
    c=h.GetBinContent(xbin,ybin)
    #if verbose:
    #  print "[SMSResultsCollector] upper limit=",xbin,ybin,c
  return c

