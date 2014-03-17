#!/usr/bin/env python

"""
.. module:: SMSResultsCollector
    :synopsis: Getting information from the results database, should not be used directly. Use SMSResults.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Doris Proschofsky <Doris.Proschofsky@assoc.oeaw.ac.at>

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

