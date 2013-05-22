#!/usr/bin/python
import sys
sys.path.append ( "../" )

import math

def ptsCount(b1,b2):
	pts=[]
	b1.extend(b2)
	for b in b1:
		pts.extend(b)
	return pts

def getT1(pts):
	if pts.count('b')==4:
		return "T1bbbb"
	elif pts.count('t+')+pts.count('t-')==4:
		return "T1tttt"
	elif pts.count('b')==2 and pts.count('t+')+pts.count('t-')==2:
		return "T1tbtb"
	else: return "T1"

def getT3(pts, b2):
        if pts.count('W+')+pts.count('W-')==1 and b2[2][0].count('b')==2:
                return "T3Wb"
        elif pts.count('W+')+pts.count('W-')==1:
                return "T3W"
	elif pts.count('Z')==1:
		return "T3Z"
	elif (pts.count('e+')==1 and pts.count('e-')==1) or (pts.count('mu+')==1 and pts.count('mu-')==1):
		return "T3lh"
	elif (pts.count('e+')+pts.count('e-')==1 or pts.count('mu+')+pts.count('mu-')==1) and pts.count('nu')==1:
		return "T3lnu"
        else: return "T3?"

def getT5(pts):
        if pts.count('photon')==2:
                return "T5gg"
        elif pts.count('W+')+pts.count('W-')==2:
                return "T5WW"
	elif pts.count('W+')+pts.count('W-')==1 and pts.count('Z')==1:
		return "T5WZ"
	elif pts.count('Z')==2:
		return "T5ZZ"
	elif pts.count('l+')+pts.count('l-')==2 and pts.count('nu')==2:
		return "T5lnu"
	elif pts.count('t+')+pts.count('t-')==4:
		return "T5tttt"
        else: return "T5?"

def getT2(pts):
        if pts.count('b')==2:
                return "T2bb"
	elif pts.count('b')==1 and pts.count('W+')+pts.count('W-')==1:
		return "T2bW" #gibt es so nicht, sollte T6bbww geben?
        elif pts.count('t+')==1 and pts.count('t-')==1:
                return "T2tt"
        elif pts.count('e+')+pts.count('e-')==2 or pts.count('mu+')+pts.count('mu-')==2:
                return "TSlepSlep"
        elif pts.count('W+')+pts.count('W-')==1 and pts.count('Z')==1:
                return "TChiWZ"
        elif pts.count('W+')==1 and pts.count('W-')==1:
                return "TChiWW"
        elif pts.count('Z')==2:
                return "TChiZZ"

	else: return "T2"

def getT4(pts):
	return "T4"

def getT6(pts):
        if pts.count('b')==2 and pts.count('Z')==2:
                return "T6bbZZ"
	elif pts.count('l+')+pts.count('l-')==3 and pts.count('nu')==1:
                return "TChiNuSlep" #same as TChiChipmSnuSlep?
        elif pts.count('l+')+pts.count('l-')==2 and pts.count('nu')==2:
                return "TChipmSnuSlep"
	elif pts.count('t')==2 and pts.count('W')==2:
		return "T6ttWW"
	elif pts.count('b')==2 and pts.count('W')==2:
		return "T6bbWW"
	elif pts.count('t')==2 and pts.count('Z')==2:
		return "T6ttZZ"
        else: return "T6?"



def getTx(element):
	"""takes EElement and returns Tx-Name"""
	import copy
	E=copy.deepcopy(element)
	Branch1=E.B[0]
	Branch2=E.B[1]
	d=E.getEinfo()
	v1=d["vertnumb"][0]
	v2=d["vertnumb"][1]
	p1=d["vertparts"][0]
	p2=d["vertparts"][1]
#Sort branches
	if v1 < v2 or (v1==v2  and sum(p1)<sum(p2)):
		Branch1, Branch2 = Branch2, Branch1
		v1, v2 = v2, v1
		p1, p2 = p2, p1
#read branch infortmation
	b1=[ v1,p1,Branch1.particles ]
	b2=[ v2,p2,Branch2.particles ]
#find SMS
	if b1[0]==1 or b2[0]==1:
		return "bad input"
	if b1[1]==[2,0] and b2[1]==[2,0]:
		return getT1(ptsCount(b1[2],b2[2]))
	if (b1[1]==[2,2,0] or b1[1]==[2,1,0]) and b2[1]==[2,0]:
		return getT3(ptsCount(b1[2],b2[2]), b2)
	if (b1[1]==[2,2,0] and b2[1]==[2,2,0]) or (b1[1]==[2,1,0] and b2[1]==[2,1,0]):
       		return getT5(ptsCount(b1[2],b2[2]))
	if b1[1]==[1,0] and b2[1]==[1,0]:
                return getT2(ptsCount(b1[2],b2[2]))
	if b1[1]==[1,1,0] and b2[1]==[1,0]:
		return getT4(ptsCount(b1[2],b2[2]))
        if b1[1]==[1,1,0] and b2[1]==[1,1,0]:
                return getT6(ptsCount(b1[2],b2[2]))
	return "?"

'''
def particle ( pid ):
  p=int(abs(pid))
  if p==1000021: return "gluino"
  if p==1000022: return "N1"
  if p==1000023: return "N2"
  if p==1000025: return "N3"
  if p==1000035: return "N3"
  if p==1000024: return "C1"
  if p==1000037: return "C2"
  if p==1000039: return "gravitino"
  if p>=1000001 and p<=1000004: return "squark"
  if p>=2000001 and p<=2000004: return "squark"
  if p==1000005: return "sbottom"
  if p==2000005: return "sbottom"
  if p==1000006: return "stop"
  if p==2000006: return "stop"
  if p==1000011 or p==1000013 or p==1000015: return "slepton"
  if p==2000011 or p==2000013 or p==2000015: return "slepton"
  if p==1000012 or p==1000014 or p==1000016: return "sneutrino"
  if p==2000012 or p==2000014 or p==2000016: return "sneutrino"
  if p==25: return "higgs"
  if p==35: return "higgs"
  if p==36: return "higgs"
  if p==37: return "higgs"
  if p==23: return "Z"
  if p==22: return "photons"
  if p==24: return "W"
  if p==16: return "neutrino"
  if p==15: return "tau"
  if p==14: return "neutrino"
  if p==13: return "muon"
  if p==12: return "neutrino"
  if p==11: return "electron"
  if p==5: return "b"
  if p==6: return "top"
  if p>=1 and p<=4: return "jets"
  return str(p)

allproduction=[ "ss", "gg", "??", "ng", "ns", "tb", "bb", "sb", "tt", "ll" ]

empty={ "T1general":0, "??":0, "ss":0, "gg": 0, "sg": 0, "nn": 0, "ng": 0,
         "ns": 0, "tb": 0, "bb": 0, "sb": 0, "?": 0, "tt": 0, "ll":0,
         "T1":0, "T2":0, "T3w":0, "T4w": 0., "T5zz": 0, "T5ww": 0, "T6zz": 0,
         "T1tttt": 0, "T2tt": 0, "T1ttttincl": 0, "T2ttincl": 0,
         "T5zwoff": 0, "T3zoff": 0, "T4zoff" : 0, "T2general": 0,
         "T3woff":0, "T4woff":0, "T5zzoff":0, "T5wwoff":0, "T6zzoff": 0,
         "T6wwoff":0, "T6zwoff":0, "T6ww": 0., "T6zw": 0., "other": 0,
         "othercolored":0,"multicol":0, "weakino": 0, "TGQ2": 0,
         "T4bb": 0, "T6bb": 0,
         "TGQinter": 0, "coloredweak":0, "T3h": 0, "T4h": 0,
         "T5h": 0, "T6h": 0, "empty": 0, "TGQt":0, "TGQtt":0, "TGQttt":0,
         "TGQtstar":0, "TGQt": 0, "TGQ2inter": 0, "TGQ?": 0, "base":0,
         "TChiwn": 0, "TChizn": 0, "TChinn":0, "TChizh": 0, "TChiwh": 0,
         "T5zw":0, "T3z":0, "T4z":0, "noncolored":0, "TGQ": 0, "total":0 }

countingvariables=[ "jets", "electron", "muon", "gluino", "squark",
  "stop", "squarkbar", "sbottom", "N1", "photons",
  "N2", "pN2", "C1", "pC1", "pC2", "pN3", "C2", "N3",
  "Z", "W", "t", "b", "higgs", "pgluino", "psquark", "tau",
  "pN1", "slepton", "sneutrino", "msugra_m12", "msugra_m0", "?s",
  "psbottom","pstop", "pslepton", "neutrino", "phiggs", "psneutrino"
  ]
'''
