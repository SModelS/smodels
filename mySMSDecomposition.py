#!/usr/bin/python


import math
import SMSglobals, SMSanalyses, sys, SMSmethods

#PYTHIA must have MSTP(42)=0 ! no mass smearing (narrow width approximation)
#Initialize global variables:
SMSglobals.initglob()
#Creat analyses list:
SMSanalyses.load()



def getT1(b1,b2):
	pts=[]
	pts.extend(b1[2])
	pts.extend(b2[2])
	if pts.count('b')==4:
		return "T1bbbb"
	elif pts.count('t+')+pts.count('t-')==4:
		return "T1tttt"
	else: return "T1"

def getT3(b1,b2):
	pts=[]
        pts.extend(b1[2])
        pts.extend(b2[2])
        if pts.count('W+')+pts.count('W-')==1 and b2[2].count('b')==2:
                return "T3wb"
        elif pts.count('W+')+pts.count('W-')==1:
                return "T3w"
	elif pts.count('Z')==1:
		return "T3z"
        else: return "T3?"

def getT5(b1,b2):
	pts=[]
        pts.extend(b1[2])
        pts.extend(b2[2])
        if pts.count('photon')==2:
                return "T5gg"
        elif pts.count('W+')+pts.count('W-')==2:
                return "T5ww"
	elif pts.count('W+')+pts.count('W-')==1 and pts.count('Z')==1:
		return "T5wz"
	elif pts.count('Z')==2:
		return "T5zz"
	elif pts.count('l+')+pts.count('l-')==2 and pts.count('nu')==2:
		return "T5lnu"
        else: return "T5?"

def getT2(b1,b2):
        pts=[]
        pts.extend(b1[2])
        pts.extend(b2[2])
        if pts.count('b')==2:
                return "T2bb"
	elif pts.count('b')==1 and pts.count('W+')+pts.count('W-')==1:
		return "T2bw"
        elif pts.count('t+')==1 and pts.count('t-')==1:
                return "T2tt"
        else: return "T2"

def getT4(b1,b2):
	return "T4"

def getT6(b1,b2):
        pts=[]
        pts.extend(b1[2])
        pts.extend(b2[2])
        if pts.count('b')==2 and pts.count('Z')==2:
                return "T6bbzz"
        else: return "T6?"

def getTChi(b1,b2):
        pts=[]
        pts.extend(b1[2])
        pts.extend(b2[2])
        if pts.count('W+')+pts.count('W-')==1 and pts.count('Z')==1:
                return "TChiwz"
	elif pts.count('W+')==1 and pts.count('W-')==1:
		return "TChiww"
	elif pts.count('Z')==2:
		return "TChizz"
	elif pts.count('l+')+pts.count('l-')==3 and pts.count('nu')==1:
		return "TChiNuSlep" #same as TChiChipmSnuSlep?
	elif pts.count('l+')+pts.count('l-')==2 and pts.count('nu')==2:
		return "TChipmSnuSlep"
        else: return "TChi?"

def getTSlep(b1,b2):
        pts=[]
        pts.extend(b1[2])
        pts.extend(b2[2])
	if pts.count('l+')+pts.count('l-'):
		return "TSlepSlep"


def getSMS(input_lhe):
#Get branches from lhe-file
	File=open("regression/%s_1.lhe" % input_lhe)

	PList = SMSmethods.getNextEvent(File)
#weight = [pytres["xsecfb"]/float(nevts),pytres["xsecfb"]/float(nevts)]
	SMSmethods.GTop()
	SMSTop = SMSmethods.getEventTop(PList, {})
#Sort branches
	Branch1 = SMSTop[0].B[0]
	Branch2 = SMSTop[0].B[1]
	if (Branch1.vertnumb < Branch2.vertnumb) or (Branch1.vertnumb == Branch2.vertnumb and sum(Branch1.vertparts) < sum(Branch2.vertparts)):
      		SMSTop[0].B = [Branch2,Branch1]
#read branch infortmation
	b1=[ SMSTop[0].B[0].vertnumb,SMSTop[0].B[0].vertparts,SMSTop[0].B[0].ElList[0].particles,abs(SMSTop[0].B[0].ElList[0].momID) ]
	b2=[ SMSTop[0].B[1].vertnumb,SMSTop[0].B[1].vertparts,SMSTop[0].B[1].ElList[0].particles,abs(SMSTop[0].B[1].ElList[0].momID) ]
#find SMS
	if b1[0]==1 or b2[0]==1:
		return "bad input lhe"
	if b1[3]==1000021 and b2[3]==1000021:
		if b1[1]==[2,0] and b2[1]==[2,0]:
			return getT1(b1,b2)
		if (b1[1]==[2,2,0] or b1[1]==[2,1,0]) and b2[1]==[2,0]:
			return getT3(b1,b2)
		if (b1[1]==[2,2,0] and b2[1]==[2,2,0]) or (b1[1]==[2,1,0] and b2[1]==[2,1,0]):
       			 return getT5(b1,b2)
		return "T(odd)?"
	if (1000001<=b1[3]<=1000006 or 2000001<=b1[3]<=2000006) and (1000001<=b2[3]<=1000006 or 2000001<=b2[3]<=2000006):
	        if b1[1]==[1,0] and b2[1]==[1,0]:
                        return getT2(b1,b2)
		if b1[1]==[1,1,0] and b2[1]==[1,0]:
			return getT4(b1,b2)
                if b1[1]==[1,1,0] and b2[1]==[1,1,0]:
                         return getT6(b1,b2)
                return "T(even)?"
	if (1000022<=b1[3]<=1000025 or b1[3]==1000035 or b1[3]==1000037) and (1000022<=b2[3]<=1000025 or b2[3]==1000035 or b2[3]==1000037):
		return getTChi(b1,b2)
	if ((1000001<=b1[3]<=1000006 or 2000001<=b1[3]<=2000006) and b2[3]==1000021) or ((1000001<=b2[3]<=1000006 or 2000001<=b2[3]<=2000006) and b1[3]==1000021):
		return "TGQ"
	if b1[3] and b2[3] in [1000011, 1000013, 1000015, 2000011, 2000013, 2000014]:
		return getTSlep(b1,b2)
	return "Topo unknown"


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
