#!/usr/bin/python

# (c) Wolfgang Waltenberger, 2012, 2013

""" the stuff that classifies an event by SMSes """

import math

def prodmode ( n ):
  if n["pgluino"]==2: return "gg"
  nallpsquarks=n["psquark"]+n["psbottom"]+n["pstop"]
  if n["pgluino"]==1 and nallpsquarks==1: return "sg"
  if n["pstop"]==2: return "tt"
  # fixme whats tb? whats lq,le,xx,hh,ht ?
  if n["psbottom"]==2: return "bb"
  if n["psbottom"]==1 and n["stop"]==1: return "tb"
  if n["psquark"]==2: return "ss"
  if (n["pslepton"] + n["psneutrino"])==2: return "ll"
  nallweakinos=n["pC1"]+n["pN2"]+n["pN1"]+n["pC2"]+n["pN3"]
  if n["pgluino"]==1 and nallweakinos==1: return "ng"
  if n["psquark"]==1 and nallweakinos==1: return "ns"
  if nallweakinos==2: return "nn"
  if nallpsquarks==2: return "ss"
  if nallpsquarks==1 and nallweakinos==1: return "ns"
  return "??"

prodModeDescription={ "gg": "gluino-gluino", "sg": "gluino-squark (any flavor)",
  "tt": "stop-(anti)stop", "bb": "sbottom-(anti)sbottom", "tb": "stop-sbottom",
  "ss": "squark-squark (light flavors, or light+heavy flavor)", "ll": "sleptons/sneutrinos", "ng": "chargino/neutralino + gluino", "ns": "chargino/neutralino + squark (any flavor)", "nn": "charginos/neutralinos", "??": "I have no idea" } 

def extendedProdMode ( n ):
  """ currently not used """
  if n["pgluino"]==2: return "gg"
  nallpsquarks=n["psquark"]+n["psbottom"]+n["pstop"]
  if n["pgluino"]==1 and n["psquark"]==1: return "sg"
  if n["pgluino"]==1 and n["psbottom"]==1: return "bg"
  if n["pgluino"]==1 and n["pstop"]==1: return "tg"
  if n["pstop"]==2: return "tt"
  # fixme whats tb? whats lq,le,xx,hh,ht ?
  if n["psbottom"]==2: return "bb"
  if n["psbottom"]==1 and n["stop"]==1: return "tb"
  if n["psquark"]==2: return "ss"
  if n["pslepton"]==2: return "ll"
  if n["pN3"]==2: return "N3N3"
  if n["pN2"]==2: return "N2N2"
  if n["pN1"]==2: return "N1N1"
  if n["pC2"]==2: return "C2C2"
  if n["pC1"]==2: return "C1C1"
  nallweakinos=n["pC1"]+n["pN2"]+n["pN1"]+n["pC2"]+n["pN3"]
  if n["pgluino"]==1 and nallweakinos==1: return "ng"
  if n["psquark"]==1 and nallweakinos==1: return "ns"
  if nallweakinos==2: return "nn"
  if nallpsquarks==2: return "ss"
  if nallpsquarks==1 and nallweakinos==1: return "ns"
  return "??"

def gluinoSMSes ( n ):
  weakinos=n["C1"]+n["N2"]+n["C2"]+n["N3"]
  if weakinos==0 and n["t"]==4:
    return "T1tttt"
  if weakinos==0 and n["b"]==4 and n["W"]==4:
    return "T1tttt*"
  if n["t"]==4: return "T1tttt+"
  if weakinos==0 and n["b"]==4:
    return "T1bbbb"
  if weakinos==2 and n["N2"]==2:
    if n["Z"]==2:
      return "T5zz"
    if n["photon"]==2:
      return "T5gg"
    return "T5NN"
  if n["higgs"]==2 and weakinos == 2 and n["N2"]==2:
      return "T5hh"
  if n["C1"]==1 and n["N2"]==1 and weakinos==2:
    if n["Z"]==1 and n["W"]==1:
      return "T5wz"
    if n["Z"]==1:
      return "T5Cz"
    if n["W"]==1:
      return "T5wN"
    return "T5CN"
  if n["C1"]==2 and n["N2"]==0:
    if n["W"]==2:
      return "T5ww"
    return "T5CC"
  if n["N2"]==1 and n["C1"]==0 and n["higgs"]==0:
    if n["Z"]==1:
      return "T3z"
    return "T3N"
  if n["C1"]==1 and n["N2"]==0 and n["higgs"]==0:
    if n["W"]==1:
      if n["b"]==2:
        return "T3wb"
      return "T3w"
    return "T3C"
  if weakinos==1 and n["N2"]==1  and n["Z"]==0 and n["higgs"]==1:
    return "T3h"
  if weakinos==0 and n["higgs"]==0 and n["squark"]==0 and n["sbottom"]==0 and n["stop"]==0:
    return "T1"
  if weakinos==0 and n["stop"]==2:
    if n["t"]==4:
      return "T5tttt"
    return "T5tt+"
  return "T1+"

def squarkSMSes ( n ):
  nPs=n["jet"]+n["neutrino"]+n["electron"]+n["muon"]+n["tau"]
  if n["C2"]==2 and n["N3"]==0:
    return "T6C2C2"
  if n["N3"]==2 and n["C2"]==0:
    return "T6N3N3"
  if n["N3"]==1 and n["N2"]==1:
    return "T6N3N2"
  if n["C2"]==1 and n["N2"]==1:
    return "T6C2N2"
  if n["C1"]==2 and n["N2"]==0:
    if n["W"]==2:
      return "T6ww"
    if nPs==6:
      return "T6CC"
    if nPs>6:
      return "T6CC*"
    return "T6CC?"
  if n["C1"]==0 and n["N2"]==2:
    if n["Z"]==2:
      return "T6zz"
    if nPs==6:
      return "T6N2N2"
    if nPs>6:
      return "T6N2N2*"
    return "T6N2N2?"
  if ( n["C1"] + n["N2"] ) ==2 \
      and n["higgs"]==2:
    return "T6hh"
  if n["C1"]==1 and n["N2"]==1:
    if n["Z"]==1 and n["W"]==1:
      return "T6wz"
    if nPs==6:
      return "T6C1N2"
    if nPs>6:
      return "T6C1N2*"
    return "T6C1N2?"
  if n["C1"]==0 and n["N2"]==1 and n["Z"]==0 and n["higgs"]==1:
    return "T4h"
  if n["C1"]==0 and n["N2"]==1 and n["higgs"]==0:
    if n["Z"]==1:
      return "T4z"
    if nPs==4:
      return "T4N2"
    if nPs>4: 
      return "T4N2*"
    return "T4N2?"
  if n["C1"]==1 and n["N2"]==0 and n["C2"]==0 and n["N3"]==0:
    if n["W"]==1:
      return "T4w"
    if nPs==4:
      return "T4C1"
    if nPs>4: 
      return "T4C1*"
    return "T4C1?"
  weakinos=n["C1"]+ n["N2"] + n["C2"] + n["N3"]
  if weakinos==0 and n["higgs"]==0:
    return "T2"
  if weakinos == 2:
    return "T6"
  if weakinos == 1:
    return "T4"
  if weakinos > 2:
    return "T8"
  return "T2general"

def associateSMSes (n ):
  m0=0
  m12=0
  if n.has_key ( "msugra_m0" ):
    m0=n["msugra_m0"]
    m12=n["msugra_m12"]
  mgl=2.8 * m12
  msq=math.sqrt ( m0**2 + 5.5 * m12**2 )
  flavor=""
  if n["pstop"]==1: flavor="t"
  if n["psbottom"]==1: flavor="b"
  if n["C1"]==0 and n["N2"]==0:
    if mgl >= msq: return "TGQ%s" % flavor
    if mgl < msq: return "TGQ2%s" % flavor
    return "TGQ%s?" % flavor
  if ( n["C1"] + n["N2"] )>0:
    if mgl >= msq: return "TGQ%sinter" % flavor
    if mgl < msq: return "TGQ2%sinter" % flavor
  return "TGQ%sincl" % flavor

def weakinoSMSes ( n ):
  try:
    if n["pN1"]==2:
      return "TChiN1N1"
    if n["pC1"]==1 and n["pN2"]==1 and n["W"]==1 and n["Z"]==1 and n["slepton"]==0:
      return "TChiwz"
    if n["pC1"]==2 and n["pN2"]==0 and n["W"]==2 and n["slepton"]==0:
      return "TChiww"
    if n["pC1"]==0 and n["pN2"]==2 and n["Z"]==2 and n["slepton"]==0:
      return "TChizz"
    if n["pC1"]==0 and n["pN2"]==2 and n["higgs"]==2 and n["slepton"]==0:
      return "TChihh"
    if n["pC1"]==1 and n["pN2"]==1 and n["higgs"]==1 and n["W"]==1  and n["slepton"]==0:
      return "TChiwh"
    if n["pC1"]==1 and n["pN2"]==1 and n["higgs"]>=1 and n["W"]>=1  and n["slepton"]==0:
      return "TChiwh+"
    if n["pC1"]==0 and n["pN2"]==2 and n["higgs"]==1 and n["Z"]==1 and n["slepton"]==0:
      return "TChizh"
    if n["pC1"]==0 and n["pN2"]==2 and n["higgs"]>=1 and n["Z"]>=1 and n["slepton"]==0:
      return "TChizh+"
    nPs=n["jet"]+n["neutrino"]+n["electron"]+n["muon"]+n["tau"]+n["b"]+n["t"]
    if n["pC1"]==0 and n["pN2"]==2 and n["slepton"]==0 and n["sneutrino"]==0:
      if nPs==4:
        return "TChiN2N2"
      if nPs>4:
        return "TChiN2N2*"
      return "TChiN2N2?"
    if n["pC1"]==1 and n["pN2"]==1 and n["slepton"]==0 and n["sneutrino"]==0:
      if nPs==4:
        return "TChiN2C1"
      if nPs>4:
        return "TChiN2C1*"
      if n["Z"]>=1 and n["W"]>=1:
        return "TChiwz"
      if n["Z"]>=1:
        return "TChiC1z"
      if n["W"]>=1:
        return "TChiN2w"
      if n["higgs"]>=1:
        return "TChiC1h"
      return "TChiN2C1?"
    if n["pC1"]==2 and n["pN2"]==0 and n["slepton"]==0 and n["sneutrino"]==0:
      if nPs==4:
        return "TChiC1C1"
      if nPs>4:
        return "TChiC1C1*"
      if n["W"]==2:
        return "TChiww"
      if n["W"]==1:
        return "TChiC1w"
      if (n["N2"]+n["C1"])>2:
        return "TChiC1C1+"
      return "TChiC1C1?"
    if n["pN1"]==1 and n["pN2"]==1 and n["higgs"]==1  and n["slepton"]==0:
      return "TChihN1"
    if n["pC1"]==1 and n["pN1"]==1 and n["W"]==1  and n["slepton"]==0:
      return "TChiwN1"
    if n["pC1"]==0 and n["pN2"]==0 and n["pN1"]==2  and n["slepton"]==0:
      return "TChiN1N1"
    if n["pN2"]==1 and n["pN1"]==1 and n["Z"]==1  and n["slepton"]==0:
      return "TChizN1"
    if n["pC1"]==1 and n["pN1"]==1  and n["slepton"]==0 and n["sneutrino"]==0:
      #print "TChiC1N1",n,"nPs=",nPs
      if nPs==2 or nPs==0:
        return "TChiC1N1"
      if nPs==1:
        return "TChiC1N1+"
      return "TChiC1N1*"
    if n["pN2"]==1 and n["pN1"]==1  and n["slepton"]==0 and n["sneutrino"]==0:
      #print "TChiN2N1",n,"nPs=",nPs
      if nPs==2 or nPs==0:
        return "TChiN2N1"
      if nPs==1:
        return "TChiN2N1+"
      return "TChiN2N1*"
    weakinos=n["pN2"] + n["pC1"] + n["pN3"] + n["pC2"]
    sleptons=n["slepton"] + n["sneutrino"] 
    if n["pN2"]==2 and n["slepton"]==2:
      return "TChiSlepSlep"
    if n["pN2"]==2 and n["sneutrino"]==2:
      return "TChiSnuSnu"
    if n["pN2"]==2 and n["sneutrino"]==1 and n["slepton"]==1:
      return "TChiSnuSlep"
    if n["pN2"]==1 and n["pC1"]==1 and n["slepton"]==1 and n["sneutrino"]==1:
      return "TChiChipmSnuSlep"
    if n["pN2"]==1 and n["pC1"]==1 and n["slepton"]==2 and n["sneutrino"]==0:
      return "TChiChipmSlepSlep"
    if n["pC1"]==2 and n["slepton"]==2:
      return "TChipmSlepSlep"
    if n["pC1"]==2 and n["sneutrino"]==2:
      return "TChipmSnuSnu"
    if n["pC1"]==2 and sleptons==2:
      return "TChipmSnuSlep"
    if ( n["pN2"] + n["pC1"])==2 and ( n["slepton"] + n["sneutrino"] )==2: 
      return "TChiSlepSlep+"
    if n["pN3"]>0 and n["pC1"]>0:
      return "TChiN3C1"
    if n["pN3"]>0 and n["pN2"]>0:
      return "TChiN3N2"
    if n["pN2"]>0 and n["pC2"]>0:
      return "TChiN2C2"
    if n["pN3"]>0 and n["pC2"]>0:
      return "TChiN3C2"
    if n["pN1"]>0 and n["pC2"]>0:
      return "TChiN1C2"
    if n["pC1"]>0 and n["pC2"]>0:
      return "TChiC2C1"
    if n["pN2"]==1 and n["pC1"]==1 and ( n["slepton"] + n["sneutrino"] )==1: 
      return "TChiN2C1Slep"
    if n["pN2"]==1 and n["pC1"]==1:
      return "TChiN2C1+"
    if n["pC1"]==2 and ( n["slepton"] + n["sneutrino"] )==1: 
      return "TChiC1C1Slep"
    if n["pC1"]==2:
      return "TChiC1C1+"
    if n["pN3"]>0:
      if n["pN3"]==2:
        return "TChiN3N3"
      return "TChiN3+"
    if n["pC2"]>0:
      if n["pC2"]==2:
        return "TChiC2C2"
      if n["pC2"]==1:
        return "TChiC2-"
      return "TChiC2+"
    return "TChiC2%dC1%dN3%dN2%dN1%d" % ( n["pC2"], n["pC1"], n["pN3"], n["pN2"], n["pN1"] )
  except Exception:
    print "weakinoSMS, something went wrong: n=",n
    sys.exit(0)

def stopSMSes ( n ):
  if n["C1"]==0 and n["N2"]==0 and n["t"]==2:
    return "T2tt"
  if n["C1"]==0 and n["N2"]==2: ## and n["t"]==2:
    if n["Z"]==2:
      return "T6ttzz"
    else:
      return "T6ttNN"
  if n["C1"]==0 and n["N2"]==1: ## and n["t"]==1:
    if n["Z"]==1:
      return "T4ttz"
    else:
      return "T4ttN"
  if n["C1"]==1 and n["N2"]==0 and n["b"]==1:
    if n["W"]==1:
      return "T4ttw"
    else:
      return "T2ttC"
  if n["C1"]==1 and n["N2"]==1: ## and n["b"]==2 and n["t"]==1:
    if n["W"]==1 and n["Z"]==1:
      return "T6ttwz"
    else:
      return "T6ttCN"
  if n["C1"]==2 and n["N2"]==0 and n["b"]==2:
    if n["W"]==2:
      return "T6ttww"
    else:
      return "T6ttCC"
  if n["C1"]==0 and n["N2"]==0 and n["C2"]==0 and n["N3"]==0 and n["W"]==2 and n["b"]==2:
    return "T2ttWbWb" ## ~t -> W b chi10
  if n["C1"]==0 and n["N2"]==0 and n["W"]==2 and n["b"]==2:
    return "T6ttWbWb" ## ~t -> W b chi10
  ##return "T2ttincl%d%d%d%d" % ( n["N2"], n["C1"], n["t"], n["b"] )
  return "T2ttincl"

def sbottomSMSes ( n ):
  try:
    if n["C1"]==0 and n["N2"]==0 and n["higgs"]==0 and n["b"]==2:
      return "T2bb"
    if n["C1"]==0 and n["N2"]==2 and n["higgs"]==0 and n["b"]==2:
      if n["Z"]==2:
        return "T6bbzz"
      else:
        return "T6bbN2N2"
    if n["C1"]==2 and n["N2"]==0 and n["higgs"]==0 and n["b"]==4:
      if n["W"]==2:
        return "T6bbww"
      else:
        return "T6bbC1C1"
    if n["C1"]==1 and n["N2"]==1 and n["higgs"]==0 and n["b"]==3:
      if n["W"]==1 and n["Z"]==1:
        return "T6bbwz"
      else:
        return "T6bbC1N2"
    if (n["C1"] + n["N2"] + n["higgs"])==2 and n["b"]==2:
      return "T6bb*"
    if (n["N2"] + n["higgs"])==1 and n["b"]==2 and n["Z"]==1:
      return "T4bbz"
    if (n["C1"] + n["higgs"])==1 and n["b"]==2 and n["W"]==1:
      return "T4bbw"
    if n["N2"]==1 and n["b"]==2:
      return "T4bbN2"
    if (n["C1"] + n["higgs"])==1 and n["b"]==2:
      return "T4bbC1"
    if (n["C1"] + n["N2"] + n["higgs"])==1 and n["b"]==2:
      return "T4bb*"
    #if n["b"]!=2:
    #return "T2bbFlavorViolating"
    return "T2bbincl"
  except Exception:
    pass
  return "except"

def coloredweakSMSes ( n ):
  if n["pgluino"]==1 and ( n["pC1"] + n["pN2"] + n["pN1"] )==1:
    if n["W"]==1: return "TNGw"
    if n["Z"]==1: return "TNGz"
    return "TNG"
  if n["psquark"]==1 and ( n["pC1"] + n["pN2"] + n["pN1"] )==1:
    if n["W"]==1: return "TNSw"
    if n["Z"]==1: return "TNSz"
    return "TNS"
  if n["pstop"]==1 and ( n["pC1"] + n["pN2"] + n["pN1"] )==1:
    if n["W"]==1: return "TNTw"
    if n["Z"]==1: return "TNTz"
    return "TNT"
  if n["psbottom"]==1 and ( n["pC1"] + n["pN2"] + n["pN1"] )==1:
    if n["W"]==1: return "TNBw"
    if n["Z"]==1: return "TNBz"
    return "TNB"
  return "othercoloredweak"

def sleptonSMSes ( n ):
  pweakinos=n["pC1"]+n["pN2"]+n["pN1"]+n["pC2"]+n["pN3"]
  psleptons=n["pslepton"]+n["psneutrino"]
  pcolored=n["pgluino"]+n["psquark"]+n["psbottom"]+n["pstop"]
  leptons=n["electron"]+n["muon"]+n["tau"]
  if n["pslepton"]==2 and pweakinos==0 and pcolored==0 and leptons==2:
    return "TSlepSlep"
  if n["pslepton"]==2:
    return "TSlepSlep+"
  if n["psneutrino"]==2:
    return "TSnuSnu+"
  if n["psneutrino"]==1 and n["pslepton"]==1:
    return "TSlepSnu+"
  return "TSlep*"

def getSMS ( n, lsp="N1" ):
  if lsp!="N1":
    return "N1NotLSP"
  if not n.has_key ( lsp ) or n[lsp]!=2: return "NotTwoLSPs"
  pweakinos=n["pC1"]+n["pN2"]+n["pN1"]+n["pC2"]+n["pN3"]
  psleptons=n["pslepton"]+n["psneutrino"]
  pcolored=n["pgluino"]+n["psquark"]+n["psbottom"]+n["pstop"]
  pallsquarks=n["psquark"]+n["psbottom"]+n["pstop"]
  allsusy=pcolored+pweakinos+psleptons
  if allsusy==0: return "nonsusy"
  if n["pgluino"]==2 and n["psquark"]==0: return gluinoSMSes ( n )
  if n["psquark"]==2 and n["pgluino"]==0: return squarkSMSes ( n )
  if psleptons==2: return sleptonSMSes ( n )
  if pallsquarks==1 and n["pgluino"]==1: return associateSMSes ( n )
  if n["pgluino"]==0 and n["pstop"]==2: return stopSMSes (n )
  if n["pgluino"]==0 and n["psbottom"]==2: return sbottomSMSes (n )
  if pcolored==0 and pweakinos==2: return weakinoSMSes(n)
  if pcolored>2: return "multicol"
  if pcolored==0: return "noncolored"
  if pcolored>0 and pweakinos > 0: return coloredweakSMSes(n)
  if pcolored > 0: return "othercolored"
  return "othersusy"

def describe ( n ):
    s=""
    for (key,value) in n.items():
      if value>0:
        if key[:1]=="p":
          s="%s:%d %s" % ( key,value,s)
        else:
          s+="%s:%d " % ( key,value)
    return s

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
  if p==22: return "photon"
  if p==24: return "W"
  if p==16: return "neutrino"
  if p==15: return "tau"
  if p==14: return "neutrino"
  if p==13: return "muon"
  if p==12: return "neutrino"
  if p==11: return "electron"
  if p==5: return "b"
  if p==6: return "t"
  if p==21: return "gluon"
  if p>=1 and p<=4: return "jet"
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

countingvariables=[ "jet", "electron", "muon", "gluino", "squark",
  "gluon", "stop", "squarkbar", "sbottom", "N1", "photon",
  "N2", "pN2", "C1", "pC1", "pC2", "pN3", "C2", "N3",
  "Z", "W", "t", "b", "higgs", "pgluino", "psquark", "tau",
  "pN1", "slepton", "sneutrino", "msugra_m12", "msugra_m0", "?s",
  "psbottom","pstop", "pslepton", "neutrino", "phiggs", "psneutrino"
  ]
