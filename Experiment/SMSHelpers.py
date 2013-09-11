#!/usr/bin/env python

"""
.. module:: SMSHelpers
        :synopsis: Some helper functions that are not to be used by the end user.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>, Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>, Doris Proschofsky <Doris.Proschofsky@assoc.oeaw.ac.at>

"""


import os
from Tools.PhysicsUnits import rmvunit
from Experiment.experimentExceptions import MetaInfoError

Base = "/afs/hephy.at/user/w/walten/public/sms/"

runs=[ "8TeV", "2012", "ATLAS8TeV", "2011", "RPV8", "RPV7" ]
## runs=[ "2012" ]

verbose=True

## track the open root files
openFiles={}

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

runs_={}
def getRun ( analysis, run=None ):
    """ search for an analysis, return the run,
            or return None if not found anywhere.
            if a specific run is given, then just check
            if we really have results. """
    key=analysis+str(run)
    if runs_.has_key ( key ): return runs_[key]
    if run:
        if os.path.exists ( "%s/%s/%s" % ( Base, run, analysis ) ):
            runs_[key]=run
            return run
#        else:
#            log ( "dont know about %s/%s" % (run,analysis), "warning" )
#            return None
    for trun in runs:
        if os.path.exists ( "%s/%s/%s" % ( Base, trun, analysis ) ):
            runs_[key]=trun
            return trun
    runs_[key]=None
    return None

pMI_={}
def parseMetaInfo ( analysis, run ):
    """ get all the meta information for a given analysis/run pair """
    key=analysis+run
    if pMI_.has_key ( key ): return pMI_[key]
    info="%s/%s/%s/info.txt" % ( Base, run, analysis )
    ret={}
    if not os.path.exists ( info ):
        log ( "cannot find %s" % info, "warning" )
        pMI_[key]=ret
        return ret
    f=open(info)
    lines=f.readlines()
    f.close()
    for line in lines:
        if line.find("#")>-1:
            line=line[:line.find("#")].strip()
        line=line.replace("\n","")
        line=line.strip()
        if line=="": continue
        tokens=line.split(":",1)
        if not len(tokens)==2:
            log ( "[117] cannot parse this line (1): ``%s'' in ``%s''" % (line, info) )
            continue
        if tokens[0]=="exclusions":
            # we treat these separately
            continue
        ret[tokens[0]]=tokens[1]
    pMI_[key]=ret
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
        if line.find("#")>-1:
            line=line[:line.find("#")].strip()
        if line=="": continue
        tokens=line.split(":",1)
        if not len(tokens)==2:
            log ( "[141] cannot parse this line (2): %s in %s" % (line, info) )
            continue
        if tokens[0]!="exclusions":
            # we're only interested in the exclusions
            continue
        excl=tokens[1]
        while excl[0]==" ": excl=excl[1:]
        if excl[-1]=='\n': excl=excl[:-1]
        excls=excl.split(" ")
        if len(excls)<2:
            log ( "[151] cannot parse this line (3): %s in %s" % (line, info) )
            continue
        values=map(int,excls[1:])
        if len(values)==2:
            values.insert(0,0)
        ret[excls[0]]=values
    return ret

mlines={}
def getLines ( analysis, run, label="condition" ):
    """ get all the conditions for a analysis/run pair """
    key=analysis+run+label
    if mlines.has_key ( key ): return mlines[key]
    info="%s/%s/%s/info.txt" % ( Base, run, analysis )
    ret={}
    if not os.path.exists ( info ):
        log ("cannot find %s" % info )
        mlines[key]=ret
        return ret
    f=open(info)
    lines=f.readlines()
    f.close()
    for line in lines:
        if line.find("#")>-1:
            line=line[:line.find("#")].strip()
        if line=="": continue
        tokens=line.split(":",1)
        if not len(tokens)==2:
            log ( "[172] cannot parse this line (2): %s in %s" % (line, info) )
            continue
        if tokens[0]!=label:
            # we're only interested in the conditions
            continue
        excl=tokens[1]
        while excl[0]==" ": excl=excl[1:]
        if excl[-1]=='\n': excl=excl[:-1]
        keyvalue=excl.split(" -> ")
        if len(keyvalue)!=2:
            log ( "[185] cannot parse the following line: %s" % keyvalue )
        ret[ keyvalue[0] ] = keyvalue[1]
        # ret.append(excl)
    mlines[key]=ret
    return ret

def conditions ( analysis, run ):
    """ get all the conditions for a analysis/run pair """
    return getLines ( analysis, run, "condition" )

def fuzzyconditions ( analysis, run ):
    """ get all the fuzzyconditions for an analysis/run pair """
    return getLines( analysis, run, "fuzzycondition" )

def constraints ( analysis, run ):
    """ get all the conditions for a analysis/run pair """
    return getLines ( analysis, run, "constraint" )

def getPotentialNames ( topo ):
    """ If T6bbww doesnt yield a result, try T6bbWW, etc """
    ret = [ topo ]
    ret.append ( topo.replace("w","W").replace("z","Z" ) )
    ret.append ( topo.replace("W","w").replace("Z","z" ) )
    return ret

def getCanonicalName ( topo ):
    """ define a canonical name: w and z's are uppercase letters, etc """
    topo=topo.replace("w","W").replace("z","Z" )
    return topo

def getRootFileName ( analysis, run=None ):
    """ returns the root file name, but does not check if it exists """
    if not run: run=getRun ( analysis )
    return "%s/%s/%s/sms.root" % ( Base, run, analysis )

upperLimitHisto={}
expupperLimitHisto={}
def getUpperLimitFromHisto ( analysis, topo, run, complain=False, expected=False ):
    key=analysis+topo+str(run)
    if expected:
        if expupperLimitHisto.has_key ( key ): return expupperLimitHisto[key] 
    else:
        if upperLimitHisto.has_key ( key ): return upperLimitHisto[key]
    import ROOT
    rootfile="%s/%s/%s/sms.root" % ( Base, run, analysis )
    if not os.path.exists(rootfile):
        log("root file %s doesnt exist" % rootfile, "warning" )
        upperLimitHisto[key]=None
        expupperLimitHisto[key]=None
        return None
    f=None
    if openFiles.has_key ( rootfile ):
        f=openFiles[rootfile]
    else:
        f=ROOT.TFile(rootfile)
        if not f or not f.IsOpen():
            log("root file %s cannot be opened" % rootfile )
            upperLimitHisto[key]=None
            expupperLimitHisto[key]=None
            return None
        openFiles[rootfile]=f
    if expected:
        histo=f.Get("expectedlimit_%s" % getCanonicalName(topo) )
        if not histo:
            if complain: log("expected upper limit histogram %s not found in %s" % ( topo, rootfile ))
            expupperLimitHisto[key]=None
            return None
        expupperLimitHisto[key]=histo
    else:
        histo=f.Get("limit_%s" % getCanonicalName(topo) )
        if not histo:
            if complain: log("histogram %s not found in %s" % ( topo, rootfile ))
            upperLimitHisto[key]=None
            return None
        upperLimitHisto[key]=histo
    return histo

def getUpperLimitAtPoint ( histo, mx, my, interpolate=False ):
    """ given already a histogram, we get the upper limit for mx and my,
            we return None if we're out of bounds """

    if rmvunit(my,'GeV')==None:
        log ( "inconsistent mx/my, mx=None, my=%s" % my )
        return None
    if not histo: return None ## 'no histogram'
    if interpolate: #for interpolation: set empty bins to the last finite value to avoid zeros in the interpolation, set result=None if the corresponding bin in the original histogram was empty
        hInf = histo.Clone()
        lastBinContent=0
        for bx in range(1,hInf.GetXaxis().GetNbins()+1):
            for by in range(1,hInf.GetYaxis().GetNbins()+1):
                if not hInf.GetBinContent(bx,by): hInf.SetBinContent(bx,by,lastBinContent)
                else: lastBinContent=hInf.GetBinContent(bx,by)
        if histo.GetBinContent(histo.GetXaxis().FindBin(rmvunit(mx,'GeV')),histo.GetYaxis().FindBin(rmvunit(my,'GeV'))):
            return hInf.Interpolate ( rmvunit(mx,'GeV'), rmvunit(my,'GeV') )
        else: return None
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

effhistos={}
def getEfficiencyHisto ( analysis, topo, run ):
    key=analysis+topo+str(run)
    if effhistos.has_key ( key ): return effhistos[key]
    import ROOT
    rootfile="%s/%s/%s/sms.root" % ( Base, run, analysis )
    if not os.path.exists(rootfile):
        log("root file %s doesnt exist" % rootfile )
        effhistos[key]=None
        return None
    f=None
    if openFiles.has_key ( rootfile ):
        f=openFiles[rootfile]
    else:
        f=ROOT.TFile(rootfile)
        if not f or not f.IsOpen():
            log("root file %s cannot be opened" % rootfile )
            effhistos[key]=None
            return None
        openFiles[rootfile]=f
    histo=f.Get("efficiency_%s" % topo )
    if not histo:
        log("histogram %s not found in %s" % ( topo, rootfile ))
        effhistos[key]=None
        return None
    effhistos[key]=histo
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

infoFields={}
def getMetaInfoField(analysis, field, run=None):
    """get one specific entry of the meta info
    """
    key = analysis + field + str(run)
    if infoFields.has_key(key):
        return infoFields[key]
    run = getRun(analysis, run)
    metainfo = parseMetaInfo(analysis, run)
    if not metainfo.has_key(field):
        infoFields[key]=None
        raise MetaInfoError("field ``%s'' not found for ``%s''" % (field,analysis) )
    f=metainfo[field]
    if len(f) == 0: 
        infoFields[key] = f
        return f
    while f[0] == ' ':
        f=f[1:]
    if f[-1] == '\n':
        f = f[:-1]
    infoFields[key] = f
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
     "T6bbWW": "#tilde{t}", 
     "TChiwz": "#tilde{#chi}^{0}_{2}#tilde{#chi}^{#pm}_{1}",
     #"TChiwz": "#tilde{#chi}^{0}",
     "TChizz": "#tilde{#chi}^{0}_{2}",
     "TChiSlep": "#tilde{#chi}^{0}_{2}",
     "TChiSlepSlep": "#tilde{#chi}^{0}_{2}#tilde{#chi}^{#pm}_{1}",
     "TChiNuSlep": "#tilde{#chi}^{0}_{2}",
}

def hasDictionary(analysis, run=None):
    """are the upper limits available in dictionary format?
    """
    if not run:
        run=getRun(analysis)
    dictfile="%s/%s/%s/sms.py" % ( Base, run, analysis )
    if os.path.exists(dictfile):
        return True
    log("Dictionary file %s doesnt exist" % dictfile )
    return False
#    hd = getMetaInfoField(analysis, "dictionary", run)
#    if hd:
#        if hd in [ '1', 'True' ]: return True
#    return False

def hasHistogram(analysis, run=None):
    """are the upper limits available as ROOT histogram?
    """
    if not run:
        run=getRun(analysis)
    rootfile="%s/%s/%s/sms.root" % ( Base, run, analysis )
    if os.path.exists(rootfile):
        return True
    log("ROOT file %s doesnt exist" % rootfile )
    return False

upperLimitDict={}
expupperLimitDict={}
def getUpperLimitDictionary ( analysis, topo, run, expected=False ):
    key=analysis+topo+str(run)
    if expected:
        if expupperLimitDict.has_key ( key ): return expupperLimitDict[key]
    else:
        if upperLimitDict.has_key ( key ): return upperLimitDict[key]
    dictfile="%s/%s/%s/sms.py" % ( Base, run, analysis )
    if not os.path.exists(dictfile):
        log("in getUpperLimitDictionary, dictionary file %s doesnt exist" % dictfile )
        upperLimitDict[key]=None
        return None
    Globals={}
    execfile(dictfile,Globals)
    if expected:
        Dict=Globals["ExpectedDict"]
        if not Dict.has_key ( topo ):
            log("dictionary doesnt have topology"+topo, "warning" )
            expupperLimitDict[key]=None
            return None
        expupperLimitDict[key]=Dict[topo]
    else:
        Dict=Globals["Dict"]
        if not Dict.has_key ( topo ):
            log("dictionary doesnt have topology"+topo, "warning" )
            upperLimitDict[key]=None
            return None
        upperLimitDict[key]=Dict[topo]
    return Dict[topo]

