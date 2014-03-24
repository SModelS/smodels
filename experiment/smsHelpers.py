#!/usr/bin/env python

"""
.. module:: smsHelpers
        :synopsis: Some helper functions that are not to be used by the end user.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Doris Proschofsky <Doris.Proschofsky@assoc.oeaw.ac.at>

"""

import os
from tools.physicsUnits import rmvunit
from experiment.experimentExceptions import MetaInfoError
import logging                                                                  
                                                                                
logger = logging.getLogger(__name__)

Base = "/afs/hephy.at/user/w/walten/public/sms/"

runs=[ "8TeV", "2012", "ATLAS8TeV", "2011", "RPV8", "RPV7" ]
## runs=[ "2012" ]

def close():
    """ close all open files """
    for (name,tfile) in openFiles.items():
        tfile.Close()
    openFiles={}

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
    for trun in runs:
        if os.path.exists ( "%s/%s/%s" % ( Base, trun, analysis ) ):
            runs_[key]=trun
            return trun
    runs_[key]=None
    return None

pMI_={}
def parseMetaInfo ( analysis, run ):
    """ get all the meta information for a given analysis/run pair """
    key=analysis+str(run)
    if pMI_.has_key ( key ): return pMI_[key]
    info="%s/%s/%s/info.txt" % ( Base, run, analysis )
    ret={}
    if not os.path.exists ( info ):
        logger.warn ( "cannot find %s" % info )
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
            logger.error ( "[117] cannot parse this line (1): ``%s'' in ``%s''" % (line, info) )
            continue
        if tokens[0]=="exclusions":
    # we treat these separately
            continue
        ret[tokens[0]]=tokens[1]
    pMI_[key]=ret
    return ret

mlines={}
def getLines ( analysis, run, label="condition" ):
    """ get all <label> lines in info.txt for an analysis/run pair 
    """
    key=analysis+run+label
    if mlines.has_key ( key ): return mlines[key]
    info="%s/%s/%s/info.txt" % ( Base, run, analysis )
    ret={}
    if not os.path.exists ( info ):
        logger.warn ("cannot find %s" % info )
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
            logger.warn ( "[172] cannot parse this line (2): %s in %s" % (line, info) )
            continue
        if tokens[0]!=label:
    # we're only interested in the conditions
            continue
        excl=tokens[1]
        while excl[0]==" ": excl=excl[1:]
        if excl[-1]=='\n': excl=excl[:-1]
        keyvalue=excl.split(" -> ")
        if len(keyvalue)!=2:
            logger.warn ( "[185] cannot parse the following line: %s" % keyvalue )
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
        return infoFields[key]
# raise MetaInfoError("field ``%s'' not found for ``%s''" % (field,analysis) )
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

def getErrorMessage ( Dict, mx, my ):
    logger.error("error message not implemented.")
    return ""

def hasDictionary(analysis, run=None):
    """are the upper limits available in dictionary format?
    """
    if not run:
        run=getRun(analysis)
    dictfile="%s/%s/%s/sms.py" % ( Base, run, analysis )
    if os.path.exists(dictfile):
        return True
    logger.warn("Dictionary file %s doesnt exist" % dictfile )
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
        logger.warn("in getUpperLimitDictionary, dictionary file %s doesnt exist" % dictfile )
        upperLimitDict[key]=None
        return None
    Globals={}
    execfile(dictfile,Globals)
    if expected:
        Dict=Globals["ExpectedDict"]
        if not Dict.has_key ( topo ):
            logger.warn("dictionary doesnt have topology "+topo )
            expupperLimitDict[key]=None
            return None
        expupperLimitDict[key]=Dict[topo]
    else:
        Dict=Globals["Dict"]
        if not Dict.has_key ( topo ):
            logger.warn("dictionary doesnt have topology "+topo )
            upperLimitDict[key]=None
            return None
        upperLimitDict[key]=Dict[topo]
    return Dict[topo]

def databaseVersion( astuple=False, addCodeName=True ):
    """ prints out version number of the *database* """
    f=open("%s/version" % Base )
    l=f.readline()
    f.close()
    l=l.replace("\n","")
    l.strip()
    if not astuple:
        if addCodeName: return l
        p=l.find("/")
        if p>-1: return l[:p]
    T,C=l.split("/")
    A,B=T.split(".")
    if addCodeName:
        return (int(A),int(B),C.strip())
    else:
        return (int(A),int(B))
