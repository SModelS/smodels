#!/usr/bin/env python

"""
.. module:: smsHelpers
   :synopsis: Contains helper functions that are not to be used by the end
   user.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Doris Proschofsky <Doris.Proschofsky@assoc.oeaw.ac.at>
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>

"""
from os.path import os
import logging

logger = logging.getLogger(__name__)

Base = "/afs/hephy.at/user/w/walten/public/sms/"
runs = ["8TeV", "2012", "ATLAS8TeV", "2011", "RPV8", "RPV7"]
runs_ = {}
mlines = {}
infoFields = {}
pMI_ = {}
upperLimitDict = {}
expupperLimitDict = {}


def getRun(analysis, run=None):
    """
    Search for an analysis. If a specific run is given, then just check if
    results are present

    :returns: The run or None, if not found.
    
    """
    key=analysis+str(run)
    if key in runs_: return runs_[key]
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


def constraints(analysis, run):
    """
    Get all the constraints for an analysis-run pair.
    
    """
    return getLines(analysis, run, "constraint")


def getLines(analysis, run, label="condition"):
    """
    Get all <label> lines in info.txt for an analysis-run pair.
    
    """
    key = analysis + run + label
    if key in mlines:
        return mlines[key]
    info = "%s/%s/%s/info.txt" % (Base, run, analysis)
    ret = {}
    if not os.path.exists(info):
        logger.warn("Cannot find %s." % info)
        mlines[key] = ret
        return ret
    f = open(info)
    lines = f.readlines()
    f.close()
    for line in lines:
        if line.find("#") > -1:
            line = line[:line.find("#")].strip()
        if line == "":
            continue
        tokens = line.split(":",1)
        if not len(tokens) == 2:
            logger.warn("Cannot parse line: %s in %s" % (line, info))
            continue
        if tokens[0] != label:
            # Ignore everything but conditions
            continue
        excl = tokens[1]
        while excl[0] == " ":
            excl = excl[1:]
        if excl[-1] == '\n':
            excl = excl[:-1]
        keyvalue = excl.split(" -> ")
        if len(keyvalue) != 2:
            logger.warn("Cannot parse line: %s" % keyvalue)
        ret[keyvalue[0]] = keyvalue[1]
    mlines[key] = ret
    return ret


def getMetaInfoField(analysis, field, run=None):
    """
    Get one specific entry of the meta info.
    
    """
    key = analysis + field + str(run)
    if key in infoFields:
        return infoFields[key]
    run = getRun(analysis, run)
    metainfo = parseMetaInfo(analysis, run)
    if not field in metainfo:
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


def parseMetaInfo ( analysis, run ):
    """
    Get all the meta information for a given analysis-run pair.
    
    """
    key = analysis + str(run)
    if key in pMI_:
        return pMI_[key]
    info = "%s/%s/%s/info.txt" % (Base, run, analysis)
    ret = {}
    if not os.path.exists(info):
        logger.warn("cannot find %s" % info)
        pMI_[key] = ret
        return ret
    f = open(info)
    lines = f.readlines()
    f.close()
    for line in lines:
        if line.find("#") > -1:
            line = line[:line.find("#")].strip()
        line = line.replace("\n","")
        line = line.strip()
        if line == "":
            continue
        tokens = line.split(":",1)
        if not len(tokens) == 2:
            logger.error("Cannot parse line: ``%s'' in ``%s''" % (line, info))
            continue
        if tokens[0] == "exclusions":
            # Ignore exclusions
            continue
        ret[tokens[0]] = tokens[1]
    pMI_[key] = ret
    return ret


def fuzzyconditions(analysis, run):
    """
    Get all fuzzyconditions for an analysis-run pair.
    
    """
    return getLines(analysis, run, "fuzzycondition")


def conditions(analysis, run):
    """
    Get all conditions for an analysis-run pair.
    
    """
    return getLines(analysis, run, "condition")


def getUpperLimitDictionary(analysis, topo, run, expected=False):
    """
    TODO: write docstring
    
    """
    key = analysis + topo + str(run)
    if expected:
        if key in expupperLimitDict:
            return expupperLimitDict[key]
    else:
        if key in upperLimitDict:
            return upperLimitDict[key]
    dictfile="%s/%s/%s/sms.py" % (Base, run, analysis)
    if not os.path.exists(dictfile):
        logger.warn("Dictionary file %s does not exist." % dictfile)
        upperLimitDict[key] = None
        return None
    Locals = {}
    execfile(dictfile, Locals)
    if expected:
        if not Locals.has_key("ExpectedDict"):
            logger.warn("Expected dictionary is missing for analysis " + analysis)
            return None
        Dict = Locals["ExpectedDict"]
        if not Dict.has_key(topo):
            logger.warn("Dictionary does not have topology " + topo )
            expupperLimitDict[key] = None
            return None
        expupperLimitDict[key]=Dict[topo]
    else:
        if not Locals.has_key("Dict"):
            logger.warn("Observed dictionary is missing for analysis " + analysis)
            return None
        Dict=Locals["Dict"]
        if not Dict.has_key(topo):
            logger.warn("Dictionary does not have topology " + topo )
            upperLimitDict[key] = None
            return None
        upperLimitDict[key] = Dict[topo]
    return Dict[topo]


def hasDictionary(analysis, run=None, topo=None):
    """
    Check if available upper limits are in dictionary format.
    
    """
    if not run:
        run = getRun(analysis)
    dictfile = "%s/%s/%s/sms.py" % (Base, run, analysis)
    if os.path.exists(dictfile):
        if topo == None:
            return True
        Locals = {}
        execfile(dictfile, Locals)
        Dict = Locals["Dict"]
        if Dict.has_key(topo):
            return True
        return False

    logger.warn("Dictionary file %s does not exist." % dictfile )
    return False
