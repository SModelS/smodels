#!/usr/bin/env python

"""
.. module:: experiment.smsHelpers
   :synopsis: Contains private helper functions to access the SMS results.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Doris Proschofsky <Doris.Proschofsky@assoc.oeaw.ac.at>
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>

"""

from os.path import os
import logging

logger = logging.getLogger(__name__)

base = "/database/path/not/set"
mlines = {}
infoFields = {}
pMI_ = {}
upperLimitDict = {}
expupperLimitDict = {}

sqrts = ["8TeV"]
experiments = ["ATLAS", "CMS"]
paths_ = {}

def getPath(analysis, path=None):
    """
    Search for an analysis
    """
    if analysis in paths_:
        return paths_[analysis]
    if path:
        if os.path.exists("%s/%s/%s" % (base, path, analysis)):
            paths_[analysis] = str(path)
            return paths_[analysis]
    for s in sqrts:
        for e in experiments:
            if os.path.exists("%s/%s/%s/%s" % (base, s, e, analysis)):
                paths_[analysis] = "%s/%s" % (s, e)
                return paths_[analysis]
    paths_[analysis] = None
    return None


def getLines(analysis, path, label="condition"):
    """
    Get all <label> lines in info.txt for an analysis-path pair.

    """
    key = analysis + path + label
    if key in mlines:
        return mlines[key]
    info = "%s/%s/%s/info.txt" % (base, path, analysis)
    ret = {}
    if not os.path.exists(info):
        logger.warn("Cannot find %s.", info)
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
        tokens = line.split(":", 1)
        if not len(tokens) == 2:
            logger.warn("Cannot parse line: %s in %s", line, info)
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
            logger.warn("Cannot parse line: %s", keyvalue)
        ret[keyvalue[0]] = keyvalue[1]
    mlines[key] = ret
    return ret


def getMetaInfoField(analysis, field, path):
    """
    Get one specific entry of the meta info.

    """
    key = analysis + field + str(path)
    if key in infoFields:
        return infoFields[key]
    path = getPath(analysis, path)
    metainfo = _parseMetaInfo(analysis, path)
    if not field in metainfo:
        infoFields[key] = None
        return infoFields[key]
    f = metainfo[field]
    if len(f) == 0:
        infoFields[key] = f
        return f
    while f[0] == ' ':
        f = f[1:]
    if f[-1] == '\n':
        f = f[:-1]
    infoFields[key] = f
    return f


def _parseMetaInfo(analysis, path):
    """
    Get all the meta information for a given analysis-path pair.

    """
    key = analysis + str(path)
    if key in pMI_:
        return pMI_[key]
    info = "%s/%s/%s/info.txt" % (base, path, analysis)
    ret = {}
    if not os.path.exists(info):
        logger.warn("Cannot find %s.", info)
        pMI_[key] = ret
        return ret
    f = open(info)
    lines = f.readlines()
    f.close()
    for line in lines:
        if line.find("#") > -1:
            line = line[:line.find("#")].strip()
        line = line.replace("\n", "")
        line = line.strip()
        if line == "":
            continue
        tokens = line.split(":", 1)
        if not len(tokens) == 2:
            logger.error("Cannot parse line: '%s' in '%s'", line, info)
            continue
        if tokens[0] == "exclusions":
            # Ignore exclusions
            continue
        ret[tokens[0]] = tokens[1]
    pMI_[key] = ret
    return ret


def getUpperLimitDictionary(analysis, topology, path, expected=False):
    """
    Returns a dictionary containing the raw Upper Limit data for the analysis
    and topology. 

    """
    key = analysis + topology + str(path)
    if expected:
        if key in expupperLimitDict:
            return expupperLimitDict[key]
    else:
        if key in upperLimitDict:
            return upperLimitDict[key]
    dictfile = "%s/%s/%s/sms.py" % (base, path, analysis)
    if not os.path.exists(dictfile):
        logger.warn("Dictionary file %s does not exist.", dictfile)
        upperLimitDict[key] = None
        return None
    localsDictionary = {}
    execfile(dictfile, localsDictionary)
    if expected:
        if not localsDictionary.has_key("ExpectedDict"):
            logger.warn("Expected dictionary is missing for analysis %s.",
                        analysis)
            return None
        dictionary = localsDictionary["ExpectedDict"]
        if not dictionary.has_key(topology):
            logger.warn("Dictionary does not have topology %s.", topology)
            expupperLimitDict[key] = None
            return None
        expupperLimitDict[key] = dictionary[topology]
    else:
        if not localsDictionary.has_key("Dict"):
            logger.warn("Observed dictionary is missing for analysis %s.",
                        analysis)
            return None
        dictionary = localsDictionary["Dict"]
        if not dictionary.has_key(topology):
            logger.warn("Dictionary does not have topology %s.", topology)
            upperLimitDict[key] = None
            return None
        upperLimitDict[key] = dictionary[topology]
    return dictionary[topology]


def hasDictionary(analysis, path=None, topology=None):
    """
    Check if available upper limits are in dictionary format.

    """
    if not path:
        path = getPath(analysis)
    dictfile = "%s/%s/%s/sms.py" % (base, path, analysis)
    if os.path.exists(dictfile):
        if topology == None:
            return True
        localsDictionary = {}
        execfile(dictfile, localsDictionary)
        dictionary = locals["Dict"]
        if dictionary.has_key(topology):
            return True
        return False

    logger.warn("Dictionary file %s does not exist.", dictfile)
    return False

def databaseVersion(astuple=False):
    """ prints out version number of the *database* """
    f = open("%s/version" % base)
    l = f.readline()
    f.close()
    l = l.replace("\n", "")
    l.strip()
    if not astuple:
        return l
    A, B = l.split(".")
    return (int(A), int(B))
