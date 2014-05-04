#!/usr/bin/env python

"""
.. module:: ParticleNames
   :synopsis: Provides functions for getting particle names from pdg ids, and
   other helpers.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import logging

logger = logging.getLogger(__name__)


def getName(pdg):
    """
    Convert pdg number to particle name according to the dictionaries rOdd and
    rEven.

    :type pdg: int
    :returns: particle name (e.g. gluino, mu-, ...)
    
    """
    p = int(pdg)
    if p in rOdd:
        return rOdd[p]
    if p in rEven:
        return rEven[p]
    else:
        return False


def getPdg(name):
    """
    Convert a name to the pdg number according to the dictionaries rOdd and
    rEven.

    :type name: string
    :returns: particle pdg; None, if name could not be resolved
    
    """
    for (pdg, pname) in rOdd.items():
        if name == pname:
            return abs(pdg)
    for (pdg, pname) in rEven.items():
        if name == pname:
            return abs(pdg)
    return None


def elementsInStr(instring):
    """
    Parse instring and return a list of elements appearing in instring.
    
    instring can also be a list of strings.
    
    :returns: list of elements appearing in instring in string format
    
    """
    if type(instring) == type('st'):
        outstr = instring
    elif type(instring) == type([]):
        outstr = ""
        for st in instring:
            if type(st) != type('st'):
                logger.error("Input must be a string or a list of strings")
                return False
            # Combine list of strings in a single string
            outstr += st

    elements = []
    outstr = outstr.replace(" ", "").replace("'", "")
    elStr = ""
    nc = 0
    # Parse the string and looks for matching ['s and ]'s, when the matching is
    # complete, store element
    for c in outstr:
        delta = 0
        if c == '[':
            delta = -1
        elif c == ']':
            delta = 1
        nc += delta
        if nc != 0:
            elStr += c
        if nc == 0 and delta != 0:
            elements.append(elStr + c)
            elStr = ""
            # Syntax checks
            ptclist = elements[-1].replace(']', ',').replace('[', ',').\
                    split(',')
            for ptc in ptclist:
                if not ptc:
                    continue
                if not ptc in rEven.values() and not ptc in ptcDic:
                    logger.error("Unknown particle " + ptc)
                    return False

    # Check if there are not unmatched ['s and/or ]'s in the string
    if nc != 0:
        logger.error("Wrong input (incomplete elements?) " + instring)
        return False

    return elements


def vertInStr(instring):
    """
    Parses instring (or a list of strings) and returns the list of particle
    vertices appearing in instring.
    
    """
    if type(instring) == type('st'):
        outstr = instring
    elif type(instring) == type([]):
        outstr = ""
        for st in instring:
            if type(st) != type('st'):
                logger.error("Input must be a string or a list of strings")
                return False
            # Combine list of strings in a single string
            outstr += st

    vertices = []
    outstr = outstr.replace(" ", "").replace("'", "")
    vertStr = ""
    nc = 0
    # Parse the string and looks for matching ['s and ]'s, when the matching is
    # complete, store element
    for c in outstr:
        delta = 0
        if c == '[':
            delta = -1
        elif c == ']':
            delta = 1
        nc += delta
        if c == '[':
            vertStr = ""
        if nc != 0 and c != '[' and c != ']':
            vertStr += c
        if delta > 0 and vertStr:
            vertices.append(vertStr.split(','))
            # Syntax checks:
            for ptc in vertices[-1]:
                if not ptc:
                    continue
                if not ptc in rEven.values() and not ptc in ptcDic:
                    logger.error("Unknown particle " + ptc)
                    return False
            vertStr = ""

    # Check if there are not unmatched ['s and/or ]'s in the string
    if nc != 0:
        logger.error("Wrong input (incomplete elements?) " + instring)
        return False

    return vertices


def simParticles(ptype1, ptype2, useDict=True):
    """
    Compares 2 particle names or 2 nested name arrays. Allows for dictionary
    labels (Ex: L = l, l+ = l, l = l-,...). For the last nested level ignore
    particle ordering. FIXME nesting? 
 
    :param ptype1: first (nested) list of particle names, e.g. ['l','jet']
    :param ptype2: second (nested) list of particle names 
    :param useDict: use the translation dictionary, i.e. allow e to stand for
    e+ or e-, l+ to stand for e+ or mu+, etc 
    :returns: boolean
    
    """
    import copy

    wrongFormat = False
    if type(ptype1) != type(ptype2):
        return False
    # Check for nested arrays (should be in the standard notation [[[]],[[]]])
    if type(ptype1) == type([]):
        if len(ptype1) != len(ptype2):
            return False
        for ib, br in enumerate(ptype1):
            if type(br) != type(ptype2[ib]) or type(br) != type([]):
                wrongFormat = True
            # Check number of vertices in branch
            if len(ptype1[ib]) != len(ptype2[ib]):
                return False
            for iv, vt in enumerate(br):
                if type(vt) != type(ptype2[ib][iv]) or type(vt) != type([]):
                    wrongFormat = True
                # Check number of particles in vertex
                if len(ptype1[ib][iv]) != len(ptype2[ib][iv]):
                    return False
                for ptc in ptype1[ib][iv] + ptype2[ib][iv]:
                    if not ptc in ptcDic.keys() + rEven.values():
                        wrongFormat = True
    if wrongFormat:
        logger.error("Wrong input format!" + str(ptype1) + " " + str(ptype2))
        return False

    # Put input in standard notation
    if type(ptype1) == type("str"):
        ptype1v = [[[ptype1]], [[ptype1]]]
        ptype2v = [[[ptype2]], [[ptype2]]]
    else:
        ptype1v = ptype1[:]
        ptype2v = ptype2[:]

    # Loop over branches
    for ibr, br in enumerate(ptype1v):
        # Loop over vertices
        for iv, vt in enumerate(br):
            # Check if lists match, ignoring possible dictionary entries
            pmatch = True
            for ptc in ptype1v[ibr][iv]:
                if ptype1v[ibr][iv].count(ptc) != ptype2v[ibr][iv].count(ptc):
                    pmatch = False
            if pmatch:
                continue
            elif not useDict:
                return False
            # If they do not match and useDict=True, generate all possible
            # lists from dictionary entries:
            allptcs = [[ptype1v[ibr][iv]], [ptype2v[ibr][iv]]]
            for allpt in allptcs:
                ptc0 = copy.deepcopy(allpt[0])
                for ipt in range(len(ptc0)):
                    if ptc0[ipt] in ptcDic:
                        for jpt in range(len(allpt)):
                            if allpt[jpt] == []:
                                continue
                            newptc = copy.deepcopy(allpt[jpt])
                            for ptc in ptcDic[ptc0[ipt]]:
                                newptc[ipt] = ptc
                                allpt.append(copy.deepcopy(newptc))
                            allpt[jpt] = []
                while allpt.count([]) > 0:
                    allpt.remove([])

            # Compare all possibilities
            match = False
            iA = 0
            while not match and iA < len(allptcs[0]):
                ptcA = allptcs[0][iA]
                for ptcB in allptcs[1]:
                    if len(ptcA) != len(ptcB):
                        return False
                    pmatch = True
                    for ptc in ptcA:
                        if ptcA.count(ptc) != ptcB.count(ptc):
                            pmatch = False
                    if pmatch:
                        match = True
                        break
                iA += 1
            if not match:
                return False
    # Entries are similar
    return True


rOdd = {1000021 : "gluino",
        1000022 : "N1",
        1000023 : "N2",
        1000025 : "N3",
        1000035 : "N4",
        1000024 : "C1",
        1000037 : "C2",
        1000039 : "gravitino",
        1000001 : "squark",
        1000002 : "squark",
        1000003 : "squark",
        1000004 : "squark",
        2000001 : "squark",
        2000002 : "squark",
        2000003 : "squark",
        2000004 : "squark",
        1000005 : "sbottom",
        2000005 : "sbottom",
        1000006 : "stop",
        2000006 : "stop",
        1000011 : "slepton",
        1000013 : "slepton",
        1000015 : "stau",
        2000011 : "slepton",
        2000013 : "slepton",
        2000015 : "stau",
        1000012 : "sneutrino",
        1000014 : "sneutrino",
        1000016 : "sneutrino",
        2000012 : "sneutrino",
        2000014 : "sneutrino",
        2000016 : "sneutrino",
        - 1000021 : "gluino",
        - 1000022 : "N1",
        - 1000023 : "N2",
        - 1000025 : "N3",
        - 1000035 : "N4",
        - 1000024 : "C1",
        - 1000037 : "C2",
        - 1000039 : "gravitino",
        - 1000001 : "squark",
        - 1000002 : "squark",
        - 1000003 : "squark",
        - 1000004 : "squark",
        - 2000001 : "squark",
        - 2000002 : "squark",
        - 2000003 : "squark",
        - 2000004 : "squark",
        - 1000005 : "sbottom",
        - 2000005 : "sbottom",
        - 1000006 : "stop",
        - 2000006 : "stop",
        - 1000011 : "slepton",
        - 1000013 : "slepton",
        - 1000015 : "stau",
        - 2000011 : "slepton",
        - 2000013 : "slepton",
        - 2000015 : "stau",
        - 1000012 : "sneutrino",
        - 1000014 : "sneutrino",
        - 1000016 : "sneutrino",
        - 2000012 : "sneutrino",
        - 2000014 : "sneutrino",
        - 2000016 : "sneutrino"}

rEven = {25 : "higgs",
         - 25: "higgs",
         35 : "H0",
         - 35: "H0",
         36 : "A0",
         - 36: "A0",
         37 : "H+",
         - 37: "H-",
         23 : "Z",
         - 23: "Z",
         22 : "photon",
         - 22: "photon",
         24 : "W+",
         - 24: "W-",
         16 : "nu",
         - 16: "nu",
         15 : "ta-",
         - 15: "ta+",
         14 : "nu",
         - 14: "nu",
         13 : "mu-",
         - 13: "mu+",
         12 : "nu",
         - 12: "nu",
         11 : "e-",
         - 11: "e+",
         5  : "b",
         - 5 : "b",
         6  : "t+",
         - 6 : "t-",
         1  : "jet",
         2  : "jet",
         3  : "jet",
         4  : "jet",
         21 : "jet",
         - 1 : "jet",
         - 2 : "jet",
         - 3 : "jet",
         - 4 : "jet",
         - 21: "jet"}

ptcDic = {"e" : ["e+", "e-"],
          "mu" : ["mu+", "mu-"],
          "ta" : ["ta+", "ta-"],
          "l+" : ["e+", "mu+"],
          "l-" : ["e-", "mu-"],
          "l" : ["e-", "mu-", "e+", "mu+"],
          "W" : ["W+", "W-"],
          "t" : ["t+", "t-"],
          "L+" : ["e+", "mu+", "ta+"],
          "L-" : ["e-", "mu-", "ta-"],
          "L" : ["e+", "mu+", "ta+", "e-", "mu-", "ta-"]}
