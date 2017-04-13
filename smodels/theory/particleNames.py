#!/usr/bin/env python

"""
.. module:: particleNames
   :synopsis: Provides functions for getting particle names from pdg ids, and
              other helpers.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.particles import rEven, rOdd, ptcDic, qNumbers, finalStates
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
import itertools

from smodels.tools.smodelsLogging import logger


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


def elementsInStr(instring,removeQuotes=True):
    """
    Parse instring and return a list of elements appearing in instring.
    instring can also be a list of strings.
    
    :param instring: string containing elements (e.g. "[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]")
    :param removeQuotes: If True, it will remove the quotes from the particle labels.
                         Set to False, if one wants to run eval on the output.
    
    :returns: list of elements appearing in instring in string format
    
    """
    
    outstr = ""
    if type(instring) == type('st'):
        outstr = instring
    elif type(instring) == type([]):
        for st in instring:
            if type(st) != type('st'):
                logger.error("Input must be a string or a list of strings. Type %s found:\n %s" %(type(st),str(st)))
                raise SModelSError()
            # Combine list of strings in a single string
            outstr += st
    else:
        raise SModelSError ( "syntax error in constraint/condition: ``%s''." \
              "Check your constraints and conditions in your database." % str(instring) )

    elements = []
    outstr = outstr.replace(" ", "")
    if removeQuotes:
        outstr = outstr.replace("'", "")
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
                ptc = ptc.replace("'","")
                if not ptc:
                    continue
                if not ptc in rEven.values() and not ptc in ptcDic and ptc != '*':
                    logger.error("Unknown particle. Add " + ptc + " to smodels/particles.py")
                    raise SModelSError()

    # Check if there are not unmatched ['s and/or ]'s in the string
    if nc != 0:
        logger.error("Wrong input (incomplete elements?) " + instring)
        raise SModelSError()

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
                raise SModelSError()
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
                if not ptc in rEven.values() and not ptc in ptcDic and ptc != '*':
                    logger.error("Unknown particle. Add " + ptc + " to smodels/particle.py")
                    raise SModelSError()
            vertStr = ""

    # Check if there are not unmatched ['s and/or ]'s in the string
    if nc != 0:
        logger.error("Wrong input (incomplete elements?) " + instring)
        raise SModelSError()

    return vertices

def simParticles(plist1, plist2, useDict=True):
    """
    Compares two lists of particle names. Allows for dictionary
    labels (Ex: L = l, l+ = l, l = l-,...). Ignores particle ordering inside
    the list
 
    :param plist1: first list of particle names, e.g. ['l','jet']
    :param plist2: second list of particle names 
    :param useDict: use the translation dictionary, i.e. allow e to stand for
                    e+ or e-, l+ to stand for e+ or mu+, etc 
    :returns: True/False if the particles list match (ignoring order)    
    """

    if not isinstance(plist1,list) or type(plist1) != type(plist2):
        logger.error("Input must be a list")
        raise SModelSError()
    if len(plist1) != len(plist2):
        return False
    for i,p in enumerate(plist1):
        if not isinstance(p,str) or not isinstance(plist2[i],str):
            logger.error("Input must be a list of particle strings")
            raise SModelSError()
        elif not p in list ( ptcDic.keys() ) + list ( rEven.values() ):
            logger.error("Unknow particle: %s" %p)
            raise SModelSError()
        elif not plist2[i] in list ( ptcDic.keys() ) + list ( rEven.values() ):
            logger.error("Unknow particle: %s" %plist2[i])
            raise SModelSError()
                        
        
    l1 = sorted(plist1)
    l2 = sorted(plist2)
    if not useDict:
        return l1 == l2
    
    #If dictionary is to be used, replace particles by their dictionay entries
    #e.g. [jet,mu+] -> [[q,g,c],[mu+]], [jet,mu] -> [[q,g,c],[mu+,mu-]] 
    extendedL1 = []    
    for i,p in enumerate(plist1):
        if not p in ptcDic:
            extendedL1.append([p])
        else:
            extendedL1.append(ptcDic[p])
    extendedL2 = []    
    for i,p in enumerate(plist2):
        if not p in ptcDic:
            extendedL2.append([p])
        else:
            extendedL2.append(ptcDic[p])
    
    #Generate all combinations of particle lists (already sorted to avoid ordering issues)
    #e.g. [[q,g,c],[mu+]] -> [[q,mu+],[g,mu+],[c,mu+]]
    extendedL1 = [sorted(list(i)) for i in itertools.product(*extendedL1)]
    extendedL2 = [sorted(list(i)) for i in itertools.product(*extendedL2)]

    #Now compare the two lists and see if there is a match:
    for plist in extendedL1:
        if plist in extendedL2: return True
        
    return False

def getFinalStateLabel(pid):
    """
    Given the particle PID, returns the label corresponding to its final state
    (e.g. 1000022 -> MET, 1000023 -> HSCP,...)
    :parameter pid: PDG code for particle (must appear in particles.py)
    :return: Final state string (e.g. MET, HSCP,...)
    """

    if not abs(pid) in qNumbers:
        logger.error("qNumbers are not defined for %i. Please, add it to particles.py." %pid)
        raise SModelSError
    elif not pid in qNumbers:  #Use the anti-particle info:
        pidQnumber = qNumbers[abs(pid)]
        pidQnumber[1] = -pidQnumber[1] #Flip the charge sign
    else:    
        pidQnumber = qNumbers[pid]
    for key,qnumberList in finalStates.items():
        if pidQnumber in qnumberList:
            return key
    
    logger.error("Final state for %i not found. Please, add it to particles.py." %pid)
    raise SModelSError

