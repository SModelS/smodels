#!/usr/bin/env python

"""
.. module:: particleNames
   :synopsis: Provides functions for getting particle names from pdg ids, and
              other helpers.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
import copy
#from smodels.particles import rEven, rOdd, ptcDic
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
    
    from smodels.particles import rEven, rOdd
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
    from smodels.particles import rEven, rOdd
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
    
    elements = []
    if isinstance(instring,list):
        for elStr in instring:
            elements += elementsInStr(elStr, removeQuotes)
    elif isinstance(instring,str):
        if instring.count('[') != instring.count(']'):
            raise SModelSError("Syntax error in string: ``%s''." \
                  "Check your constraints and conditions in your database." % str(instring))

        tempstring = instring.replace(" ", "")
        from smodels.particles import rEven, ptcDic
        while tempstring.find('[') != -1:
            el0 = tempstring.find('[')
            elf = el0
            brackDiff = tempstring.count('[',el0,elf+1)-tempstring.count(']',el0,elf+1)
            while brackDiff != 0:
                elf = tempstring.find(']',elf+1)
                brackDiff = tempstring.count('[',el0,elf+1)-tempstring.count(']',el0,elf+1)
            tempEl = tempstring[el0:elf+1]
            tempEl = tempEl.replace(']'," ").replace("["," ").replace(","," ")
            particles = list(set([str(ptc) for ptc in tempEl.split() if ptc]))
            particles = [ptc.replace("'","").replace('"',"").strip() for ptc in particles]
            unknownParticles = [ptc for ptc in particles if ((not ptc in list(rEven.values()))
                                                         and (not ptc in ptcDic) and ptc != '*')]
            if unknownParticles:
                raise SModelSError("Unknown particles: %s. Add missing particles to smodels/particles.py" 
                                   %str(unknownParticles))
            if removeQuotes:
                elements.append(tempstring[el0:elf+1].replace("'","").replace('"',""))
            else:
                elements.append(tempstring[el0:elf+1])
            tempstring = tempstring[:el0]+tempstring[elf+1:]
        return elements
    else:
        raise SModelSError("Input must be a string or a list of strings. Type %s found:\n %s" 
                           %(type(instring),instring))



def vertInStr(instring,removeQuotes=True):
    """
    Parses instring (or a list of strings) and returns the list of particle
    vertices appearing in instring.
    
    :param instring: string containing elements (e.g. "[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]")
    :param removeQuotes: If True, it will remove the quotes from the particle labels.
                         Set to False, if one wants to run eval on the output.
    
    :returns: list of elements appearing in instring in string format
    
    """
    
    vertices = []
    if isinstance(instring,list):
        for elStr in instring:
            vertices += elementsInStr(elStr, removeQuotes)
    elif isinstance(instring,str):
        if instring.count('[') != instring.count(']'):
            raise SModelSError("Syntax error in string: ``%s''." \
                  "Check your constraints and conditions in your database." % str(instring))

        tempstring = instring.replace(" ", "")
        from smodels.particles import rEven, ptcDic
        while tempstring.find(']') != -1:
            elf = tempstring.find(']')
            el0 = tempstring.rfind('[',0,elf)
            tempEl = tempstring[el0:elf+1]
            tempEl = tempEl.replace(']'," ").replace("["," ").replace(","," ")
            if not tempEl.strip():
                tempstring = tempstring[:el0]+tempstring[elf+1:]
                continue
            particles = list(set([str(ptc) for ptc in tempEl.split() if ptc]))
            particles = [ptc.replace("'","").replace('"',"").strip() for ptc in particles]
            unknownParticles = [ptc for ptc in particles if ((not ptc in list(rEven.values()))
                                                         and (not ptc in ptcDic) and ptc != '*')]
            if unknownParticles:
                raise SModelSError("Unknown particles: %s. Add missing particles to smodels/particles.py" 
                                   %str(unknownParticles))
            if removeQuotes:
                vertices.append(tempstring[el0:elf+1].replace("'","").replace('"',""))
            else:
                vertices.append(tempstring[el0:elf+1])
            tempstring = tempstring[:el0]+tempstring[elf+1:]
        return vertices
    else:
        raise SModelSError("Input must be a string or a list of strings. Type %s found:\n %s" 
                           %(type(instring),instring))
        
        
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
    from smodels.particles import rEven, ptcDic
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

