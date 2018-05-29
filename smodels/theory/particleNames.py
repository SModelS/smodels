#!/usr/bin/env python

"""
.. module:: particleNames
   :synopsis: Provides functions for getting particle objects from pdg ids or name, and
              other helpers.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
from smodels.particleDefinitions import SMList, BSMList, SMnames, particleLists, SM
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger



def getObjectFromPdg(pdg):
    """
    Convert pdg number to particle object according to the Particle class.
    :type pdg: int
    :returns: Particles object 
    """
    found = False
    for particle in SMList+BSMList:
        if particle.pdg==pdg: 
            found = True
            p = particle
    if found: return p
    else: logger.warning("Particle %i not defined in particleClass.py" %(pdg))
    
    
    
def getObjectFromName(name):
    """
    Convert particle label to object according to the Particle class.
    :type name: str
    :returns: Particles object or ParticleList object
    """    
    found = False
    for particle in SM + BSMList:
        if particle.label == name: 
            found = True
            p = particle

    if found: return p
    else: logger.warning("Particle %s not defined in particleClass.py" %(name))



def getNamesList(particleList):
	""" 
	Convert list of particles to list of particle names according to the Particle class.
	:type particleList: list of instances of Particles class
	:returns: list of str
	"""
	NamesList = [particle.label for particle in particleList]
	return NamesList


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
    if isinstance(instring,str):
        outstr = instring
    elif isinstance(instring,list):
        for st in instring:
            if type(st) != type('st'):
                logger.error("Input must be a string or a list of strings")
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
                #if ptc == '*':
                #    ptc = StrWildcard()                            
                if not ptc in SMnames and not ptc in getNamesList(particleLists):
                    logger.error("Unknown particle. Add " + ptc + " to smodels/particleDefinitions.py")
                    raise SModelSError()

    # Check if there are not unmatched ['s and/or ]'s in the string
    if nc != 0:
        raise SModelSError("Wrong input (incomplete elements?) " + instring)

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
                #if ptc == '*':
                #    ptc = StrWildcard()                    
                if not ptc in SMnames and not ptc in getNamesList(particleLists):
                    logger.error("Unknown particle. Add " + ptc + " to smodels/particleDefinitions.py")
                    raise SModelSError()
            vertStr = ""

    # Check if there are not unmatched ['s and/or ]'s in the string
    if nc != 0:
        logger.error("Wrong input (incomplete elements?) " + instring)
        raise SModelSError()

    return vertices




class StrWildcard(str):
    """
    A string wildcard class. It will return True when compared to any other string.
    """
    
    def __init__(self):
        str.__init__(self)
        
    def __str__(self):
        return '*'    

    def __repr__(self):
        return self.__str__()

    def __cmp__(self,other):
        if isinstance(other,str):
            return 0
        else:
            return -1

    def __eq__(self,other):
        return self.__cmp__(other) == 0  
    
    def __ne__(self,other):
        return self.__cmp__(other) != 0
     





