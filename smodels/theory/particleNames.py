#!/usr/bin/env python

"""
.. module:: particleNames
   :synopsis: Provides functions for getting particle objects from pdg ids or name, and
              other helpers.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
import copy
from smodels.particleDefinitions import SMList, BSMList, SMparticles, particleLists, SM
from smodels.theory.particleClass import particleInList
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
import itertools
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
    if type(instring) == type('st'):
        outstr = instring
    elif type(instring) == type([]):
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
                if not ptc in SMparticles and not ptc in getNamesList(particleLists):
                    logger.error("Unknown particle. Add " + ptc + " to smodels/particleDefinitions.py")
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
                if not ptc in SMparticles and not ptc in getNamesList(particleLists):
                    logger.error("Unknown particle. Add " + ptc + " to smodels/particleDefinitions.py")
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
        #if hasattr(plist2[i], 'label'):
        if not isinstance(p.label,str) or not isinstance(plist2[i].label,str) :
            logger.error("Input must be two lists of particle or particle list objects")
            raise SModelSError()
        """    
        if type(plist2[i]) == str:            
            if not isinstance(p.label,str) or not isinstance(plist2[i],str) :
                logger.error("Input must be two lists of str")
                raise SModelSError()
        """    
    l1 = sorted(getNamesList(plist1))       

    if type(plist2[0])==str: l2 = sorted(plist2)
    else: l2 = sorted(getNamesList(plist2))    
    
    if not useDict:
        return l1 == l2
        
    #If dictionary is to be used, replace particles by their dictionay entries
    #e.g. [jet,mu+] -> [[q,g,c],[mu+]], [jet,mu] -> [[q,g,c],[mu+,mu-]] 

    extendedL1 = []    
    for i,p in enumerate(plist1):
        if not particleInList(p,particleLists, inlist=False): 
            extendedL1.append([p.label])
        else:
            for pList in particleLists:

                if pList.label == p.label: extendedL1.append( pList.getNames() )
                
    extendedL2 = []    
    for i,p in enumerate(plist2):
        if not particleInList(p,particleLists, inlist=False):
            if type(p)== str: extendedL2.append([p])
            else: extendedL2.append([p.label])

        else:    
            for pList in particleLists:
                if ( type(p) == str and pList.label == p ) or ( type(p) != str and pList.label == p.label): 
                    extendedL2.append( pList.getNames() )           

    #Generate all combinations of particle lists (already sorted to avoid ordering issues)
    #e.g. [[q,g,c],[mu+]] -> [[q,mu+],[g,mu+],[c,mu+]]
    extendedL1 = [sorted(list(i)) for i in itertools.product(*extendedL1)]
    extendedL2 = [sorted(list(i)) for i in itertools.product(*extendedL2)]
        
    
    #Now compare the two lists and see if there is a match:
    for plist in extendedL1:
        if plist in extendedL2: return True
        
    return False
    
    
    

