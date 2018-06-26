#!/usr/bin/env python

"""
.. module:: particleNames
   :synopsis: Provides functions for getting particle objects from pdg ids or name, and
              other helpers.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.particleDefinitions import BSM
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger



def getObjectFromPdg(pdg):
    """
    Convert pdg number to particle object according to the Particle class.
    :type pdg: int
    :returns: Particles object 
    """

    for particle in BSM:
        if particle.pdg==pdg:
            return particle
    
    raise SModelSError("Particle %i not defined in particle.py" %(pdg))
    
    
    
def getObjectFromLabel(label):
    """
    Convert particle label to object according to the Particle class.
    :type name: str
    :returns: Particles object or ParticleList object
    """    

    for particle in BSM:         
        if particle.label == label:
            return particle 

    raise SModelSError("Unknown particle. Add %s to particle definitions" %label)



def getNamesList(particleList=BSM):
    """ 
	Convert list of particles to list of particle names according to the Particle class.
	If particleList is not define, use all particles defined.
	
	:type particleList: list of instances of Particles class
	:returns: list of str
	"""
    
    namesList = [particle.label for particle in particleList]
    namesList = list(set(namesList))
    
    return namesList


def getPDGList(particleList=BSM):
    """ 
    Convert list of particles to list of particle PDGs according to the Particle class.
    If particleList is not define, use all particles defined.
    
    :type particleList: list of instances of Particles class
    :returns: list of str
    """
    
    pdgsList = [particle.pdg for particle in particleList]
    pdgsList = list(set(pdgsList))
    
    return pdgsList


def elementsInStr(instring,removeQuotes=True):
    """
    Parse instring and return a list of elements appearing in instring.
    instring can also be a list of strings.
    
    :param instring: string containing elements (e.g. "[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]")
    :param removeQuotes: If True, it will remove the quotes from the particle labels.
                         Set to False, if one wants to run eval on the output.
    
    :returns: list of elements appearing in instring in string format
    
    """
    
    particleNames = getNamesList()
    
    outstr = ""
    if isinstance(instring,str):
        outstr = instring
    elif isinstance(instring,list):
        for st in instring:
            if not isinstance(st,str):
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
                if not ptc in particleNames:
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
    
    particleNames = getNamesList()
    
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
                if not ptc in particleNames:
                    logger.error("Unknown particle. Add " + ptc + " to smodels/particleDefinitions.py")
                    raise SModelSError()
            vertStr = ""

    # Check if there are not unmatched ['s and/or ]'s in the string
    if nc != 0:
        logger.error("Wrong input (incomplete elements?) " + instring)
        raise SModelSError()

    return vertices

