"""
.. module:: particleComparison
   :synopsis: Provides functions to compare particles.

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

from smodels.theory.particleNames import getNamesList
from smodels.tools.smodelsLogging import logger
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from particle import Particle,ParticleList
import itertools    
    


def simParticles(plist1, plist2):
    """
    Compares two lists of particles. Allows for particleList objects
    (Ex: L = l, l+ = l, l = l-,...). Ignores particle ordering inside
    the list.

    :param plist1: first list of Particle or ParticleList objects
    :param plist2: first list of Particle or ParticleList objects
    :param useLists: use the particle lists, i.e. allow e to stand for
                    e+ or e-, l+ to stand for e+ or mu+, etc 
    :returns: True/False if the particles lists match (ignoring order)    
    """    
    
    particleNames = getNamesList()
    
    if not isinstance(plist1,list) or type(plist1) != type(plist2):
        logger.error("Input must be a list")
        raise SModelSError()
    if len(plist1) != len(plist2):
        return False
    
    for p in plist1+plist2:              
        if not p.label in particleNames:
            raise SModelSError("Unknown particle: %s" %p)
  
    #Expand inclusive particles (ParticleList objects) 
    extendedL1 = []    
    for p in plist1:
        if isinstance(p,Particle):
            extendedL1.append([p])
        elif isinstance(p,ParticleList):
            extendedL1.append(p.particles)

     
    extendedL2 = []    
    for p in plist2:
        if isinstance(p,Particle):
            extendedL2.append([p])
        elif isinstance(p,ParticleList):
            extendedL2.append(p.particles)

    extendedL1 = [sortParticleList(list(i)) for i in itertools.product(*extendedL1)]
    extendedL2 = [sortParticleList(list(i)) for i in itertools.product(*extendedL2)]

    #Now compare the two lists and see if there is a match:
    for plist in extendedL1:
        if plist in extendedL2:
            return True  
    return False
     
     
     
def sortParticleList(ptcList):
    """
    sorts a list of particle or particle list objects by their label
    :param ptcList: list to be sorted containing particle or particle list objects
    :return: sorted list of particles
    """
    
    newPtcList = sorted(ptcList, key=lambda x: x.label) 
    return newPtcList    
