"""
.. module:: particleComparison
   :synopsis: Provides functions to compare particles.

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
"""

from smodels.particleDefinitions import particleLists
from smodels.theory.particleClass import particleInList, sortParticleList
from smodels.tools.smodelsLogging import logger
import itertools    
    
def compareBSMparticles(BSMptcs1, BSMptcs2):
    """
    Compares intermediate particles of two branches.
    Comparison is done based on masses.
    :return: True/False, intermediate particles are/ are not equal
             comp representing whether the first masses are larger than the second
    """
    
    plist1 = getBSMparticles(BSMptcs1)    
    plist2 = getBSMparticles(BSMptcs2)    
    
    masses1 = [ particle.mass for particle in plist1 ]      
    masses2 = [ particle.mass for particle in plist2 ] 
    
    m1m2eq = ( masses1 == masses2 )
    comp = ( masses1 > masses2 )
        
    return m1m2eq, comp  


def getBSMparticles(BSMptcs):  
    """
    Reduces the nested list structure of the BSM particles in branch to a list of particles
    :return: list of particle objects
    """
      
    if len(BSMptcs) == 1:
        particles = [particle for particleList in BSMptcs for particle in particleList ]
    elif len(BSMptcs) > 1:
        particles = [particle for particle in BSMptcs[0] ]
    else: particles = []
    return particles  
    



def simParticles(plist1, plist2, useLists=True):
    """
    Compares two lists of particles. Allows for particleList objects
    (Ex: L = l, l+ = l, l = l-,...). Ignores particle ordering inside
    the list.

    :param plist1: first list of particles/ particle lists, e.g. [l,jet]
    :param plist2: second list of particles/ particle lists
    :param useLists: use the particle lists, i.e. allow e to stand for
                    e+ or e-, l+ to stand for e+ or mu+, etc 
    :returns: True/False if the particles lists match (ignoring order)    
    """    
    
    if not isinstance(plist1,list) or type(plist1) != type(plist2):
        logger.error("Input must be a list")
        raise SModelSError()
    if len(plist1) != len(plist2):
        return False
    
    
                        
    for i,p in enumerate(plist1):
        if not isinstance(p.label,str) or not isinstance(plist2[i].label,str) :
            logger.error("Input must be two lists of particle or particle list objects which have a label")
            raise SModelSError()
  
    l1 = sortParticleList(plist1) 
    l2 = sortParticleList(plist2)
    if not useLists:
        return l1 == l2
        
    #If lists are to be used, replace particles by their list entries (list.particles)
    #e.g. [jet,mu+] -> [[q,g,c],[mu+]], [jet,mu] -> [[q,g,c],[mu+,mu-]] 

    extendedL1 = []    
    for i,p in enumerate(plist1):
        if not particleInList(p,particleLists, inlist=False): 
            extendedL1.append([p])
        else:
            extendedL1.append( p.particles )

     
    extendedL2 = []    
    for i,p in enumerate(plist2):
        if not particleInList(p,particleLists, inlist=False):
            extendedL2.append([p])
        else:    
            extendedL2.append( p.particles )  

    extendedL1 = [sortParticleList(list(i)) for i in itertools.product(*extendedL1)]
    extendedL2 = [sortParticleList(list(i)) for i in itertools.product(*extendedL2)]

    #Now compare the two lists and see if there is a match:
    for plist in extendedL1:
        if plist in extendedL2: return True  
    return False
    

    
