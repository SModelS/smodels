"""
.. module:: particles
   :synopsis: Updates masses, width and branches from SLHA file

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
"""

import pyslha
import copy
from math import exp
from smodels.particleDefinitions import SMpdgs, BSMList, BSMpdgs
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import GeV, MeV, m, mm, fm


def updateParticles(slhafile, BSMList):
    getMassWidthBranches(slhafile, BSMList)
    getPromptDecays(BSMList)
    

def getMassWidthBranches(slhafile, BSMList):
    """
    Update mass, total width and branches of BSM particles from input slha file. 
        
    :param BSMList: list of particle objects containing all BSM particles defined in particleDefinitions
        
    """


    res = pyslha.readSLHAFile(slhafile)

    rOdd = BSMpdgs
    rEven = SMpdgs

    writeIgnoreMessage ( res.decays.keys(), rEven, rOdd )

	    
    for pid in res.decays.keys():
        for particle in BSMList:
            if abs(particle.pdg) == pid:
                particle.mass = abs(res.blocks['MASS'][pid])*GeV
                particle.width = abs(res.decays[pid].totalwidth)*GeV

                if not pid in rOdd:
                    logger.warning("Particle %s is not in BSMList, its mass, width and branches cannot be included") %pid
                    continue
                brs = []
                for decay in res.decays[pid].decays:
                    nEven = nOdd = 0.
                    for pidd in decay.ids:
                        if pidd in rOdd: nOdd += 1
                        elif pidd in rEven: nEven += 1
                        else:
                            logger.warning("Particle %i not defined in particleDefinitions.py,decay %i -> [%s] will be ignored" %(pidd,pid,decay.ids))
                            break
                    if nOdd + nEven == len(decay.ids) and nOdd == 1:
                        brs.append(decay)
                    else:
                        logger.info("Ignoring decay: %i -> [%s]",pid,decay.ids)
                
                if particle.pdg == pid: particle.branches = brs
                elif (-1)*particle.pdg == pid:
                    brsConj = copy.deepcopy(brs)
                    for br in brsConj:
                        br.ids = [-x for x in br.ids]                              
                    particle.branches = brsConj

    #getPromptDecays(BSMList)


def writeIgnoreMessage ( keys, rEven, rOdd ):
    msg = ""
    for pid in keys:
        if not pid in list(rEven) + list(rOdd):
            logger.warning("Particle %i not defined in particleDefinitions.py, its decays will be ignored" %(pid))
            continue
        if pid in rEven:
            msg += "%i, " %pid # used to be string of particle name (getObjectFromPdg(pid).label causes circular dependence) -- pid not so pretty 
            continue         
    if len(msg)>0:
            logger.info ( "Ignoring %s decays" % msg[:-2] )
            
            
            
def getPromptDecays(BSMList):
    """
    Update br including whether a particle decays promptly or is long-lived.
    """
    l_inner = 10.*mm
    l_outer = 10.*m
    gb_inner = 10.
    gb_outer = 0.6    
    
    hc = 197.327*MeV*fm  #hbar * c

    for particle in BSMList:    
        F_prompt = 1. - exp( -1*particle.width * l_inner /(gb_inner*hc) )       
        F_long = exp( -1*particle.width * l_outer /(gb_outer*hc) )

        for branch in particle.branches:
            branch.br *= F_prompt
        
        if F_long:
            stable = pyslha.Decay(br=F_long, nda=0 ,ids=[])
            particle.branches.append(stable)
            
        if (F_long+F_prompt) > 1.:
            logger.error("Sum of decay fractions > 1 for "+str(particle.pdg))
            return False    
            

