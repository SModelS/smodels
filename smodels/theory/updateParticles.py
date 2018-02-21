"""
.. module:: particles
   :synopsis: Updates masses, width and branches from SLHA file

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
"""

import pyslha
from smodels.particleDefinitions import SMpdgs, BSMList, BSMpdgs
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import GeV
#from smodels.theory.particleNames import getObjectFromPdg



def updateParticles(slhafile, BSMList):
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
                particle.mass = round(abs(res.blocks['MASS'][pid]),1)*GeV
                particle.width = round(abs(res.decays[pid].totalwidth),1)*GeV

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
                particle.branches = brs



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
            

