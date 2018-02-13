import pyslha
import copy
from smodels.particleClass_nomass import SMList, SMpdgs, SMparticles, BSMList, BSMpdgs, BSMparticles, ptcDict
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import GeV


def UpdateParticles(slhafile, BSMList):
    """
	updates mass, total width and branches of BSM particles from input slha file
	(replaces the former _getDictionariesFromSLHA)
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
					continue
				brs = []
				for decay in res.decays[pid].decays:
					nEven = nOdd = 0.
					for pidd in decay.ids:
						if pidd in rOdd: nOdd += 1
						elif pidd in rEven: nEven += 1
						else:
						    logger.warning("Particle %i not defined in particleClass.py,decay %i -> [%s] will be ignored" %(pidd,pid,decay.ids))
						    break
					if nOdd + nEven == len(decay.ids) and nOdd == 1:
						brs.append(decay)
					else:
						logger.info("Ignoring decay: %i -> [%s]",pid,decay.ids)
				particle.branches = brs

			else: print("Particle %s is not in BSMList, its mass, width and branches cannot be included") %pid



def writeIgnoreMessage ( keys, rEven, rOdd ):
    msg = ""
    for pid in keys:
        if not pid in list(rEven) + list(rOdd):
            logger.warning("Particle %i not defined in particleClass.py, its decays will be ignored" %(pid))
            continue
        if pid in rEven:
            msg += "%s, " % getObjectFromPdg(pid).label
            continue         
    if len(msg)>0:
            logger.info ( "Ignoring %s decays" % msg[:-2] )


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
    Convert particle label to particle object according to the Particle class.

    :type name: str
    :returns: Particles object 
    """    
    found = False
    for particle in SMList+BSMList:
		if particle.label == name: 
			found = True
			p = particle
    for ptcs in ptcDict:
	    if ptcs == name: 
	        found = True
            p = ptcDict[ptcs]
    if found: return p
    else: logger.warning("Particle %s not defined in particleClass.py" %(name))

