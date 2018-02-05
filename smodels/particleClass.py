import pyslha
import copy
from smodels.particleClass_nomass import SMList, SMpdgs, SMparticles, BSMList, BSMpdgs, BSMparticles, ptcDic
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import GeV


def UpdateParticles(inputfile, BSMList):
	brDic, massDic, widthDic = _getDictionariesFromSLHA(inputfile)	
		
	for pdg in massDic.keys():
		for particle in BSMList:
			if particle.pdg == pdg:
				particle.mass = massDic[pdg]
				particle.width = widthDic[pdg]
				particle.branches = brDic[pdg]
			else: print("Particle %s is not in BSMList, its mass cannot be included") %pdg
	



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

def _getDictionariesFromSLHA(slhafile):
    """
    Create mass and BR dictionaries from an SLHA file.
    Ignore decay blocks with R-parity violating or unknown decays
    """

    res = pyslha.readSLHAFile(slhafile)

    rOdd = BSMpdgs
    rEven = SMpdgs
    
    # Get mass and branching ratios for all particles
    brDic = {}
    writeIgnoreMessage ( res.decays.keys(), rEven, rOdd )

    for pid in res.decays.keys():
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

        brsConj = copy.deepcopy(brs)
        for br in brsConj:
            br.ids = [-x for x in br.ids]
        brDic[pid] = brs
        brDic[-pid] = brsConj
    # Get mass list for all particles
    massDic = dict(res.blocks['MASS'].items())
    for pid in list ( massDic.keys() )[:]:
        massDic[pid] = round(abs(massDic[pid]),1)*GeV
        if not -pid in massDic: massDic[-pid] = massDic[pid]    

    # Get total width for all particles
    widthlist = [res.decays[i].totalwidth for i in res.decays.keys()]
    widthDic = dict(zip(res.decays.keys(), widthlist))

    for pid in list ( widthDic.keys() )[:]:
        widthDic[pid] = round(abs(widthDic[pid]),1)*GeV
        if not -pid in widthDic: widthDic[-pid] = widthDic[pid]  

    return brDic, massDic, widthDic


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
    else: return False

