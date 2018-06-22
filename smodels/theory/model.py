"""
.. module:: model
   :synopsis: Create model from file 

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
"""

import pyslha
import copy
from math import exp
from smodels.particleDefinitions import SMpdgs
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import GeV, MeV, m, fm


class Model(object):
    """
    An instance of this class represents a BSM model
    This class contains the input file and the particles of the model
    """
    def __init__(self, inputFile, BSMparticleList):
        """
        Initializes the model
        :parameter inputFile: input file
        :parameter BSMparticleList: ParticleList object containing the particles of the model
        """
        self.inputFile = inputFile
        self.BSMparticleList = BSMparticleList
        
    def __str__(self):
        return self.inputFile
    

    def updateParticles(self):
        self.getMassWidthBranches()
        self.addLongLived()
        return self
        
    
    def getMassWidthBranches(self):
        """
        Update mass, total width and branches of BSM particles from input slha file. 
            
        :param BSMparticleList: ParticleList object containing all BSM particles defined in particleDefinitions
            
        """
    
        res = pyslha.readSLHAFile(self.inputFile)
        BSMpdgs = self.BSMparticleList.getPdgs()
    
        writeIgnoreMessage ( res.decays.keys(), SMpdgs, BSMpdgs )
    
            
        for pid in res.decays.keys():
            for particle in self.BSMparticleList.particles:
                if abs(particle.pdg) == pid:
                    particle.mass = abs(res.blocks['MASS'][pid])*GeV
                    particle.width = abs(res.decays[pid].totalwidth)*GeV
    
                    if not pid in BSMpdgs:
                        logger.warning("Particle %s is not in BSMparticleList.particles, its mass, width and branches cannot be included") %pid
                        continue
                    brs = []
                    for decay in res.decays[pid].decays:
                        nEven = nOdd = 0.
                        for pidd in decay.ids:
                            if pidd in BSMpdgs: nOdd += 1
                            elif pidd in SMpdgs: nEven += 1
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
                
                 
                 
    def addLongLived(self):
        """
        Update BSM particles including whether a particle can be long-lived.
        
        Using the widths of the particles (see getMassWidthBranches), reweights the BR's with the fraction
        of decays within the detector and add the fraction of "long-lived decays".
        The fraction of decays within the detector and "long-lived decays" are defined as:
      
        F_long = exp(-width*l_outer/gb_outer)
        F_decay = 1 - F_long
        
        where l_outer is the outer radius of the detector
        and gb_outer is the estimate for the kinematical factor gamma*beta.
        We use l_outer=10.m and gb_outer= 0.6    
        """
        
        l_outer = 10.*m
        gb_outer = 0.6    
        
        hc = 197.327*MeV*fm  #hbar * c
        
        for particle in self.BSMparticleList.particles:       
            F_long = exp( -1*particle.width * l_outer /(gb_outer*hc) )
            """
            F_decay = (1. - F_long) 
            
            
            for branch in particle.branches:
                branch.br *= F_decay            
            """
            if F_long:
                stable = pyslha.Decay(br=1., nda=0 ,ids=[]) #F_long
                particle.branches.append(stable)
                
                
                
def writeIgnoreMessage ( keys, SMpdgs, BSMpdgs ):
    msg = ""
    for pid in keys:
        if not pid in list(SMpdgs) + list(BSMpdgs):
            logger.warning("Particle %i not defined in particleDefinitions.py, its decays will be ignored" %(pid))
            continue
        if pid in SMpdgs:
            msg += "%i, " %pid # used to be string of particle name (getObjectFromPdg(pid).label causes circular dependence) -- pid not so pretty 
            continue         
    if len(msg)>0:
            logger.info ( "Ignoring %s decays" % msg[:-2] )                
