"""
.. module:: model
   :synopsis: Create model from file 

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
"""

import pyslha
from smodels.particleDefinitions import allParticles, SM
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import GeV
from smodels.theory.particleNames import getObjectFromPdg, getPDGList
from smodels.theory import lheReader, crossSection


class Model(object):
    """
    An instance of this class holds all the relevant information from the input model.
    This class contains the input file, the particles of the model and their production
    cross-sections.
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
        self.getParticleData()
    
    def getParticleData(self, promptWidth = 1e-8*GeV, stableWidth = 1e-25*GeV):
        """
        Update mass, total width and branches of allParticles particles from input SLHA or LHE file. 
            
        :param promptWidth: Maximum width for considering particles as decaying prompt
        :param stableWidth: Minimum width for considering particles as stable
            
        """
    
        try:
            res = pyslha.readSLHAFile(self.inputFile)
            massDict = res.blocks['MASS'] 
            decaysDict = res.decays
            self.xsections = crossSection.getXsecFromSLHAFile(self.inputFile)
        except:
            massDict,decaysDict = lheReader.getDictionariesFrom(self.inputFile)
            self.xsections = crossSection.getXsecFromLHEFile(self.inputFile)
            
        SMpdgs = getPDGList(SM)
        allpdgs = getPDGList()
        BSMpdgs = [pdg for pdg in allpdgs if not abs(pdg) in SMpdgs]
        evenPDGs = [particle.pdg for particle in allParticles if particle.Z2parity == 'even']
        oddPDGs = [particle.pdg for particle in allParticles if particle.Z2parity == 'odd']

        for pdg in massDict:
            if not pdg in allpdgs:
                logger.info("Particle %s is not in defined and will be ignored") %pdg
                continue
            
            if not pdg in BSMpdgs:
                logger.debug("PDG %i belongs to the SM and will be ignored.")
                continue
            
            if not pdg in oddPDGs:
                logger.debug("PDG %i is Z2-even and will be ignored.")
                continue
                
            particles = [getObjectFromPdg(pdg),getObjectFromPdg(-pdg)]
            for particle in particles:
                if not particle:
                    continue
                particle.mass = abs(massDict[pdg])*GeV

            particle.decays = []
            if not pdg in decaysDict:
                logger.debug("No decay found for %i. It will be considered stable" %particle.pdg)                
                particle.width = 0.*GeV
                continue
            
            decays = decaysDict[pdg]
            particle.totalwidth = abs(decays.totalwidth)*GeV
            if particle.totalwidth < stableWidth:
                particle.totalwidth = 0.*GeV  #Treat particle as stable
                logger.debug("Particle %s has width below the threshold and will be assumed as stable" %particle.pdg)
                continue

            if particle.totalwidth > promptWidth:
                particle.totalwidth = float('inf')  #Treat particle as prompt
                logger.debug("Particle %s has width above the threshold and will be assumed as prompt" %particle.pdg)
            
            chargeConj = pdg/particle.pdg # = +1 for particle and -1 for anti-particle            
            for decay in decays:
                pids = decay.ids
                missingIDs = set(pids).difference(set(allpdgs))
                if missingIDs:
                    logger.info("Particle(s) %s is not define. Decay %s will be ignored" %(missingIDs,decay))
                    continue
                oddPids = [pid for pid in decay.ids if abs(pid) in oddPDGs]
                evenPids = [pid for pid in decay.ids if abs(pid) in evenPDGs]
                if len(oddPids) != 1 or len(evenPids+oddPids) != len(decays.ids):
                    logger.info("Decay %s is not of the form Z2-odd -> Z2-odd + [Z2-even particles] and will be ignored" %(decay))
                    continue
                particle.decays.append(decay)
                particle.decays[-1].ids = [pid*chargeConj for pid in decay.ids] #Conjugated decays
                 
        
