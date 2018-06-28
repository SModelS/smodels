"""
.. module:: model
   :synopsis: Create model from file 

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
"""

import pyslha,copy
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import GeV
from smodels.theory import lheReader, crossSection
from smodels.theory.particle import ParticleList,ParticleWildcard
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.particleDefinitions import SM

class Model(object):
    """
    An instance of this class holds all the relevant information from the input model.
    This class contains the input file, the particles of the model and their production
    cross-sections.
    """
    
    def __init__(self, BSMparticles, inputFile=None):
        """
        Initializes the model
        :parameter inputFile: input file
        :parameter particles: ParticleList object containing the particles of the model
        """
        self.inputFile = inputFile
        self.particles = [copy.deepcopy(particle) for particle in BSMparticles]
        
    def __str__(self):
        return self.inputFile

    def updateParticles(self, promptWidth = 1e-8*GeV, stableWidth = 1e-25*GeV):        
        self.getParticleData(promptWidth, stableWidth)
    
    def getParticleData(self, promptWidth = 1e-8*GeV, stableWidth = 1e-25*GeV):
        """
        Update mass, total width and branches of allParticles particles from input SLHA or LHE file. 
            
        :param promptWidth: Maximum width for considering particles as decaying prompt
        :param stableWidth: Minimum width for considering particles as stable
            
        """
        
        #Trick to suppress pyslha error messages:
        import sys
        storeErr = sys.stderr
        try:
            sys.stderr = None
            res = pyslha.readSLHAFile(self.inputFile)
            massDict = res.blocks['MASS'] 
            decaysDict = res.decays
            self.xsections = crossSection.getXsecFromSLHAFile(self.inputFile)                        
        except:
            massDict,decaysDict = lheReader.getDictionariesFrom(self.inputFile)
            self.xsections = crossSection.getXsecFromLHEFile(self.inputFile)
            logger.info("Using LHE input. All unstable particles will be assumed to have prompt decays.")
        sys.stderr = storeErr
        
        evenPDGs = []
        oddPDGs = []
        for particle in (self.particles+SM):
            if not hasattr(particle,'pdg') or not hasattr(particle,'Z2parity'):
                continue
            if isinstance(particle.pdg,list):
                continue
            if particle.Z2parity == 'even' and not particle.pdg in evenPDGs:
                evenPDGs.append(particle.pdg)
            elif particle.Z2parity == 'odd' and not particle.pdg in oddPDGs:
                oddPDGs.append(particle.pdg)
        allPDGs = list(set(evenPDGs + oddPDGs))

        
        #Remove cross-sections for even particles or particles which do not belong to the model:
        modelPDGs = [particle.pdg for particle in self.particles if not isinstance(particle.pdg,list)]
        for xsec in self.xsections.xSections[:]:
            
            Move pdg->ParticleOBj
            
            for pid in xsec.pid:
                if not pid in modelPDGs:
                    logger.debug("Cross-section for %s includes particles not belonging to model and will be ignored" %str(xsec.pid))
                    self.xsections.delete(xsec)                    
                if not pid in oddPDGs:
                    logger.debug("Cross-section for %s includes even particles and will be ignored" %str(xsec.pid))
                    self.xsections.delete(xsec)
        


        for particle in self.particles:
            if isinstance(particle,(ParticleList,ParticleWildcard)):
                continue

            if not hasattr(particle,'pdg') or not hasattr(particle,'Z2parity'):
                raise SModelSError("PDG and/or Z2-parity for particle %s has not been defined" %particle.label)

            pdg = particle.pdg            
            if abs(pdg) in massDict:
                particle.mass = abs(massDict[abs(pdg)])*GeV
            else:
                logger.debug("No mass found for %i. Its mass will be set to None." %particle.pdg)
                particle.mass = None

            particle.decays = []
            if not pdg in decaysDict and not -pdg in decaysDict:
                logger.debug("No decay found for %i. It will be considered stable" %particle.pdg)                
                particle.width = 0.*GeV
                continue
            
            if pdg in decaysDict:
                particleData = decaysDict[pdg]
                chargeConj = 1
            else:
                particleData = decaysDict[-pdg]
                chargeConj = -1
                
            particle.totalwidth = abs(particleData.totalwidth)*GeV
            if particle.totalwidth < stableWidth:
                particle.totalwidth = 0.*GeV  #Treat particle as stable
                logger.debug("Particle %s has width below the threshold and will be assumed as stable" %particle.pdg)
                continue

            if particle.totalwidth > promptWidth:
                particle.totalwidth = float('inf')  #Treat particle as prompt
                logger.debug("Particle %s has width above the threshold and will be assumed as prompt" %particle.pdg)
            else:
                particle.decays.append(None) #Include possibility for particle being long-lived (non-prompt)
            
            for decay in particleData.decays:
                pids = decay.ids
                missingIDs = set(pids).difference(set(allPDGs))
                if missingIDs:
                    logger.info("Particle(s) %s is not defined within model. Decay %s will be ignored" %(missingIDs,decay))
                    continue
                oddPids = [pid for pid in decay.ids if abs(pid) in oddPDGs]
                evenPids = [pid for pid in decay.ids if abs(pid) in evenPDGs]
                if len(oddPids) != 1 or len(evenPids+oddPids) != len(decay.ids):
                    logger.info("Decay %s is not of the form Z2-odd -> Z2-odd + [Z2-even particles] and will be ignored" %(decay))
                    continue
                
                Move pdg->ParticleOBj
                
                newDecay = pyslha.Decay(br=decay.br,nda=decay.nda,parentid=decay.parentid,ids=decay.ids[:])
                #Conjugated decays if needed:
                if chargeConj == -1:
                    newDecay.ids = [pid*chargeConj if pid in allPDGs else pid for pid in decay.ids] 
                particle.decays.append(newDecay)
                 
        
