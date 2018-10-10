"""
.. module:: model
   :synopsis: Create model from file 

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
"""

import pyslha,copy
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import GeV
from smodels.theory import lheReader, crossSection
from smodels.theory.particle import ParticleList,InclusiveParticle
from smodels.theory.exceptions import SModelSTheoryError as SModelSError

class Model(object):
    """
    An instance of this class holds all the relevant information from the input model.
    This class contains the input file, the particles of the model and their production
    cross-sections.
    """
    
    def __init__(self, BSMparticles, SMparticles, inputFile=None):
        """
        Initializes the model
        :parameter inputFile: input file (SLHA or LHE)
        :parameter BSMparticles: list with BSM particle objects
        :parameter SMparticles: list with SM particle objects
        """
        
        
        self.inputFile = inputFile
        self.oddParticles = [copy.deepcopy(particle) for particle in BSMparticles]
        self.SMparticles = SMparticles[:]
        
    def __str__(self):
        
        return self.inputFile

    def getParticlesWith(self,**kwargs):
        """
        Return the particle objects with the listed attributes.
        
        :returns: List of particle objects 
        """
    
        particleList = []
        for particle in self.oddParticles + self.SMparticles:            
            if any(not hasattr(particle,attr) for attr in kwargs.keys()):
                continue
            if any(getattr(particle,attr) != value for attr,value in kwargs.items()):
                continue

            #Avoid double counting:
            if not any(particle is ptc for ptc in particleList):
                particleList.append(particle)
        
        if not particleList:
            logger.warning("Particle with attributes %s not found in models" %str(kwargs))
        
        return particleList
    
    def getValuesFor(self,attributeStr):
        """ 
        Returns a list with all the values for attribute appearing in the model.
         
        :param attributeStr: String for the desired attribute
         
        :returns: list of values
        """
        
        valueList = []
        for particle in self.oddParticles+self.SMparticles:
            if not hasattr(particle,attributeStr):
                continue
            value = getattr(particle,attributeStr)
            if not value in valueList:
                valueList.append(value)
        
        return valueList
    

    def updateParticles(self, promptWidth = 1e-8*GeV, stableWidth = 1e-25*GeV, roundMasses = 1):        
        """
        Update mass, total width and branches of allParticles particles from input SLHA or LHE file. 
            
        :param promptWidth: Maximum width for considering particles as decaying prompt
        :param stableWidth: Minimum width for considering particles as stable
        :param roundMasses: If set, it will round the masses to this number of digits (int)  
            
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
            logger.info("Using LHE input. All particles not appearing in the events will be removed from the model.")
            BSMparticlesInEvents = []
            for particle in self.oddParticles:
                if not particle.pdg in massDict and not -particle.pdg in massDict:
                    continue
                BSMparticlesInEvents.append(particle)
            self.oddParticles = BSMparticlesInEvents
                
        sys.stderr = storeErr
        
        evenPDGs = list(set([ptc.pdg for ptc in self.getParticlesWith(Z2parity='even') if isinstance(ptc.pdg,int)]))
        oddPDGs = list(set([ptc.pdg for ptc in self.getParticlesWith(Z2parity='odd') if isinstance(ptc.pdg,int)]))
        allPDGs = list(set(evenPDGs + oddPDGs))
        
        #Remove cross-sections for even particles or particles which do not belong to the model:
        modelPDGs = [particle.pdg for particle in self.oddParticles if isinstance(particle.pdg,int)]
        for xsec in self.xsections.xSections[:]:
            for pid in xsec.pid:
                if not pid in modelPDGs:
                    logger.debug("Cross-section for %s includes particles not belonging to model and will be ignored" %str(xsec.pid))
                    self.xsections.delete(xsec)                    
                if not pid in oddPDGs:
                    logger.debug("Cross-section for %s includes even particles and will be ignored" %str(xsec.pid))
                    self.xsections.delete(xsec)
        


        for particle in self.oddParticles:
            if isinstance(particle,(ParticleList,InclusiveParticle)):
                continue

            if not hasattr(particle,'pdg') or not hasattr(particle,'Z2parity'):
                raise SModelSError("PDG and/or Z2-parity for particle %s has not been defined" %particle.label)

            pdg = particle.pdg  
            
            if abs(pdg) in massDict.keys():
                if roundMasses and int(roundMasses) > 0:
                    particle.mass = round(abs(massDict[abs(pdg)]),int(roundMasses))*GeV
                else:
                    particle.mass = abs(massDict[abs(pdg)])*GeV
            else:
                raise SModelSError("No mass found for %i in input file %s." %(particle.pdg,self.inputFile))

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
                particle.totalwidth = float('inf')*GeV  #Treat particle as prompt
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

                newDecay = pyslha.Decay(br=decay.br,nda=decay.nda,parentid=decay.parentid,ids=decay.ids[:])
                #Conjugated decays if needed:
                if chargeConj == -1:
                    newDecay.ids = [pid*chargeConj if pid in allPDGs else pid for pid in decay.ids] 
                    
                    
                #Convert PDGs to particle objects:
                newDecay.daughters = [self.getParticlesWith(pdg=pdg)[0] for pdg in newDecay.ids]
                particle.decays.append(newDecay)
                 
        
