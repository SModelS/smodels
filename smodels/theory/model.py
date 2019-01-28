"""
.. module:: model
   :synopsis: Create model from file 

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
"""

import pyslha,copy
from smodels.tools.smodelsLogging import logger
from smodels.tools.physicsUnits import GeV
from smodels.theory import lheReader, crossSection
from smodels.theory.particle import MultiParticle
from smodels.theory.exceptions import SModelSTheoryError as SModelSError

class Model(object):
    """
    An instance of this class holds all the relevant information from the input model.
    This class contains the input file, the particles of the model and their production
    cross-sections.
    """
    
    def __init__(self, BSMparticles, SMparticles):
        """
        Initializes the model
        :parameter BSMparticles: list with BSM particle objects
        :parameter SMparticles: list with SM particle objects
        """
        
        
        self.inputFile = None
        self.BSMparticles = [copy.deepcopy(particle) for particle in BSMparticles]
        self.SMparticles = SMparticles[:]

        #Check if for each PDG there is a unique particle object defined  
        allPDGs = self.getValuesFor('pdg')
        for pdg in allPDGs:
            p = self.getParticlesWith(pdg=pdg)
            if len(p) > 1:
                raise SModelSError("PDG = %s has been defined for multiple particles (%s). Check your model definitions." %(pdg,p))


    def __str__(self):
        
        return self.inputFile

    def getParticlesWith(self,**kwargs):
        """
        Return the particle objects with the listed attributes.
        MultiParticle objects are added if any of its particles matches the listed attributes.
        In order to avoid double counting, MultiParticle objects are only included 
        if they do not contain any of the particles already in the list.
        For instance, if MP.particles = [e-,e+] and [e-] already appears in the list,
        MP will not be added.
        
        :returns: List of particle objects 
        """
    
        particleList = []
        for particle in self.BSMparticles + self.SMparticles:
            pList = [particle]
            if isinstance(particle,MultiParticle):
                pList += particle.particles
            for p in pList:
                if any(not hasattr(p,attr) for attr in kwargs.keys()):
                    continue
                if any(getattr(p,attr) != value for attr,value in kwargs.items()):
                    continue

                particleList.append(particle)

        #Remove repeated entries. If the list contains a Particle and a MultiParticle
        #which contains the Particle, remove the MultiParticle
        cleanList = []
        #First include all Particle objects:
        for particle in particleList:
            if isinstance(particle,MultiParticle):
                continue
            if any(particle is ptc for ptc in cleanList):
                continue
            cleanList.append(particle)
        #Now only include MultiParticle objects, if they do not contain
        #any of the particles already in the list:
        for particle in particleList:
            if not isinstance(particle,MultiParticle):
                continue
            if any(particle is ptc for ptc in cleanList):
                continue
            if any((particle.contains(p) and not particle is p) for p in particleList):
                continue
            cleanList.append(particle)


        if not cleanList:
            logger.warning("Particle with attributes %s not found in models" %str(kwargs))
        
        return cleanList
    
    def getValuesFor(self,attributeStr):
        """ 
        Returns a list with all the values for attribute appearing in the model.
         
        :param attributeStr: String for the desired attribute
         
        :returns: list of values
        """
        
        valueList = []
        for particle in self.BSMparticles+self.SMparticles:
            if not hasattr(particle,attributeStr):
                continue
            value = getattr(particle,attributeStr)
            if isinstance(particle,MultiParticle) and isinstance(value,list):
                for v in value:
                    if not v in valueList:
                        valueList.append(v)
            else:
                valueList.append(value)

        return valueList
    
    def updateParticles(self, inputFile, promptWidth = None, stableWidth = None,
                        roundMasses = 1, erasePrompt=['spin','eCharge','colordim']):
        """
        Update mass, total width and branches of allParticles particles from input SLHA or LHE file. 

        :param inputFile: input file (SLHA or LHE)
            
        :param promptWidth: Maximum width for considering particles as decaying prompt. If None, it
                            will be set 1e-8 GeV.
        :param stableWidth: Minimum width for considering particles as stable. If None, it
                            will be set 1e-25 GeV.
        :param roundMasses: If set, it will round the masses to this number of digits (int)

        :param erasePromptQNs: If set, all particles with prompt decays (totalwidth > promptWidth)
                               will have the corresponding properties (quantum numbers). So all particles with the same
                               mass and Z2parity will be considered as equal when combining elements.
            
        """
        
        self.inputFile = inputFile
        if promptWidth is None:
            promptWidth = 1e-8*GeV
        if stableWidth is None:
            stableWidth = 1e-25*GeV

        #Download input file, if requested
        if self.inputFile.startswith("http") or self.inputFile.startswith("ftp"):
            logger.info("Asked for remote slhafile %s. Fetching it." % self.inputFile)
            import requests
            import os.path
            r = requests.get(self.inputFile)
            if r.status_code != 200:
                logger.error("Could not retrieve remote file %d: %s" %(r.status_code, r.reason))
                raise SModelSError()
            basename = os.path.basename(self.inputFile)
            f = open(basename, "w")
            f.write(r.text)
            f.close()
            self.inputFile = basename

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
            for particle in self.BSMparticles:
                if not particle.pdg in massDict and not -particle.pdg in massDict:
                    continue
                BSMparticlesInEvents.append(particle)
            self.BSMparticles = BSMparticlesInEvents
        finally:
            sys.stderr = storeErr
        
        logger.debug("Loaded %i BSM particles" %len(self.BSMparticles))
        
        allPDGs = list(set(self.getValuesFor('pdg')))
        evenPDGs = []
        oddPDGs = []
        for pdg in allPDGs:
            if all(p.Z2parity == 'even' for p in self.getParticlesWith(pdg=pdg)):
                evenPDGs.append(pdg)
            elif all(p.Z2parity == 'odd' for p in self.getParticlesWith(pdg=pdg)):
                oddPDGs.append(pdg)
        
        #Remove cross-sections for even particles or particles which do not belong to the model:
        modelPDGs = [particle.pdg for particle in self.BSMparticles if isinstance(particle.pdg,int)]
        for xsec in self.xsections.xSections[:]:
            for pid in xsec.pid:
                if not pid in modelPDGs:
                    logger.debug("Cross-section for %s includes particles not belonging to model and will be ignored" %str(xsec.pid))
                    self.xsections.delete(xsec)                    
                if not pid in oddPDGs:
                    logger.debug("Cross-section for %s includes even particles and will be ignored" %str(xsec.pid))
                    self.xsections.delete(xsec)
        


        for particle in self.BSMparticles:
            if isinstance(particle,MultiParticle):
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
                particle.totalwidth = 0.*GeV
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
                logger.debug("Particle %s has width above the threshold and will be assumed as prompt." %particle.pdg)
                if erasePrompt and particle.Z2parity == 'odd':
                    logger.debug("Erasing quantum numbers of (prompt) particle %s." %particle.pdg)
                    for attr in erasePrompt:
                        delattr(particle,attr)
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
                    logger.debug("Decay %i -> %s is not of the form Z2-odd -> Z2-odd + [Z2-even particles] and will be ignored" %(pdg,pids))
                    continue

                newDecay = pyslha.Decay(br=decay.br,nda=decay.nda,parentid=decay.parentid,ids=decay.ids[:])
                #Conjugated decays if needed:
                if chargeConj == -1:
                    newDecay.ids = [pid*chargeConj if pid in allPDGs else pid for pid in decay.ids] 
                    
                    
                #Convert PDGs to particle objects:
                newDecay.daughters = []
                for pdg in newDecay.ids:
                    daughter = self.getParticlesWith(pdg=pdg)
                    if not daughter:
                        raise SModelSError("Particle with PDG = %i was not found in model. Check the model definitions." %pdg)
                    elif len(daughter) > 1:
                        raise SModelSError("Multiple particles defined with PDG = %i. PDG ids must be unique." %pdg)
                    else:
                        daughter = daughter[0]
                    newDecay.daughters.append(daughter)
                particle.decays.append(newDecay)
                
        #Reset particle equality tracking:
        for p in self.SMparticles+self.BSMparticles:
            equals = [id(p)]
            if isinstance(p,MultiParticle):
                equals += [id(ptc) for ptc in p.particles]            
            p._equals = set(equals)
            p._differs = set([])
