"""
.. module:: model
   :synopsis: Create model from file

.. moduleauthor:: Alicia Wongel <alicia.wongel@gmail.com>
"""

import pyslha
import os
from smodels.base.smodelsLogging import logger
from smodels.base.physicsUnits import GeV
from smodels.base import lheReader, crossSection
from smodels.base.particle import Particle, MultiParticle
from smodels.base.exceptions import SModelSBaseError as SModelSError

class Model(object):
    """
    An instance of this class holds all the relevant information from the input model.
    This class contains the input file, the particles of the model and their production
    cross-sections.
    """

    def __init__(self, BSMparticles, SMparticles, label=None):
        """
        Initializes the model
        :parameter BSMparticles: list with BSM particle objects
        :parameter SMparticles: list with SM particle objects
        :parameter label: Optional string to label the model
        """

        self.inputFile = None
        # allBSMparticles stores all possible particles for a given model:
        self.allBSMparticles = BSMparticles[:]
        # BSMparticles is used to store the model particles defined by the inputFile.
        # When a model is updated, a copy of the BSM particles are stored in BSMparticles
        # and these are used to store masses, decays,...
        self.BSMparticles = BSMparticles[:]  # at initialization assume all particles will be used
        self.SMparticles = SMparticles[:]
        for p in self.SMparticles:
            p.isSM = True
        for p in self.BSMparticles:
            p.isSM = False

        self.label = label

        #  Check if for each PDG there is a unique particle object defined
        allPDGs = self.getValuesFor('pdg')
        for pdg in allPDGs:
            p = self.getParticlesWith(pdg=pdg)
            if len(p) > 1:
                raise SModelSError(f"PDG = {pdg} has been defined for multiple particles ({p}). Check your model definitions.")

    def __str__(self):
        if not self.label:
            return f'Model: {self.inputFile}'
        else:
            return f'Model: {self.label}'

    def __repr__(self):

        return str(self)

    def __eq__(self, other):
        """
        Simple model comparison.

        Check if all the particles are the same (including their ordering)
        and if the model label is the same
        """

        if self.label != other.label:
            return False
        if len(self.SMparticles) != len(other.SMparticles):
            return False
        if len(self.BSMparticles) != len(other.BSMparticles):
            return False
        if self.SMparticles != other.SMparticles:
            return False
        if self.BSMparticles != other.BSMparticles:
            return False

        return True

    def getParticle(self, **kwargs):
        """
        Return a single particle object with the listed attributes.
        If no particle is found or more than one particle is found, raise an error.

        :returns: Particle object
        """

        particleList = self.getParticlesWith(**kwargs)
        if not particleList:
            raise SModelSError("Particle with attributes %s has not been found in model %s"
                               % (kwargs, self))
        elif len(particleList) > 1:
            raise SModelSError("Multiple particles with attributes %s found in model %s"
                                   % (kwargs, self))
        else:
            return particleList[0]

    def getParticlesWith(self, **kwargs):
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
            if isinstance(particle, MultiParticle):
                pList += particle.particles
            for p in pList:
                if any(not hasattr(p, attr) for attr in kwargs.keys()):
                    continue
                if any(getattr(p, attr) != value for attr, value in kwargs.items()):
                    continue

                particleList.append(particle)

        # Remove repeated entries. If the list contains a Particle and a MultiParticle
        # which contains the Particle, remove the MultiParticle
        cleanList = []
        # First include all Particle objects:
        for particle in particleList:
            if isinstance(particle, MultiParticle):
                continue
            if any(particle is ptc for ptc in cleanList):
                continue
            cleanList.append(particle)
        # Now only include MultiParticle objects, if they do not contain
        # any of the particles already in the list:
        for particle in particleList:
            if not isinstance(particle, MultiParticle):
                continue
            if any(particle is ptc for ptc in cleanList):
                continue
            if any((particle.contains(p) and not particle is p) for p in particleList):
                continue
            cleanList.append(particle)

        if not cleanList:
            logger.warning(f"Particle with attributes {str(kwargs)} not found in {self}")

        return cleanList

    def getValuesFor(self, attributeStr):
        """
        Returns a list with all the values for attribute appearing in the model.

        :param attributeStr: String for the desired attribute

        :returns: list of values
        """

        valueList = []
        for particle in self.BSMparticles+self.SMparticles:
            if not hasattr(particle, attributeStr):
                continue
            value = getattr(particle, attributeStr)
            if isinstance(particle, MultiParticle) and isinstance(value, list):
                for v in value:
                    if not v in valueList:
                        valueList.append(v)
            else:
                valueList.append(value)

        return valueList

    def getModelDataFrom(self, inputFile):
        """
        Reads the input file (LHE or SLHA) and extract the relevant information
        (masses, widths, BRs and cross-sections). If a http address is given, it will
        attempt to download the file.

        :param inputFile: input file (SLHA or LHE), can also be a string containing the SLHA file

        :return: dictionary with masses, dictionary with decays and XSectionList object
        """

        # Download input file, if requested
        if inputFile.startswith("http") or inputFile.startswith("ftp"):
            logger.info(f"Asked for remote slhafile {inputFile}. Fetching it.")
            import requests
            r = requests.get(inputFile)
            if r.status_code != 200:
                logger.error("Could not retrieve remote file %d: %s" % (r.status_code, r.reason))
                raise SModelSError()
            baseName = os.path.basename(inputFile)
            f = open(baseName, "w")
            f.write(r.text)
            f.close()
            inputFile = baseName

        #  Trick to suppress pyslha error messages:
        import sys
        storeErr = sys.stderr
        #  Try to read file assuming it is an SLHA file:
        try:
            sys.stderr = None
            if os.path.isfile(inputFile):
                res = pyslha.readSLHAFile(inputFile)
            else:
                res = pyslha.readSLHA(inputFile)
            massDict = res.blocks['MASS']
            #  Make sure both PDG signs appear in massDict
            for pdg, mass in massDict.items():
                if not -pdg in massDict:
                    massDict[-pdg] = abs(mass)
            decaysDict = res.decays
            xsections = crossSection.getXsecFromSLHAFile(inputFile)
        #  If fails assume it is an LHE file:
        except (IOError, AttributeError, KeyError):
            massDict, decaysDict = lheReader.getDictionariesFrom(inputFile)
            xsections = crossSection.getXsecFromLHEFile(inputFile)
            logger.info("Using LHE input. All unstable particles will be assumed to have prompt decays.")
            logger.info("Using LHE input. All particles not appearing in the events will be removed from the model.")
        finally:
            sys.stderr = storeErr

        return massDict, decaysDict, xsections

    def getSMandBSMList(self):
        """
        Get the list of SM and BSM particles, according to the isSM value
        defined for each particle.

        :return: list with PDGs of even particles, list with PDGs of odd particles
        """

        allPDGs = list(set(self.getValuesFor('pdg')))
        smPDGs = []
        bsmPDGs = []
        for pdg in allPDGs:
            if all(p.isSM for p in self.getParticlesWith(pdg=pdg)):
                smPDGs.append(pdg)
            elif all(not p.isSM for p in self.getParticlesWith(pdg=pdg)):
                bsmPDGs.append(pdg)

        return smPDGs, bsmPDGs

    def filterCrossSections(self):
        """
        Remove cross-sections for even particles or particles which do not belong to the model.
        Valid cross-sections are stored in self.xsections

        :return: Number of cross-sections after filter.
        """

        smPDGs, bsmPDGs = self.getSMandBSMList()
        allPDGs = smPDGs+bsmPDGs
        for xsec in self.xsections.xSections[:]:
            if any(pid not in allPDGs for pid in xsec.pid):
                logger.debug(f"Cross-section for {str(xsec.pid)} includes particles not belonging to model and will be ignored")
                self.xsections.delete(xsec)
            if all(pid in smPDGs for pid in xsec.pid):
                logger.debug(f"Cross-section for {str(xsec.pid)} includes only SM particles and will be ignored")
                self.xsections.delete(xsec)

        return len(self.xsections.xSections)

    def setMasses(self,massDict,roundMasses,minMass):
        """
        Define particle masses using massDict.

        :param roundMasses: If set, it will round the masses to this number of digits (int)
        :param massDict: dictionary with PDGs as keys and masses as values.
        :param minMass: Minimal mass for BSM particles in the model. Any particle with mass below minMass will
                        have its mass set to minMass.
        """

        for particle in self.BSMparticles:
            if isinstance(particle, MultiParticle):
                continue

            if not hasattr(particle, 'pdg'):
                raise SModelSError(f"PDG for particle {particle.label} has not been defined")

            pdg = particle.pdg
            if abs(pdg) in massDict.keys():
                if roundMasses and int(roundMasses) > 0:
                    particle.mass = round(abs(massDict[abs(pdg)]), int(roundMasses))*GeV
                else:
                    particle.mass = abs(massDict[abs(pdg)])*GeV
                if minMass > 0*GeV:
                    particle.mass = max(particle.mass,minMass)
            else:
                raise SModelSError("No mass found for %i in input file %s." % (particle.pdg, self.inputFile))

    def setDecays(self, decaysDict, promptWidth, stableWidth, ignorePromptQNumbers):

        allPDGs = list(set(self.getValuesFor('pdg')))

        for particle in self.BSMparticles:
            if isinstance(particle, MultiParticle):
                continue

            if not hasattr(particle, 'pdg'):
                raise SModelSError(f"PDG for particle {particle.label} has not been defined")

            pdg = particle.pdg
            particle.decays = []
            if pdg in decaysDict:
                particleData = decaysDict[pdg]
                chargeConj = 1
            elif -pdg in decaysDict:
                particleData = decaysDict[-pdg]
                chargeConj = -1
            else:
                logger.error("Decay information for particle %i could not be found" % pdg)
                raise SModelSError()

            particle.totalwidth = abs(particleData.totalwidth)*GeV
            if particle.totalwidth < stableWidth:
                particle._isStable = True  # Treat particle as stable
                logger.debug(f"Particle {particle.pdg} has width below the threshold and will be assumed as stable")
                continue

            if particle.totalwidth > promptWidth:
                particle._isPrompt = True  # Treat particle as prompt
                logger.debug(f"Particle {particle.pdg} has width above the threshold and will be assumed as prompt.")
            else:
                particle.decays.append(None)  # Include possibility for particle being long-lived (non-prompt)

            for decay in particleData.decays:
                pids = decay.ids
                missingIDs = set(pids).difference(set(allPDGs))
                if missingIDs:
                    logger.info(f"Particle(s) {missingIDs} is not defined within model. Decay {decay} will be ignored")
                    continue

                #  Conjugated decays if needed
                #  (if pid*chargeConj is not in model, assume the particle is its own anti-particle)
                decayIDs = [pid*chargeConj if pid*chargeConj in allPDGs else pid for pid in decay.ids]
                newDecay = pyslha.Decay(br=decay.br, nda=decay.nda, parentid=decay.parentid, ids=decayIDs)

                #  Convert PDGs to particle objects:
                daughters = []
                for pid in newDecay.ids:
                    daughter = self.getParticlesWith(pdg=pid)
                    if not daughter:
                        raise SModelSError("Particle with PDG = %i was not found in model. Check the model definitions." % pid)
                    elif len(daughter) > 1:
                        raise SModelSError("Multiple particles defined with PDG = %i. PDG ids must be unique." % pid)
                    else:
                        daughter = daughter[0]
                    daughters.append(daughter)
                newDecay.daughters = daughters
                particle.decays.append(newDecay)

            # Check for broad width resonances
            broadWidth = 0.02
            if particle.totalwidth > broadWidth*particle.mass:
                for decay in particle.decays:
                    # Check if the particle can decay to SM only:
                    if decay is not None and all(daughter.isSM for daughter in decay.daughters):
                    # if all(daughter.isSM for daughter in decay.daughters):
                        logger.warning(f"Particle {str(particle)} has a total width/mass = {float(particle.totalwidth / particle.mass):1.2f}. Some results may not be valid for broad resonances!")
                        break


        #  Check if all unstable particles have decay channels defined:
        for p in self.BSMparticles:
            if p.totalwidth < stableWidth:
                continue
            ndecays = len([dec for dec in p.decays if dec is not None])
            if ndecays == 0:
                if not p.isSM:
                    logger.warning(f"No valid decay found for {p}. It will be considered stable.")
                p.totalwidth = 0.*GeV

        # Finally erase attributes of prompt particles
        if ignorePromptQNumbers:
            for particle in self.BSMparticles:
                 if not particle.isSM and particle.totalwidth > promptWidth:
                    logger.debug(f"Erasing quantum numbers of (prompt) particle {particle.pdg}.")
                    for attr in ignorePromptQNumbers:
                        if hasattr ( particle, attr ):
                            delattr(particle, attr)

    def updateParticles(self, inputFile, promptWidth = None, stableWidth = None,
                        roundMasses = 1, ignorePromptQNumbers=[],
                        minMass=1.0*GeV):
        """
        Update mass, total width and branches of allParticles particles from input SLHA or LHE file.

        :param inputFile: input file (SLHA or LHE), can also be the file content given as a long string

        :param promptWidth: Maximum width for considering particles as decaying prompt. If None, it
                            will be set 1e-11 GeV.
        :param stableWidth: Minimum width for considering particles as stable. If None, it
                            will be set 1e-25 GeV.
        :param roundMasses: If set, it will round the masses to this number of digits (int)

        :param ignorePromptQNumbers: If set, all particles with prompt decays (totalwidth > promptWidth)
                               will have the corresponding properties (quantum numbers). So all particles with the same
                               mass and isSM=False will be considered as equal when combining elements.

        :param minMass: Minimal mass for BSM particles in the model. Any particle with mass below minMass will
                        have its mass set to minMass.

        """

        self.inputFile = inputFile
        if not os.path.isfile(inputFile):
            self.inputFile = "string"
        if promptWidth is None:
            promptWidth = 1e-11*GeV
        if stableWidth is None:
            stableWidth = 1e-25*GeV

        massDict, decaysDict, xsections = self.getModelDataFrom(inputFile)
        self.xsections = xsections

        #  Restric BSM particles to the overlap between allBSMparticles and massDict:
        self.BSMparticles = []
        for particle in self.allBSMparticles:
            if particle.pdg not in massDict and -particle.pdg not in massDict:
                continue
            self.BSMparticles.append(particle.copy())
        if len(self.BSMparticles) == len(self.allBSMparticles):
            logger.info(f"Loaded {len(self.BSMparticles)} BSM particles")
        else:
            logger.info("Loaded %i BSM particles (%i particles not found in %s)"
                        % (len(self.BSMparticles), len(self.allBSMparticles)-len(self.BSMparticles),
                           self.inputFile))

        #  Remove cross-sections for even particles:
        nXsecs = self.filterCrossSections()
        if nXsecs == 0:
            msg = f"No cross-sections found in {self.inputFile} for the BSM particles. "
            msg += "Check if the model is compatible with the input file."
            logger.error(msg)
            raise SModelSError(msg)

        #Set particle masses:
        self.setMasses(massDict,roundMasses,minMass)

        #  Set particle decays
        self.setDecays(decaysDict, promptWidth, stableWidth, ignorePromptQNumbers)

        #  Reset particle equality from all particles:
        for p in Particle.getinstances():
            p._comp = {p._id: 0}
            if isinstance(p, MultiParticle):
                for ptc in p.particles:
                    p._comp[ptc._id] = 0
