#!/usr/bin/env python3

"""
.. module:: ioObjects
   :synopsis: Definitions of input/output parameters which are read from parameter.in.
    
.. moduleauthor:: Ursula Laa <ursula.laa@lpsc.in2p3.fr>    
.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import os
from smodels.theory import lheReader
from smodels.tools.physicsUnits import GeV, fb
from smodels import installation
import pyslha
from smodels.share.models.SMparticles import SMList, SMparticleList
from smodels.theory.model import Model
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger
from smodels.tools import runtime

SMpdgs = SMparticleList.pdg


class OutputStatus(object):
    """
    Object that holds all status information and has a predefined printout.    
    """
    
    def __init__( self, status, inputFile, parameters, databaseVersion):
        """
        Initialize output. If one of the checks failed, exit.
        
        :parameter status: status of input file
        :parameter inputFile: input file name
        :parameter parameters: input parameters
        :parameter databaseVersion: database version (string)        
        """

        try:
            filename=os.path.join ( installation.installDirectory(),
                                    'smodels/version' )
            with open( filename, 'r') as versionFile:
                version = versionFile.readline()
            self.smodelsVersion = version.replace('\n','')
        except IOError:
            self.smodelsVersion = None

        self.inputfile = inputFile.replace("//","/")
        self.parameters = parameters
        self.filestatus = status[0]
        self.warnings = status[1]
        self.databaseVersion = databaseVersion
        self.statusStrings = {-4: "#database not found",
                              -3: "#no topology passed cut on production cross section",
                              -2: "#bad input file, did not run decomposition",
                              -1: "#could not run the decomposition",
                               0: "#no matching experimental results",
                               1: "#decomposition was successful"}

        self.status = 0
        if not self.databaseVersion:
            self.status = -4
        if self.filestatus < 0:
            self.status = -2


    def updateStatus(self, status):
        """
        Update status.
        
        :parameter status: new status flag
        
        """
        self.status = status

    def updateSLHAStatus(self, status):
        """
        Update SLHA status.
        
        :parameter status: new SLHA status flag
        
        """
        self.slhastatus = status
        return

    def addWarning(self, warning):
        """
        Append warning to warnings.
        
        :parameter warning: warning to be appended
        
        """
        self.warnings += warning
        return



class FileStatus(object):
    """
    Object to run several checks on the input file.
    It holds an LheStatus (SlhaStatus) object if inputType = lhe (slha)
    """

    def __init__(self):

        self.filestatus = None
        self.status = 0, "File not checked\n"


    def checkFile(self, inputFile):
        """
        Run checks on the input file.
        
        :parameter inputFile: path to input file   
        """
        
        inputType = runtime.filetype( inputFile )

        if inputType == 'lhe':
            self.filestatus = LheStatus(inputFile)
            self.status = self.filestatus.status
        elif inputType == 'slha':
            self.filestatus = SlhaStatus(inputFile)
            self.status = self.filestatus.status
        else:
            self.filestatus = None
            self.status = -5, 'Unknown input type: %s' % inputType 


class LheStatus(object):
    """
    Object to check if input lhe file contains errors.
    
    :ivar filename: path to input LHE file
    
    """

    def __init__(self, filename):
        self.filename = filename
        self.status = self.evaluateStatus()

    def evaluateStatus(self):
        """
        run status check
        """
        if not os.path.exists(self.filename):
            # set status flag to -3, as in slha checks for missing input file
            return -3, "Inputfile %s not found" % self.filename
        lhe = lheReader.LheReader(self.filename)
        nevents = lhe.metainfo["nevents"]
        totxsec = lhe.metainfo["totalxsec"]
        sqrts = lhe.metainfo["sqrts"]
        if (not type(sqrts) == type(1 * GeV)) or (not sqrts.asNumber()):
            return -1, "Center-of-mass energy not found in the input LHE file %s" % self.filename
        elif not nevents:
            return -1, "No events found in the input LHE file %s" % self.filename
        elif (not type(totxsec) == type(1 * fb)) or (not totxsec.asNumber()):
            return -1, "Total cross section not found in the input LHE file %s" % self.filename
        return 1, "Input file ok"



class SlhaStatus(object):
    """
    An instance of this class represents the status of an SLHA file.
    The output status is:
    = 0 : the file is not checked,
    = 1: the check is ok
    = -1: case of a physical problem, e.g. charged LSP,
    = -2: case of formal problems, e.g. no cross sections
        
    """
    def __init__(self, filename,
                 findMissingDecayBlocks=True,
                 findIllegalDecays=False, checkXsec=True):
        
        """
        :parameter filename: path to input SLHA file
        :parameter findMissingDecayBlocks: if True add a warning for missing decay blocks
        :parameter findIllegalDecays: if True check if all decays are kinematically allowed
        :parameter checkXsec: if True check if SLHA file contains cross sections
        :parameter findLonglived: if True find stable charged particles and displaced vertices        
        """
        
        self.filename = filename
        self.slha = self.read()
        
        from smodels.particlesLoader import BSMList
        
        if not self.slha:
            self.status = -3, "Could not read input SLHA file"
            return
        try:
            model = Model(BSMList,SMList)
            model.updateParticles(filename)
            self.model = model
            self.illegalDecays = self.findIllegalDecay(findIllegalDecays)
            self.xsec = self.hasXsec(checkXsec)
            self.decayBlocksStatus = self.findMissingDecayBlocks(findMissingDecayBlocks)
            self.status = self.evaluateStatus()
        except (SModelSError,TypeError,IOError,ValueError,AttributeError) as e:
            self.status = -4, "Error checking SLHA file: "+str(e)


    def read(self):
        """
        Get pyslha output object.
        
        """
        try: ret = pyslha.readSLHAFile(self.filename)
        except (pyslha.ParseError,IOError): 
            return None
        if not ret.blocks["MASS"]: return None
        return ret


    def evaluateStatus(self):
        """
        Get status summary from all performed checks.

        :returns: a status flag and a message for explanation

        """
        if not self.slha:
            return -3, "Could not read input slha file"
        ret = 0
        warning = None
        retMes = "#Warnings:\n"
        st , message = self.decayBlocksStatus  # add only warning, no negative staus flag in case of missing decay blocks
        if st < 0:
            warning = True
            retMes += message + "\n"
        for st, message in [self.xsec]:
            if st < 0:
                ret = -2
                retMes = retMes + "#" + message + ".\n"
            elif st == 1 and not ret == -2:
                ret = 1
        for st, message in [self.illegalDecays]:
            if st < 0:
                ret = -1
                retMes = retMes + "#" + message + "\n"
            elif st == 1 and ret >= 0: ret = 1
        if ret == 0:
            return 0, "No checks performed"
        if ret == -1:
            return -1, "#ERROR: special signatures in this point.\n" + retMes
        if ret == -2:
            return -2, retMes
        if not warning: retMes = "Input file ok"
        return ret, retMes

    def emptyDecay(self, pid):
        """
        Check if any decay is missing for the particle with pid
        
        :parameter pid: PID number of particle to be checked
        :returns: True if the decay block is missing or if it is empty, None otherwise
        
        """
        if not abs(pid) in self.slha.decays: return True  # consider missing decay block as empty
        if not self.slha.decays[abs(pid)].decays: return True
        return None


    def findMissingDecayBlocks(self, findMissingBlocks):
        """
        For all non-SMpdgs particles listed in mass block, check if decay block is written
        
        :returns: status flag and message
        
        """
        if not findMissingBlocks:
            return 0, "Did not check for missing decay blocks"
        st = 1
        missing = []
        pids = self.slha.blocks["MASS"].keys()
        for pid in pids:
            if pid in SMpdgs:
                continue
            if not pid in self.slha.decays:
                missing.append(pid)
                st = -1
        if st == 1:
            msg = "No missing decay blocks"
        else: msg = "# Missing decay blocks for %s" % str(missing)
        return st, msg


    def findIllegalDecay(self, findIllegal):
        """
        Find decays for which the sum of daughter masses excels the mother mass
        
        :parameter findIllegal: True if check should be run
        :returns: status flag and message
        
        """
        if not findIllegal:
            return 0, "Did not check for illegal decays"
        st = 1
        badDecay = "Illegal decay for PIDs "
        for particle, block in self.slha.decays.items():
            if particle in SMpdgs : continue
            if not particle in self.slha.blocks["MASS"].keys(): continue
            mMom = abs(self.slha.blocks["MASS"][particle])
            for dcy in block.decays:
                mDau = 0.
                for ptc in dcy.ids:
                    ptc = abs(ptc)
                    if ptc in SMpdgs:
                        smParticle = self.model.getParticlesWith(pdg=ptc)
                        if not smParticle:
                            raise SModelSError("Particle with PDG = %i could not be found." %ptc)
                        elif len(smParticle) != 1:
                            raise SModelSError("Multiple particles defined with PDG = %i in model" %ptc)
                        else:
                            smParticle = smParticle[0]
                        mDau += smParticle.mass/GeV
                    elif ptc in self.slha.blocks["MASS"].keys(): 
                        mDau += abs(self.slha.blocks["MASS"][ptc])
                    else:
                        return -2, "Unknown PID %s in decay of %s" % (str(ptc), str(particle) + ". Add " + str(ptc) + " to smodels/particle.py")
                if mDau > mMom:
                    st = -1
                    if not str(particle) in badDecay: badDecay += str(particle) + " "
        if st == 1:
            badDecay = "No illegal decay blocks"
            
        return st, badDecay


    def hasXsec(self, checkXsec):
        """
        Check if XSECTION table is present in the slha file.
        
        :parameter checkXsec: set True to run the check
        :returns: status flag, message
        
        """
        if not checkXsec:
            return 0, "Did not check for missing XSECTION table"
        with open(self.filename) as f:
            for line in f:
                if "XSECTION" in line:
                    return 1, "XSECTION table present"

        msg = "XSECTION table is missing. Please include the cross section information and try again.\n"
        msg += "\n\t For MSSM models, it is possible to compute the MSSM cross sections"
        msg += " using Pythia through the command:\n\n"
        msg += "\t  ./smodelsTools.py xseccomputer -p -f " + self.filename + " \n\n"
        msg += "\t For more options and information run: ./smodelsTools.py xseccomputer -h\n"
        logger.error(msg)
        return -1, msg
