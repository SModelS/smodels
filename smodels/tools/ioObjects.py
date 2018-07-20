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
from smodels.theory import crossSection
from smodels.theory.theoryPrediction import TheoryPrediction
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger
from smodels.tools import runtime

class ResultList(object):
    """
    Class that collects a list of theory predictions plus the corresponding upper limits.
    """
    
    
    def __init__(self, theoPredictionsList=[],maxcond = 1.):
        """
        :parameter theoryPredictionList: list of TheoryPrediction objects
        """
        
        self.theoryPredictions = []
        if theoPredictionsList:
            for theoPred in theoPredictionsList:
                self.addTheoPrediction(theoPred,maxcond)
            self.sort()

    def addTheoPrediction(self, theoPred, maxcond):
        """
        Add a result to the theoryPredictions, unless it violates maxcond.
        
        :parameter theoPred: a Theory Prediction object to be added to ResultList
        :parameter maxcond: maximum condition violation
        
        """
        
        if not isinstance(theoPred,TheoryPrediction):
            logger.error("Only TheoryPrediction objects can be added to ResultList")
            raise SModelSError()
                    
        mCond = theoPred.getmaxCondition()
        if mCond == 'N/A' or mCond > maxcond:
            return False
        
        self.theoryPredictions.append(theoPred)
        return True
    
    def getR(self, theoPred, expected = False):
        """
        Calculate R value.
        
        :parameter theoPred: Theory Prediction object
        :returns: R value = weight / upper limit        
        """
        return theoPred.getRValue( expected )

    def _getRNone(self,theoPred, expected = False ):
        """
        Simple helper function to sort also with None values.
        None is replaced with -1.
        """
        ret = self.getR ( theoPred, expected )
        if ret == None: return -1.
        return ret

    def sort(self):
        """
        Reverse sort theoryPredictions by R value.
        
        """
        self.theoryPredictions = sorted( self.theoryPredictions, key=self._getRNone, 
                                         reverse=True )

    def getBestExpected(self):
        """
        Find EM result with the highest expected R vaue.
        :returns: Theory Prediction object
        """
        rexpMax = -1.
        bestExp = None
        for tP in self.theoryPredictions:
            expResult = tP.expResult
            datasetID = tP.dataset.getID()
            dataType = expResult.datasets[0].getType()
            if dataType != 'efficiencyMap':
                continue
            ulExp = expResult.getUpperLimitFor(dataID=datasetID, expected = True)
            rexp = tP.value[0].value/ulExp
            if rexp > rexpMax:
                rexpMax = rexp
                bestExp = tP
        return bestExp


    def isEmpty(self):
        """
        Check if outputarray is empty.
        
        """
        return len(self.theoryPredictions) == 0


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
        self.statusStrings = {-1: "#could not run the decomposition",
                              - 3: "#no cross sections above sigmacut found",
                              - 4: "#database not found",
                              - 2: "#bad input file, did not run decomposition",
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


    def checkFile(self, inputFile, sigmacut=None):
        """
        Run checks on the input file.
        
        :parameter inputFile: path to input file
        :parameter sigmacut: sigmacut in fb        
        """
        
        
        
        inputType = runtime.filetype( inputFile )

        if inputType == 'lhe':
            self.filestatus = LheStatus(inputFile)
            self.status = self.filestatus.status
        elif inputType == 'slha':
            self.filestatus = SlhaStatus(inputFile, sigmacut=sigmacut)
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
    def __init__(self, filename, maxDisplacement=.01, sigmacut=.03 * fb,
                 checkLSP=True, findMissingDecayBlocks=True,
                 findIllegalDecays=False, checkXsec=True, findLonglived=True):
        
        """
        :parameter filename: path to input SLHA file
        :parameter maxDisplacement: maximum c*tau for promt decays in meters
        :parameter sigmacut: sigmacut in fb
        :parameter checkLSP: if True check if LSP is neutral
        :parameter findMissingDecayBlocks: if True add a warning for missing decay blocks
        :parameter findIllegalDecays: if True check if all decays are kinematically allowed
        :parameter checkXsec: if True check if SLHA file contains cross sections
        :parameter findLonglived: if True find stable charged particles and displaced vertices        
        """
        
        self.filename = filename
        self.maxDisplacement = maxDisplacement
        self.sigmacut = sigmacut
        self.slha = self.read()
        if not self.slha:
            self.status = -3, "Could not read input SLHA file"
            return
        try:
            self.lsp = self.findLSP()
            self.lspStatus = self.testLSP(False)
            self.illegalDecays = self.findIllegalDecay(findIllegalDecays)
            self.xsec = self.hasXsec(checkXsec)
            self.decayBlocksStatus = self.findMissingDecayBlocks(findMissingDecayBlocks)
            self.longlived = self.findLonglivedParticles(False)
            self.status = self.evaluateStatus()
        ## except Exception,e:
        except (SModelSError,TypeError,IOError,ValueError,AttributeError) as e:
            self.status = -4, "Error checking SLHA file: "+str(e)


    def read(self):
        """
        Get pyslha output object.
        
        """
        try: ret = pyslha.readSLHAFile(self.filename)
        except (pyslha.AccessError,pyslha.ParseError,IOError): 
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
        for st, message in [self.lspStatus,
                            self.longlived, self.illegalDecays]:
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
        Check if any decay is listed for the particle with pid
        
        :parameter pid: PID number of particle to be checked
        :returns: True if the decay block is missing or if it is empty, None otherwise
        
        """
        if not abs(pid) in self.slha.decays: return True  # consider missing decay block as empty
        if not self.slha.decays[abs(pid)].decays: return True
        return None


    def findMissingDecayBlocks(self, findMissingBlocks):
        """
        For all non-rEven particles listed in mass block, check if decay block is written
        
        :returns: status flag and message
        
        """
        if not findMissingBlocks:
            return 0, "Did not check for missing decay blocks"
        st = 1
        missing = []
        pids = self.slha.blocks["MASS"].keys()
        from smodels.particlesLoader import rEven
        for pid in pids:
            if pid in rEven:
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
        from smodels.particlesLoader import rEven
        for particle, block in self.slha.decays.items():
            if particle in rEven : continue
            if not particle in self.slha.blocks["MASS"].keys(): continue
            mMom = abs(self.slha.blocks["MASS"][particle])
            for dcy in block.decays:
                mDau = 0.
                for ptc in dcy.ids:
                    ptc = abs(ptc)
                    if ptc in SMmasses: mDau += SMmasses[ptc]
                    elif ptc in self.slha.blocks["MASS"].keys(): mDau += abs(self.slha.blocks["MASS"][ptc])
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


    def testLSP(self, checkLSP):
        """
        Check if LSP is charged.
        
        :parameter checkLSP: set True to run the check
        :returns: status flag, message
        
        """
        if not checkLSP:
            return 0, "Did not check for charged lsp"
        qn = Qnumbers(self.lsp)
        if qn.pid == 0:
            return -1, "lsp pid " + str(self.lsp) + " is not known\n"
        if qn.charge3 != 0 or qn.cdim != 1:
            return -1, "lsp has 3*electrical charge = " + str(qn.charge3) + \
                       " and color dimension = " + str(qn.cdim) + "\n"
        return 1, "lsp is neutral"


    def findLSP(self, returnmass=None):
        """
        Find lightest particle (not in rEven).
        
        :returns: pid, mass of the lsp, if returnmass == True
        
        """
        pid = 0
        minmass = None
        from smodels.particlesLoader import rEven
        for particle, mass in self.slha.blocks["MASS"].items():
            if particle in rEven:
                continue
            mass = abs(mass)
            if minmass == None:
                pid, minmass = particle, mass
            if mass < minmass:
                pid, minmass = particle, mass
        if returnmass:
            return pid, minmass * GeV
        return pid


    def getLifetime(self, pid, ctau=False):
        """
        Compute lifetime from decay-width for a particle with pid.
        
        :parameter pid: PID of particle
        :parameter ctau: set True to multiply lifetime by c
        :returns: lifetime
        
        """
        widths = self.getDecayWidths()
        try:
            if widths[abs(pid)]: lt = (1.0 / widths[abs(pid)]) / 1.51926778e24
            else:
                # Particle is stable
                return -1
            if self.emptyDecay(pid): return -1  # if decay block is empty particle is also considered stable
            if ctau:
                return lt * 3E8
            else:
                return lt
        except KeyError:
            logger.warning("No decay block for %s, consider it as a stable particle" % str(pid) )
            return -1


    def sumBR(self, pid):
        """
        Calculate the sum of all branching ratios for particle with pid.
        
        :parameter pid: PID of particle
        :returns: sum of branching ratios as given in the decay table for pid
        
        """
        decaylist = self.slha.decays[pid].decays
        totalBR = 0.
        for entry in decaylist:
            totalBR += entry.br
        return totalBR


    def deltaMass(self, pid1, pid2):
        """
        Calculate mass splitting between particles with pid1 and pid2.
        
        :returns: mass difference
        
        """
        m1 = self.slha.blocks["MASS"][pid1]
        m2 = self.slha.blocks["MASS"][pid2]
        dm = abs(abs(m1) - abs(m2))
        return dm


    def findNLSP(self, returnmass=None):
        """
        Find second lightest particle (not in rEven).
        
        :returns: pid ,mass of the NLSP, if returnmass == True
        
        """
        lsp = self.findLSP()
        pid = 0
        minmass = None
        from smodels.particlesLoader import rEven
        for particle, mass in self.slha.blocks["MASS"].items():
            mass = abs(mass)
            if particle == lsp or particle in rEven:
                continue
            if minmass == None:
                pid, minmass = particle, mass
            if mass < minmass:
                pid, minmass = particle, mass
        if returnmass:
            return pid, minmass * GeV
        return pid


    def getDecayWidths(self):
        """
        Get all decay-widths as a dictionary {pid: width}.
        
        """
        widths = {}
        for particle, block in self.slha.decays.items():
            widths[particle] = block.totalwidth
        return widths


    def getDecayWidth(self, pid):
        """
        Get the decay-width for particle with pid, if it exists.
        
        """
        widths = self.getDecayWidths()
        try:
            return widths[pid]
        except KeyError:
            print("%s is no valid PID" % pid)


    def massDiffLSPandNLSP(self):
        """
        Get the mass difference between the lsp and the nlsp.
        
        """
        lsp = self.findLSP()
        nlsp = self.findNLSP()
        return self.deltaMass(lsp, nlsp)


    def findLonglivedParticles(self, findLonglived):
        """
        Find meta-stable particles that decay to visible particles
        and stable charged particles.
        
        :returns: status flag, message
        
        """
        if not findLonglived:
            return 0, "Did not check for long lived particles"

        # Get list of cross sections:
        xsecList = crossSection.getXsecFromSLHAFile(self.filename)
        # Check if any of particles being produced have visible displaced vertices
        # with a weight > sigmacut
        chargedList = []
        missingList = []
        ltstr = ""
        from smodels.particlesLoader import rEven
        for pid in xsecList.getPIDs():
            if pid in rEven: continue
            if pid == self.findLSP(): continue
            xsecmax = xsecList.getXsecsFor(pid).getMaxXsec()
            if xsecmax < self.sigmacut: continue
            lt = self.getLifetime(pid, ctau=True)
            if lt < 0:
                # error for stable charged particles
                if self.visible(abs(pid)):
                    if not abs(pid) in chargedList:
                        chargedList.append(abs(pid))
                        if not str(abs(pid)) in ltstr: ltstr += "#%s : c*tau = inf\n" % str(abs(pid))
            if lt < self.maxDisplacement: continue
            brvalue = 0.
            daughters = []
            # Sum all BRs which contain at least one visible particle
            for decay in self.slha.decays[abs(pid)].decays:
                for pidb in decay.ids:
                    if self.visible(abs(pidb), decay=True):
                        brvalue += decay.br
                        daughters.append(decay.ids)
                        break
                    elif self.visible(abs(pidb), decay=True) == None:
                        if not abs(pidb) in missingList:
                            missingList.append(abs(pidb))
            if xsecmax * brvalue > self.sigmacut:
                if not abs(pid) in chargedList:
                    chargedList.append(abs(pid))
                    if not str(abs(pid)) in ltstr: ltstr += "#%s : c*tau = %s\n" % (str(abs(pid)), str(lt))
        if not chargedList and not missingList: return 1, "no long lived particles found"
        else:
            msg = ""
            if chargedList:
                msg += "#Visible decays of longlived particles / stable charged particles: %s\n%s" % (str(chargedList), ltstr)
            if missingList:
                msg += "#Missing decay blocks of new r-Even particles appearing in displaced vertices: %s\n" % (str(missingList))
        return -1, msg

    def degenerateChi(self):
        """
        Check if chi01 is lsp and chipm1 is NLSP. If so, check mass splitting.
        This function is not used, the limit is arbitrary.

        """
        lsp, m1 = self.findLSP(returnmass=True)
        nlsp, m2 = self.findNLSP(returnmass=True)
        if lsp == 1000022 and nlsp == 1000024:
            if abs(m2) - abs(m1) < 0.18:
                return True
        return None

    def visible(self, pid, decay=None):
        """
        Check if pid is detectable.
        If pid is not known, consider it as visible.
        If pid not SM particle and decay = True, check if particle or decay products are visible.
        
        """
        if pid in SMvisible: return True
        if pid in SMinvisible: return False
        qn = Qnumbers(pid)
        if qn.pid == 0:
            return True
        if qn.charge3 != 0 or qn.cdim != 1:
            return True
        if decay:
            if not pid in self.slha.decays:
                logger.warning("Missing decay block for pid %s" % (str(pid)))
                return None  # Note: purposely distinguished from False so I can propagate the information to the output file
            for decay in self.slha.decays[pid].decays:
                for pids in decay.ids:
                    if self.visible(abs(pids), decay=True): return True
        return False


class Qnumbers:
    """
    An instance of this class represents quantum numbers.
    
    Get quantum numbers (spin*2, electrical charge*3, color dimension) from qNumbers.
    
    """
    def __init__(self, pid):
        self.pid = pid
        from smodels.particlesLoader import qNumbers
        if not pid in qNumbers.keys():
            self.pid = 0
        else:
            self.l = qNumbers[pid]
            self.spin2 = self.l[0]
            self.charge3 = self.l[1]
            self.cdim = self.l[2]

SMmasses = {1: 4.8e-3 , 2: 2.3e-3 , 3: 95e-2 , 4: 1.275 , 5: 4.18 , 6: 173.21 , 11: 0.51099e-3 , 12: 0, 13: 105.658e-3 , 14: 0 , 15: 1.177682 , 16: 0 , 21: 0 , 22: 0 , 23: 91.1876 , 24: 80.385 , 25: 125.5, 111: 0.135, 211: 0.140}

SMvisible = [1, 2, 3, 4, 5, 6, 11, 13, 15, 21, 22, 23, 24, 25, 211, 111]
SMinvisible = [12, 14, 16]
