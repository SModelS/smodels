#!/usr/bin/env python

"""
    .. module:: ioObjects
    :synopsis: Definitions of input/output parameters which are read from parameter.in
    
    .. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>    
    .. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
"""

import os, sys
from smodels.theory import lheReader
from smodels.theory.printer import Printer
from smodels.tools.physicsUnits import GeV, fb
from smodels.tools import modpyslha as pyslha
from smodels.particles import qNumbers, rEven
from smodels.theory import crossSection
import logging

logger = logging.getLogger(__name__)

class ResultList(Printer):
    """
    Class that collects experimental constraints and has a predefined printout
    """
    def __init__(self, outputarray = [], bestresultonly = None, describeTopo = None):
        self.outputarray = outputarray
        self.bestresultonly = bestresultonly
        self.describeTopo = describeTopo

    def addResult(self, res, maxcond):
        mCond = res.getmaxCondition()
        if mCond == 'N/A': return
        if mCond > maxcond: return
        self.outputarray.append(res)
        return

    def getR(self, res):
        #calculate R value
        return res.value[0].value / res.analysis.getUpperLimitFor(res.mass)

    def sort(self):
        #reverse sort outputarray by R value
        self.outputarray = sorted(self.outputarray, key=self.getR, reverse=True)
        return

    def getBestResult(self):
        #return best result
        self.sort()
        return self.outputarray[0]

    def isEmpty(self):
        #check if outputarray is empty
        return len(self.outputarray)==0

    def formatData(self,outputLevel):
        """
        to access printout format
        :param outputLevel: general control for the output depth to be printed 
                            (0 = no output, 1 = basic output, 2 = detailed output,...
        """
        return self.formatResultsData(outputLevel)

class OutputStatus(Printer):
    """
    Object that holds all status information and has a predefined printout 
    """
    def __init__(self,inputStatus,parameters,databaseVersion,outputfile):
        """
        Initialize output. If one of the checks failed, exit.
        """
        
        self.outputfile = outputfile
        self.inputfile = inputStatus.filestatus.filename
        self.parameters = parameters
        self.filestatus = inputStatus.status[0]
        self.warnings = inputStatus.status[1]
        self.databaseVersion = databaseVersion
        self.statusStrings = {-1: "#could not run the decomposition",
                              -3: "#no cross sections above sigmacut found",
                              -4: "#database not found",
                              -2: "#bad input slha, did not run decomposition",
                               0: "#no matching experimental results",
                               1: "#decomposition was successful"}

        self.status = 0
        if not self.databaseVersion or self.filestatus < 0:
            self.status = min(self.databaseVersion,self.filestatus)
        self.checkStatus()
    
    def checkStatus(self):
        if self.status < 0:
            self.printout("file",self.outputfile)
            sys.exit()


    def updateStatus(self, status):
        self.status = status
        self.checkStatus()

    def updateSLHAStatus(self, status):
        self.slhastatus = status
        return

    def addWarning(self, warning):
        self.warnings += warning
        return

    def formatData(self,outputLevel):
        """
        to access printout format
        :param outputLevel: general control for the output depth to be printed 
                            (0 = no output, 1 = basic output, 2 = detailed output,...
        """
        return self.formatStatusData(outputLevel)


class FileStatus(Printer):
    """
    Object to run several checks on the input file.
    It holds an LheStatus (SlhaStatus) object if inputType = lhe (slha)
    """
    
    def __init__(self, inputType, inputFile, sigmacut=None):

        self.inputType = inputType
                
        if self.inputType == 'lhe':
            self.filestatus = LheStatus(inputFile)            
            self.status = self.filestatus.status
        elif self.inputType == 'slha':
            self.filestatus = SlhaStatus(inputFile,sigmacut=sigmacut)            
            self.status = self.filestatus.status
        else:
            self.filestatus = None
            self.status = -5, 'Unknown input type: %s' % self.inputType
            

class LheStatus(Printer):
    """
    Object to check if input lhe file contains errors
    """

    def __init__(self, filename):
        self.filename = filename
        self.status = self.evaluateStatus()

    def evaluateStatus(self):
        if not os.path.exists(self.filename):
            #set status flag to -3, as in slha checks for missing input file
            return -3, "Inputfile %s not found" %self.filename
        lhe = lheReader.LheReader(self.filename)
        if not lhe.metainfo["nevents"]:
            return -1, "No events found in inputfile %s" %self.filename
        elif not lhe.metainfo["totalxsec"]:
            return -1, "Total cross-section not found in inputfile %s" %self.filename
        elif not lhe.metainfo["sqrts"]:
            return -1, "Center-of-mass energy not found in inputfile %s" %self.filename
        return 1, "Input file ok"
    
    
    
class SlhaStatus(Printer):
    """
    An instance of this class represents the status of an SLHA file.
    The output status is:
    = 0 : the file is not checked,
    = 1: the check is ok
    = -1: case of a physical problem, e.g. charged LSP,
    = -2: case of formal problems, e.g. missing decay blocks, in the file
    The parameter maxDisplacement is specified in meters. 
    """
    def __init__(self, filename, maxDisplacement=.01, sigmacut=.01*fb,
                 checkLSP = True,
                 findIllegalDecays = False, checkXsec = True, findLonglived = True):
        self.filename = filename
        self.maxDisplacement = maxDisplacement
        self.sigmacut = sigmacut
        self.slha = self.read()
        if not self.slha:
            self.status = -3, "Could not read input slha"
            return
        self.lsp = self.findLSP()
        self.lspStatus = self.testLSP(checkLSP)
        self.illegalDecays = self.findIllegalDecay(findIllegalDecays)
        self.xsec = self.hasXsec(checkXsec)
        self.longlived = self.findLonglivedParticles(findLonglived)
        
        self.status = self.evaluateStatus()
    
    
    def read(self):
        """
        Get pyslha output object.
        
        """
        try: ret = pyslha.readSLHAFile(self.filename)
        except: return None
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
        retMes = "#Warnings:\n"
        for st, message in [
                            self.xsec]:
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
        return ret, "Input file ok" 

    def emptyDecay(self, pid):
        """
        Check if any decay is listed for the particle with pid
        returns True if the decay block is missing or if it is empty
        """
        if not abs(pid) in self.slha.decays: return True #consider missing decay block as empty
        if not self.slha.decays[abs(pid)].decays: return True
        return None


    def findIllegalDecay(self, findIllegal):
        """
        Find decays for which the sum of daughter masses excels the mother mass
        """
        if not findIllegal:
            return 0, "Did not check for illegal decays"
        st = 1
        badDecay = "Illegal decay for PIDs "
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
                    else: return -2, "Unknown PID %s in decay of %s" %(str(ptc),str(particle)) # FIXME unknown pid, what to do??
                if mDau > mMom:
                    st = -1
                    if not str(particle) in badDecay: badDecay += str(particle)+ " "
        if st == 1:
            badDecay = "No illegal decay blocks"
        return st, badDecay


    def hasXsec(self, checkXsec):
        """
        Check if XSECTION table is present in the slha file.
        
        """
        if not checkXsec:
            return 0, "Did not check for missing XSECTION table"
        f = open(self.filename)
        for line in f:
            if "XSECTION" in line:
                return 1, "XSECTION table present"
        
        msg = "XSECTION table is missing. Please include the cross-section information and try again.\n"
        msg += "\n\t For MSSM models, it is possible to compute the MSSM cross-sections"
        msg += " using Pythia through the command:\n\n"
        msg += "\t  ./smodels.py xseccomputer "+ self.filename+" -p \n\n"        
        msg += "\t For more options and information run: ./smodels.py xseccomputer -h\n"
        logger.error(msg)
        sys.exit()


    def testLSP(self, checkLSP):
        """
        Check if lsp is charged.
        
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
        for particle, mass in self.slha.blocks["MASS"].items():
            if particle in rEven:
                continue
            mass = abs(mass)
            if minmass == None:
                pid, minmass = particle, mass
            if mass < minmass:
                pid, minmass = particle, mass
        if returnmass:
            return pid, minmass*GeV 
        return pid


    def getLifetime(self, pid, ctau=False):
        """
        Compute lifetime from decay-width for a particle with pid.
        
        """
        widths = self.getDecayWidths()
        try:
            if widths[abs(pid)]: lt = (1.0 / widths[abs(pid)]) / 1.51926778e24
            else:
                # Particle is stable
                return -1
            if self.emptyDecay(pid): return -1 # if decay block is empty particle is also considered stable 
            if ctau:
                return lt*3E8
            else:
                return lt
        except KeyError:
            logger.warning("No decay block for %s, consider it as a stable particle")
            return -1 
            

    def sumBR(self, pid):
        """
        Calculate the sum of all branching ratios for particle with pid.
        
        """
        decaylist = self.slha.decays[pid].decays
        totalBR = 0.
        for entry in decaylist:
            totalBR += entry.br
        return totalBR
    

    def deltaMass(self, pid1, pid2):
        """
        Calculate mass splitting between particles with pid1 and pid2.
        
        """
        m1 = self.slha.blocks["MASS"][pid1]
        m2 = self.slha.blocks["MASS"][pid2]
        dm = abs(abs(m1) - abs(m2))
        return dm
    

    def findNLSP(self, returnmass = None):
        """
        Find second lightest particle (not in rEven).
        
        :returns: pid ,mass of the NLSP, if returnmass == True
        
        """
        lsp = self.findLSP()
        pid = 0
        minmass = None
        for particle, mass in self.slha.blocks["MASS"].items():
            mass = abs(mass)
            if particle == lsp or particle in rEven:
                continue
            if minmass == None:
                pid, minmass = particle, mass
            if mass < minmass:
                pid, minmass = particle, mass
        if returnmass:
            return pid, minmass*GeV
        return pid
    

    def getDecayWidths(self):
        """
        Get all decay-widths as a dictionary {pid: width}.
        
        """
        widths = {}
        for particle,block in self.slha.decays.items():
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
        find meta-stable particles that decay to visible particles
        and stable charged particles
        """
        if not findLonglived:
            return 0, "Did not check for long lived particles"
        
        #Get list of cross-sections:
        xsecList = crossSection.getXsecFromSLHAFile(self.filename)
        #Check if any of particles being produced have visible displaced vertices
        #with a weight > sigmacut
        chargedList = []
        for pid in xsecList.getPIDs():
            if pid in rEven: continue
            if pid == self.findLSP(): continue
            xsecmax = xsecList.getXsecsFor(pid).getMaxXsec()            
            if xsecmax < self.sigmacut: continue           
            lt = self.getLifetime(pid, ctau=True)
            if lt < 0:
                #error for stable charged particles
                if self.visible(abs(pid)):
                    if not abs(pid) in chargedList: chargedList.append(abs(pid))
            if lt < self.maxDisplacement: continue
            brvalue = 0.
            daughters = []
            #Sum all BRs which contain at least one visible particle
            for decay in self.slha.decays[abs(pid)].decays:
                for pidb in decay.ids:
                    if self.visible(abs(pidb), decay=True):
                        brvalue += decay.br
                        daughters.append(decay.ids)
                        break
            if xsecmax*brvalue > self.sigmacut:
                if not abs(pid) in chargedList: chargedList.append(abs(pid))
        if not chargedList: return 1,"no long lived particles found"
        else: return -1, "#Visible decays / stable charged particles: %s" %str(chargedList)


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
        Check if pid is detectable
        If pid is not known, consider it as visible
        If pid not SM particle and decay = True, check if particle or decay products are visible
        """
        if pid in SMvisible: return True
        if pid in SMinvisible: return False
        qn = Qnumbers(pid)
        if qn.pid == 0:
            return True
        if qn.charge3 != 0 or qn.cdim != 1:
            return True
        if decay:
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
        if not pid in qNumbers.keys():
            self.pid = 0
        else:
            self.l = qNumbers[pid]
            self.spin2 = self.l[0]
            self.charge3 = self.l[1]
            self.cdim = self.l[2]

SMmasses = {1: 4.8e-3 , 2: 2.3e-3 , 3: 95e-2 ,4: 1.275 ,5: 4.18 , 6: 173.21 ,11: 0.51099e-3 ,12: 0,13: 105.658e-3 ,14: 0 ,15: 1.177682 , 16: 0 , 21: 0 , 22: 0 ,23: 91.1876 ,24: 80.385 , 25: 125.5}

SMvisible = [1,2,3,4,5,6,11,13,15,21,22,23,24,25]
SMinvisible = [12,14,16]
