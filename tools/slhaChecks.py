#!/usr/bin/python

"""
.. module:: slhaChecks
   :synopsis: Check SLHA file for integrity.

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>

"""

from __future__ import print_function
import setPath
from smodels.tools import modpyslha as pyslha


class SlhaStatus:
    """
    An instance of this class represents the status of an SLHA file.
    The output status is 0 if the file is not checked, 1 if the check is ok
    The final status is -1 in case of a physical problem, e.g. charged LSP,
    and -2 in case of formal problems, e.g. missing decay blocks, in the file
    The parameter maxFlightlength is specified in meters. 
    """
    def __init__(self, filename, maxFlightlength=3., findMissingDecays=True,
                 findEmptyDecays=True, checkXsec=True, checkLSP=True,
                 checkFlightlength=True, model="MSSM"):
        self.filename = filename
        self.model = model
        self.slha = self.read()
        self.maxFlightlength = maxFlightlength
        self.lsp = self.findLSP()
        self.lspStatus = self.testLSP(checkLSP)
        self.ctauStatus = self.checkCtau(checkFlightlength)
        self.missingDecays = self.checkDecayBlock(findMissingDecays)
        self.emptyDecays = self.findEmptyDecay(findEmptyDecays)
        self.xsec = self.hasXsec(checkXsec)
        self.status = self.evaluateStatus()


    def read(self):
        """
        Get pyslha output object.
        
        """
        return pyslha.readSLHAFile(self.filename)


    def evaluateStatus(self):
        """
        Get status summary from all performed checks.

        :returns: a status flag and a message for explanation

        """
        ret = 0
        retMes = "#Warnings:\n"
        for st, message in [self.missingDecays, self.emptyDecays,
                            self.xsec]:
            if st < 0:
                ret = -2
                retMes = retMes + "#" + message + ".\n"
            elif st == 1 and not ret == -2:
                ret = 1
        for st, message in [self.lspStatus, self.ctauStatus]:
            if st < 0:
                ret = -1
                retMes = retMes + "#" + message + ".\n"
        if ret == 0:
            return 0, "No checks performed"
        if ret == -1:
            return -1, "#ERROR: charged lsp or long-lived NLSP.\n" + retMes
        if ret == -2:
            return -2, retMes
        return ret, "Input file ok" 


    def checkDecayBlock(self, findMissingDecays):
        """
        Check if there is a decay table for each particle with pid > 50 given
        in the mass block.
        
        """
        if not findMissingDecays :
            return 0, "Did not check for missing decay blocks"
        st = 1
        missing = "Missing decay block for PIDs"
        pids = self.slha.blocks["MASS"].keys()
        for pid in pids:
            if pid <= 50:
                continue
            if not pid in self.slha.decays:
                missing = missing + ", " + str(pid)
                st = -1
        if st == 1:
            missing = "No missing decay blocks"
        return st, missing


    def checkCtau(self, checkFlightlength):
        """
        Check if c*tau of NLSP is larger than given maximum.
        Report error for long lived, charged NLSP.

        """
        if not checkFlightlength:
            return 0, "Did not check c*tau of NLSP"
        ct = self.getLifetime(self.findNLSP(), ctau = True)
        qn = Qnumbers(self.findNLSP(), self.model)
        if qn.pid == 0:
            return -1, "NLSP pid " + str(self.lsp) + " is not known"
        if qn.charge3 != 0 or qn.cdim != 1:
            c = "charged"
        else:
            c = "neutral"
        #if NLSP is charged it should decay prompt
        if ct < 0:
            if c == "neutral":
                 return 1, "Neutral NLSP is stable"
            return -1, "Charged NLSP is stable"
        if ct > self.maxFlightlength:
            if c == "neutral":
                 return 1, "Neutral NLSP is long-lived, c*tau = " + str(ct)
            return -1, "Charged NLSP is long-lived, c*tau = " + str(ct)
        return 1, "Prompt NLSP (%s) decay" % c


    def findEmptyDecay(self, findEmpty):
        """
        Find all particles that are considered stable in SModelS.
        
        Stable particle means that no branching ratios are given in the
        decay table of the particle in the slha file.
        
        """
        if not findEmpty:
            return 0, "Did not check for empty decay blocks"
        st = 1
        stableParticles = "Empty decay block for PIDs"
        for particle, block in self.slha.decays.items():
            if particle == self.lsp:
                continue
            if not block.decays:
                stableParticles = stableParticles + ", " + str(particle)
                st = -1
        if st == 1:
            stableParticles = "No empty decay blocks"
        return st, stableParticles


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
        return -1, "XSECTION table missing"


    def testLSP(self, checkLSP):
        """
        Check if lsp is charged.
        
        """
        if not checkLSP:
            return 0, "Did not check for charged lsp"
        qn = Qnumbers(self.lsp, self.model)
        if qn.pid == 0:
            return -1, "lsp pid " + str(self.lsp) + " is not known"
        if qn.charge3 != 0 or qn.cdim != 1:
            return -1, "lsp has 3*electrical charge = " + str(qn.charge3) + \
                       " and color dimension = " + str(qn.cdim)
        return 1, "lsp is neutral"


    def findLSP(self, returnmass=None):
        """
        Find lightest particle with pid > 50.
        
        :returns: pid, mass of the lsp, if returnmass == True
        
        """
        pid = 0
        minmass = None
        for particle, mass in self.slha.blocks["MASS"].items():
            if particle <= 50:
                continue
            mass = abs(mass)
            if minmass == None:
                pid, minmass = particle, mass
            if mass < minmass:
                pid, minmass = particle, mass
        if returnmass:
            return pid, minmass
        return pid


    def getLifetime(self, pid, ctau=False):
        """
        Compute lifetime from decay-width for a particle with pid.
        
        """
        widths = self.getDecayWidths()
        try:
            if widths[pid]: lt = (1.0 / widths[pid]) / 1.51926778e24
            else:
                # Particle is stable
                return -1
            if ctau:
                return lt*3E8
            else:
                return lt
        except KeyError:
            print("%s is no valid PID" % pid)
            

    def sumBR(self, pid):
        """
        Calculate the sum of all branching ratios for particle with pid.
        
        """
        decaylist = self.slha[1][pid].decays
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
        Find second lightest particle with pid >= 50.
        
        :returns: pid ,mass of the NLSP, if returnmass == True
        
        """
        lsp = self.findLSP()
        pid = 0
        minmass = None
        for particle, mass in self.slha.blocks["MASS"].items():
            mass = abs(mass)
            if particle == lsp or particle <= 50:
                continue
            if minmass == None:
                pid, minmass = particle, mass
            if mass < minmass:
                pid, minmass = particle, mass
        if returnmass:
            return pid, minmass
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


class Qnumbers:
    """
    An instance of this class represents quantum numbers.
    
    Get quantum numbers (spin*2, electrical charge*3, color dimension) from "model"_qNumbers.
    
    """
    def __init__(self, pid, model="MSSM"):
        exec("from smodels.tools.%s_qNumbers import qNumbers" %model)
        self.pid = pid
        self.model = model
        if not pid in qNumbers.keys():
            self.pid = 0
        else:
            self.l = qNumbers[pid]
            self.spin2 = self.l[0]
            self.charge3 = self.l[1]
            self.cdim = self.l[2]


if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser() # pylint: disable-msg=C0103
    argparser.add_argument('-f', '--filename',
                           help = 'filename of input slha')
    argparser.add_argument('-maxfl', '--flightlength',
                           help = 'maximum c*tau in m',
                           default = 3.)
    argparser.add_argument('-mD', '--decays',
                           help = 'find missing decay blocks',
                           action = 'store_false')
    argparser.add_argument('-fE', '--empty',
                           help = 'find empty decay blocks',
                           action = 'store_false')
    argparser.add_argument('-xS', '--xsec',
                           help = 'check if file contains xsection block',
                           action = 'store_false')
    argparser.add_argument('-lsp', '--lsp',
                           help = 'check if lsp is neutral and colorless',
                           action = 'store_false')
    argparser.add_argument('-ctau', '--ctau',
                           help = 'check if nlsp has prompt decay',
                           action = 'store_false')
    argparser.add_argument('-m', '--model',
                           help = 'give input model, e.g. MSSM',
                           default = 'MSSM')
    args = argparser.parse_args() # pylint: disable-msg=C0103
    status = SlhaStatus(args.filename, args.flightlength, args.decays,
                        args.empty, args.xsec, args.lsp, args.ctau,
                        args.model)
    # pylint: disable-msg=C0103
    print(status.status)
