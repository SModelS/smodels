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
from smodels.theory.printer import Printer
from smodels.tools.physicsUnits import addunit

class SlhaStatus(Printer):
    """
    An instance of this class represents the status of an SLHA file.
    The output status is 0 if the file is not checked, 1 if the check is ok
    The final status is -1 in case of a physical problem, e.g. charged LSP,
    and -2 in case of formal problems, e.g. missing decay blocks, in the file
    The parameter maxFlightlength is specified in meters. 
    """
    def __init__(self, filename, maxFlightlength=3., maxDisplacement=.1, sigmacut=.01,
                 findMissingDecays=True, findDisplaced=True, massgap=5., maxcond=.2,
                 findEmptyDecays=True, checkXsec=True, checkLSP=True,
                 checkFlightlength=True, model="MSSM"):
        self.filename = filename
        self.model = model
        self.slha = self.read()
        if not self.slha:
            self.status = -3, "Could not read input slha"
            return
        self.maxFlightlength = maxFlightlength
        self.maxDisplacement = maxDisplacement
        self.sigmacut = sigmacut
        self.massgap = massgap
        self.maxcond = maxcond
        self.lsp = self.findLSP()
        self.lspStatus = self.testLSP(checkLSP)
        self.ctauStatus = self.checkCtau(checkFlightlength)
        self.missingDecays = self.checkDecayBlock(findMissingDecays)
        self.emptyDecays = self.findEmptyDecay(findEmptyDecays)
        self.xsec = self.hasXsec(checkXsec)
        self.vertexStatus = self.findDisplacedVertices(findDisplaced)
        self.status = self.evaluateStatus()
    
    def formatData(self):
        return self.formatSLHAData()

    def read(self):
        """
        Get pyslha output object.
        
        """
        try: return pyslha.readSLHAFile(self.filename)
        except: return None


    def evaluateStatus(self):
        """
        Get status summary from all performed checks.

        :returns: a status flag and a message for explanation

        """
        if not self.slha:
            return -3, "Could not read input slha file"
        ret = 0
        retMes = "#Warnings:\n"
        for st, message in [self.missingDecays, self.emptyDecays,
                            self.xsec]:
            if st < 0:
                ret = -2
                retMes = retMes + "#" + message + ".\n"
            elif st == 1 and not ret == -2:
                ret = 1
        for st, message in [self.lspStatus, self.ctauStatus,
                            self.vertexStatus]:
            if st < 0:
                ret = -1
                retMes = retMes + "#" + message + ".\n"
        if ret == 0:
            return 0, "No checks performed"
        if ret == -1:
            return -1, "#ERROR: special signatures in this point.\n" + retMes
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
                 return 1, "Neutral NLSP is long-lived, c*tau = " + str(addunit(ct),"m")
            return -1, "Charged NLSP is long-lived, c*tau = " + str(addunit(ct), "m")
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
        return -1, "XSECTION table missing, will be computed by SModelS"


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
            return pid, addunit(minmass, "GeV")
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
            return pid, addunit(minmass, "GeV")
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

    def findDisplacedVertices(self, findDisplaced):
        """
        find meta-stable particles that decay via leptons or higgs
        """
        if not findDisplaced:
            return 0, "Did not check for displaced vertices"

        particlelist = [1, 2, 3, 4, 5, 6, 11, 13, 15, 25, 35, 36, 37]
        ok = 1
        msg = "Found displaced vertices:\n"
        for particle,block in self.slha.decays.items():
            if particle <= 50: continue
            if particle == self.findLSP(): continue
            if self.slha.blocks["MASS"][particle] - self.slha.blocks["MASS"][self.findLSP()] < self.massgap: continue
            lt = self.getLifetime(particle, ctau=True)
            if lt<self.maxDisplacement: continue
            pcs = 0.
            brvalue = 0.
            old = 0.
            daughters = []
            for decay in block.decays:
                for pid in decay.ids:
                    if abs(pid) in particlelist:
                        if not pcs: pcs = self.getXSEC(particle)
                        brvalue += decay.br
                        daughters.append(decay.ids)
                        break
            if pcs*brvalue*brvalue > self.sigmacut:
                if lt < 1.: msg = msg + "#Displaced vertex: "
                elif lt < 3.: msg = msg + "#Metastable particle: "
                else: msg = msg + "#Longlived particle decay: "
                ok = -1
                msg = msg + "%s is decaying to %s\n" %(particle, str(daughters))
        if ok == 1: msg = "no displaced vertices found"
        return ok, msg

    def getXSEC(self, pid, sqrts=8000.):
        """
        get crosssection for pair production, read last line
        """
        inblock = None
        xsec = 0.
        pid = abs(int(pid))
        for line in open(self.filename):
            line=line.split('#')[0]
            if not line.strip(): continue
            lineelements = line.split()
            if lineelements[0]=="XSECTION":
                if float(lineelements[1])==sqrts:
                    if abs(int(lineelements[-1]))==abs(int(lineelements[-2]))==pid:
                        inblock = True
                        continue
                    else: inblock = None
                else: inblock = None
            if not inblock: continue
            if not line.strip():
                inblock=None
                continue
            xsec += float(line.split()[-3])
        return xsec

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
        exec("from tools.%s_qNumbers import qNumbers" %model)
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
    argparser.add_argument('-maxDisp', '--displacement',
                           help = 'give maximum displacement of secondary vertex in m',
                           default = .001)
    argparser.add_argument('-sigmacut','--sigmacut',
                           help = 'give sigmacut in fb', #FIXME is that true? of pb?
                           default = .01)
    argparser.add_argument('-fD', '--displaced',
                           help = 'find displaced vertices',
                           action = 'store_false')
    argparser.add_argument('-mg', '--massgap',
                           help= 'give massgap for mass compression in GeV',
                           default = 5.)
    args = argparser.parse_args() # pylint: disable-msg=C0103
    status = SlhaStatus(args.filename, args.flightlength, args.displacement, args.sigmacut, args.decays,
                        args.displaced, args.massgap, args.empty, args.xsec, args.lsp, args.ctau,
                        args.model)
    # pylint: disable-msg=C0103
    print(status.status)
