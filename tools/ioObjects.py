#!/usr/bin/env python

"""
    .. module:: ioObjects
    :synopsis: Definitions of input/output parameters which are read from parameter.in
    
    .. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>    
    .. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
"""

from smodels.theory.printer import Printer


class ResultList(Printer):
    """
    Class that collects experimental constraints and has a predefined printout
    """
    def __init__(self, outputarray = [], bestresultonly = None, describeTopo = None):
        self.outputarray = outputarray
        self.bestresultonly = bestresultonly
        self.describeTopo = describeTopo

    def addResult(self,res):
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
        self.sort()
        return self.outputarray[0]

    def isEmpty(self):
        #check if outputarray is empty
        return len(self.outputarray)==0

    def formatData(self):
        """
        to access printout format
        """
        return self.formatResultsData()

class OutputStatus(Printer):
    """
    Object that holds all status information and has a predefined printout 
    """
    def __init__(self, slhastatus, warnings, databaseVersion):
        self.slhastatus = slhastatus
        self.warnings = warnings
        self.databaseVersion = databaseVersion
        self.statusStrings = {-1: "#could not run the decomposition",
                              -3: "#no cross sections above sigmacut found",
                              -4: "#database not found",
                              -2: "#bad input slha, did not run decomposition",
                               0: "#no matching experimental results",
                               1: "#decomposition was successful"}
        self.status = self.initiateStatus()

    def initiateStatus(self):
        """
        evaluate initial status from slhastatus
        """
        if not self.databaseVersion: return -4
        if self.slhastatus == -1 or self.slhastatus == -3: return -2
        return 0

    def updateStatus(self, status):
        self.status = status
        return

    def formatData(self):
        """
        to access printout format
        """
        return self.formatStatusData()

class MissingTopo():
    """
    Object to describe one missing topology result
    """
    def __init__(self, topo, weights):
        self.topo = topo
        self.weights = weights
        self.value = None

class MissingTopoList(Printer):
    """
    Object to find and collect MissingTopo objects, plus printout functionality
    """
    def __init__(self, sqrts):
        self.sqrts = sqrts
        self.topos = []
        
    def formatData(self):
        return self.formatMissingData()

    def addToTopos(self, el):
        name = self.orderbranches(self.generalName(el.__str__()))
        for topo in self.topos:
            if name == topo.topo:
                topo.weights.__add__(el.weight)
                return
        self.topos.append(MissingTopo(name, el.weight))
        return

    def generalName(self,instr):
        from smodels.theory.particleNames import ptcDic
        exch = ["W", "l", "t", "ta"]
        for pn in exch:
            for on in ptcDic[pn]: instr = instr.replace(on, pn)
        return instr

    def orderbranches(self,instr):
        from smodels.theory.element import Element
        li = Element(instr).getParticles()
        li.sort()
        return str(li).replace("'","").replace(" ","")

    def findMissingTopos(self, smstoplist, listOfAnalyses, sigmacut, minmassgap):
        from smodels.tools.physicsUnits import rmvunit, addunit
        for el in smstoplist:
            for sel in el.elementList:
                 if sel.compressElement(True, True, minmassgap): continue
                 covered = None
                 for ana in listOfAnalyses:
                     if not ana.getEfficiencyFor(sel) == 0: covered = True
                 if not covered: self.addToTopos(sel)
        for topo in self.topos:
            topo.value = rmvunit(topo.weights.getXsecsFor(self.sqrts)[0].value,"fb")
        return
        

