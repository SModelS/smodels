#!/usr/bin/env python

"""
    .. module:: ioObjects
    :synopsis: Definitions of input/output parameters which are read from parameter.in
    
    .. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>    
    .. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
"""

from smodels.theory.printer import Printer


class ExptResults:
    """
    A class to store all relevant information for one result
    """
    def __init__(self, aname, topo, sqrts, cond, tval, exptval, topo_description, mmother, mlsp):
        self.aname = aname
        self.topo = topo
        self.sqrts = sqrts
        self.cond = cond
        self.tval = tval
        self.exptval = exptval
        self.rval = tval/exptval
        self.topo_description = topo_description
        self.mmother = mmother
        self.mlsp = mlsp

class ResultList(Printer):
    """
    Class that collects ExptResults objects and has a predefined printout
    """
    def __init__(self, outputarray = [], bestresult = [], bestresultonly = None, describeTopo = None):
        self.outputarray = outputarray
        self.bestresult = bestresult
        self.bestresultonly = bestresultonly
        self.describeTopo = describeTopo

    def addResult(self,res):
        self.outputarray.append(res)
        return

    def findBest(self):
        best = None
        for res in self.outputarray:
            if not best or best.rval<res.rval:
                best = res
        if best: self.bestresult = [best]
        return

    def formatData(self):
        """
        to access printout format
        """
        return self.formatResultsData()

class OutputStatus(Printer):
    """
    Object that holds all status information and has a predefined printout 
    """
    def __init__(self, status, slhastatus, warnings):
        self.status = status
        self.slhastatus = slhastatus
        self.warnings = warnings
        self.statusStrings = {-1: "#could not run the decomposition",
                              -3: "#no cross sections above sigmacut found",
                              -2: "#bad input slha, did not run decomposition",
                               0: "#no matching experimental results",
                               1: "#decomposition was successful"}

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
            if name == topo.topo: #FIXME need to give correct format of el, plus need general name function!
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
        


