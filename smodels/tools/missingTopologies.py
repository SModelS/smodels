#!/usr/bin/env python

"""
.. module:: tools.missingTopologies
   :synopsis: Definitions of classes used to find, format missing topologies
    
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>    
.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>

"""

from smodels.theory.printer import Printer

class MissingTopo():
    """
    Object to describe one missing topology result
    :ivar topo: topology description
    :ivar weights: weights dictionary
    """
    def __init__(self, topo, weights):
        self.topo = topo
        self.weights = weights
        self.value = None

class MissingTopoList(Printer):
    """
    Object to find and collect MissingTopo objects, plus printout functionality
    :ivar sqrts: center of mass energy for which missing topologies should be evaluated
    """
    def __init__(self, sqrts):
        self.sqrts = sqrts
        self.topos = []

    def formatData(self,outputLevel):
        return self.formatMissingData(outputLevel)

    def addToTopos(self, el, sumL=None):
        """
        adds an element to the list of missing topologies
        if the element contributes to a missing topology that is already
        in the list, add weight to topology
        :parameter el: element to be added
        """
        name = self.orderbranches(self.generalName(el.__str__(), sumL))
        for topo in self.topos:
            if name == topo.topo:
                topo.weights += el.weight
                return
        self.topos.append(MissingTopo(name, el.weight))
        return

    def generalName(self, instr, sumL=None):
        """
        generalize by summing over charges
        e, mu are combined to l
        :parameter instr: element as string
        :returns: string of generalized element
        """
        from smodels.theory.particleNames import ptcDic
        if sumL: exch = ["W", "l", "t", "ta"]
        else: exch = ["W", "e", "mu", "t", "ta"]
        for pn in exch:
            for on in ptcDic[pn]:
                instr = instr.replace(on, pn)
        return instr

    def orderbranches(self, instr):
        """
        unique ordering of branches
        :parameter instr: element as string
        :returns: string of ordered element
        """
        from smodels.theory.element import Element
        li = Element(instr).getParticles()
        for be in li:
            for ve in be:
                ve.sort()
        li.sort()
        return str(li).replace("'", "").replace(" ", "")

    def findMissingTopos(self, smstoplist, listOfAnalyses, minmassgap, doCompress, doInvisible, sumL=None):
        """
        Loops over all the elements in smstoplist and checks if the elements
        are tested by any of the analysis in listOfAnalysis.
        
        :parameter smstoplist: list of topologies (TopologyLis object)
        :parameter listOfAnlysis: a list of ULanalysis objects
        :parameter minmassgap: the parameter for mass compression (Unum object)
        :parameter doCompress: if set to True will ignore elements which can be mass compressed (True/Fals)
        :parameter doInvisible: if set to True will ignore elements which can be invisibly compressed (True/False)
        :parameter sumL: if True, missing topologies will not distinguish e and mu
        """
        
        from smodels.tools.physicsUnits import fb
        for top in smstoplist:
            for el in top.elementList:
                if el.compressElement(doCompress, doInvisible, minmassgap):
                    continue
                covered = None
                for ana in listOfAnalyses:
                    if not ana.getEfficiencyFor(el) == 0:
                        covered = True
                if not covered:
                    self.addToTopos(el, sumL)
        for topo in self.topos:
            if not topo.weights.getXsecsFor(self.sqrts): continue
            topo.value = topo.weights.getXsecsFor(self.sqrts)[0].value / fb
        return


