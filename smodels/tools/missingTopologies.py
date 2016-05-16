#!/usr/bin/env python

"""
.. module:: tools.missingTopologies
   :synopsis: Definitions of classes used to find, format missing topologies
    
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>    
.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>

"""
from smodels.tools.physicsUnits import fb


class MissingTopo():
    """
    Object to describe one missing topology result
    :ivar topo: topology description
    :ivar weights: weights dictionary
    """
    def __init__(self, topo, weights, contributingElements=[]):
        self.topo = topo
        self.weights = weights
        self.value = None
        self.contributingElements = contributingElements

class MissingTopoList(object):
    """
    Object to find and collect MissingTopo objects, plus printout functionality
    
    :ivar sqrts: center of mass energy for which missing topologies should be evaluated
    
    """
    def __init__(self, sqrts):
        self.sqrts = sqrts
        self.topos = []

    def formatData(self,outputLevel):
        return self.formatMissingData(outputLevel)

    def addToTopos(self, el, sumL=None, sumJet=None):
        """
        adds an element to the list of missing topologies
        if the element contributes to a missing topology that is already
        in the list, add weight to topology
        :parameter el: element to be added
        """
        name = self.orderbranches(self.generalName(el.__str__(), sumL, sumJet))
        for topo in self.topos:
            if name == topo.topo:
                topo.weights += el.weight
                topo.contributingElements.append(el)
                return
        self.topos.append(MissingTopo(name, el.weight, [el]))
        return

    def generalName(self, instr, sumL=None, sumJet=None):
        """
        generalize by summing over charges
        e, mu are combined to l
        :parameter instr: element as string
        :returns: string of generalized element
        """
        from smodels.theory.particleNames import ptcDic
        if sumL: exch = ["W", "l", "t", "ta"]
        else: exch = ["W", "e", "mu", "t", "ta"]
        if sumJet: exch.append("jet")
        for pn in exch:
            for on in ptcDic[pn]:
                instr = instr.replace(on, pn).replace("hijetjets","higgs")
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
    
    
    def findMissingTopos(self, smstoplist, listOfAnalyses, minmassgap, doCompress, doInvisible, sumL=None, sumJet=None):
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
        allMothers = []
        for el in smstoplist.getElements():
            for mEl in el.motherElements:
                cID = mEl[-1].elID
                if not cID in allMothers: allMothers.append(cID)
        for el in smstoplist.getElements():
            if el.elID in allMothers: continue
            if el.covered > 0: continue
            self.addToTopos(el, sumL, sumJet)
        for topo in self.topos:
            if not topo.weights.getXsecsFor(self.sqrts): continue
            topo.value = topo.weights.getXsecsFor(self.sqrts)[0].value / fb
        return    

#     def findMissingTopos(self, smstoplist, resultList, sumL=None, sumJet=None):
#         """
#         Loops over all the elements in smstoplist and checks if the elements
#         are tested by any of the analysis in listOfAnalysis.
#         
#         :parameter smstoplist: list of topologies (TopologyLis object)
#         :parameter resultList: a ResultList object
#         :parameter minmassgap: the parameter for mass compression (Unum object)
#         :parameter doCompress: if set to True will ignore elements which can be mass compressed (True/False)
#         :parameter doInvisible: if set to True will ignore elements which can be invisibly compressed (True/False)
#         :parameter sumL: if True, missing topologies will not distinguish e and mu
#         """
# 
#         #First collect all element IDs from decomposition
#         allIDs = []
#         idDict = {}
#         motherIDs = []
#         for el in smstoplist.getElements():
#             idDict[el.elID] = el
#             allIDs.append(el.elID)
#             for mEl in el.motherElements:
#                 cID = mEl[-1].elID
#                 motherIDs.append(cID)
#         motherIDs = set(motherIDs)
#         allIDs = set(allIDs)
#         allIDs = allIDs.difference(motherIDs)  #Remove mother elements to avoid double counting (???)
#         
#         #Now group the results by sqrts
#         sqrtsDict = {}
#         for tp in resultList.theoryPredictions:
#             sqrts = tp.expResult.globalInfo.sqrts
#             if sqrts in sqrtsDict:
#                 sqrtsDict[sqrts].append(tp)
#             else:
#                 sqrtsDict[sqrts] = [tp]
#         
#         for sqrts,tps in sqrtsDict.items():
#             usedIDs = []
#             for tp in tps:
#                 usedIDs += tp.IDs
#             usedIDs = set(usedIDs)
#             missedIDs = allIDs.difference(usedIDs)
#             missedEls = [idDict[elID] for elID in missedIDs]
#             for el in missedEls:
#                 self.addToTopos(el, sumL, sumJet)
#             for topo in self.topos:
#                 if not topo.weights.getXsecsFor(self.sqrts): continue
#                 topo.value = topo.weights.getXsecsFor(self.sqrts)[0].value/fb