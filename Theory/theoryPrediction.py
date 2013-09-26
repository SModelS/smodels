#!/usr/bin/env python

"""
.. module:: TheoryPrediction
        :synopsis: Classes encapsulating the results of the computation of reference cross sections and related methods
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

import copy
import clusterTools
from element import Element
import topology
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class PredictionForAnalysis:
    """
    Holds all the relevant information to compute the theory prediction for a single analysis/plot
    """
    
    def __init__(self,Analysis=None):
        self.analysis = Analysis
        self.analysisTopo = topology.Topology()  #stores the topology and theory elements matching the analysis
        self.analysisEls = None   #stores the original analysis elements (with empty masses and weights)
        self.maxMassDist = 0.2  #Maximum relative mass distance allowed between two masses (for clustering)
        self.massDict = {}  #a dictionary to store the list of cluster masses (keys) and the original masses in each cluster (values)
        
        if Analysis:
            anaElements = getAnalysisElements(Analysis)
            if anaElements:
                self.analysisEls = [anaEl.copy() for anaEl in anaElements]                       
                self.analysisTopo.vertnumb = anaElements[0].getEinfo()['vertnumb']
                self.analysisTopo.vertparts = anaElements[0].getEinfo()['vertparts']                
                #Consistency check:
                self.analysisTopo.ElList = anaElements
                if not self.analysisTopo.checkConsistency():
                    logger.error('[PredictionForAnalysis] : elements for analysis '+Analysis.label+' are inconsistent')
                    return False
                self.analysisTopo.ElList = []
        
        

    def getTheoryXSecs(self,SMSTopList):
        """
        Adds all elements in SMSTopList matching one or more elements in analysisEls
        """                
        for it,topo in enumerate(SMSTopList):            
            if topo != self.analysisTopo: continue  #Skip unmatching topologies
            for ie,element in enumerate(topo.ElList):
                element_b = element.switchBranches()
                for anaElement in self.analysisEls:
                    same_order = element.particlesMatch(anaElement,order=True)      #same branch ordering matches
                    oppos_order = element_b.particlesMatch(anaElement,order=True)   #opposite branch ordering matches                    
                    if not same_order and not oppos_order: continue         #Skip unmatching particles, independent of branch ordering
                    newelement = anaElement.copy()                   #Replace the new element particles by the analysis particles
                    for xsec in element.weight.XSections:                        
                        if xsec.info.sqrts == self.analysis.sqrts:
                            newelement.weight.XSections.append(copy.deepcopy(xsec))  #Restrict cross-sections to analysis sqrts                                            
                    newelement.setMasses(element.getMasses(),same_order,oppos_order)   #Set masses according to correct branch ordering
                    self.analysisTopo.addElement(newelement)
                    
        #If some analysis elements were not found in SMSTopList, add them with empty weights and masses:
#         for anaElement in self.analysisEls:
#             if not anaElement.isInList(self.analysisTopo.ElList,igmass=True,useDict=False):
#                 self.analysisTopo.ElList.append(anaElement.copy())
                    
    def clusterElements(self, keepMassInfo=False):
        """
        Replaces the original elements in analysisTopo by their clustered values.
        If keepMassInfo, saves the original masses and their cluster value in massDict
        """
       
        #Delete elements with bad masses and replace the masses by their 'good' value:
        self.setGoodMasses()
        goodMasses = self.getMassList()              
        
        #Copy topology with good masses before clustering:
        oldElements = [el.copy() for el in self.analysisTopo.ElList]
        
        #Cluster masses:
        clusterTools.clusterAnalysis = self.analysis  #set analysis to be used for clustering
        clusters = clusterTools.doCluster(goodMasses,self.maxMassDist)
                
        self.analysisTopo.ElList = []  #remove all the original elements
        for cluster in clusters:
            massCluster = []
            for ic in cluster: massCluster.append(goodMasses[ic])
            avgMass = clusterTools.massAvg(massCluster,"harmonic")
            newElements = self.getClusterElements(avgMass,massCluster,oldElements)
            self.analysisTopo.ElList.extend(newElements)
            if keepMassInfo: self.massDict[avgMass] = massCluster
            
        
    def getClusterElements(self,avgMass,massCluster,elementList):
        """
        Replaces all masses elements in elementList by the cluster mass. If the original element mass
        does not belong to the cluster, give it zero cross-sections
        """
        
        newList = []
        #Generate list of analysis elements with the cluster mass and no weights
        for el in self.analysisEls:
            newEl = el.copy()
            newEl.setMasses(avgMass)
            newList.append(newEl)
        for element in elementList:
            if not element.getMasses() in massCluster: continue   #skip elements outside cluster
            for newEl in newList:
                if element.getParticles() == newEl.getParticles():
                    newEl.weight.combineWith(element.weight)       #Combine cross-sections
                                
        return newList
    
    
    def getMassList(self):
        """ Gets the list of masses appearing in analysisTopo.ElList"""
        massList = []
        for el in self.analysisTopo.ElList:
            if not el.getMasses() in massList: massList.append(el.getMasses())
        return massList

    def setGoodMasses(self):
        """
        Check if the masses appearing in the elements are 'good'. If they are, replace the element mass
        by its good value, otherwise remove the element.
        """
        
        dmin = self.maxMassDist
        for iel,element in enumerate(self.analysisTopo.ElList):
            mass = element.getMasses()
            goodmass = clusterTools.goodMass(mass,dmin)            
            if goodmass:
                element.setMasses(goodmass)
                if not goodmass in goodMasses: goodMasses.append(goodmass)
            else:
                self.analysisTopo.ElList.pop(iel)  #Remove elements with bad mass


def getAnalysisElements(ana):
    """ Generates a list of all elements appearing in the result and conditions of Analysis with mass arrays = None
    and empty cross-sections"""
                 

    stringList = []
    if len(ana.listOfPlots) > 1:
        logger.warning('[getAnalysisElements]: Analysis '+ana.label+' has more than one plot. Using the first one')
     
    plot = ana.listOfPlots[0]    
    
#Get all element strings:
    allElements = plot.result
    if plot.condition: allElements += ';'+plot.condition      
    while "[" in allElements:        #String still contains an element                                
        elementStr = allElements[allElements.find("[[["):allElements.find("]]]")+3] #Get element string
        allElements = allElements.replace(elementStr,"")    # Remove string
        stringList.append(elementStr)
            
#Remove repeated elements:
    stringList = set(stringList)
    
#Now add all elements to the element list
    elementList = []
    for elStr in stringList:        
        newElement = Element(elStr)
        for branch in newElement.branches:
            branch.masses.append(None)
            for vertex in branch.particles:  branch.masses.append(None)
        if newElement.checkConsistency(): elementList.append(newElement)

    return elementList