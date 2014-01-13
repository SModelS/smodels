#!/usr/bin/env python

"""
.. module:: analysis
        :synopsis: Encapsulates all data around one result of one analysis, i.e.
            the association with one plot, one reference cross section result,
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""
        
import logging

class ULanalysis:
    """Class to store upper limit-type analyses. Stores the conditions and the elements constrained by the analysis
    as well as basic analysis info """    
    def __init__(self):
        self.label = ""
        self.sqrts = 0
        self.lum = 0
        self.conditions = None
        self.constrainedEls = {}  #Dictionary with elements as keys and efficiencies as values
       
    def __str__(self):
        return self.label

    def getEfficiencyFor(self,element):
        """Get (simple) efficiency for element. Equals zero if element is not constrained by the analysis or
        the element multiplicative factor if it is."""            
        for el in self.constrainedEls:
            if element.particlesMatch(el): return self.constrainedEls[el]
        return 0.       #Return zero, if element is not found
        
                
    def split(self):
        """ if the analysis contains more than one result or plot, splits in a list of simple analyses
        with a single result/plot. Returns a list of simple analyses. If the analysis is already
        simple, return the a one element list with itself"""
        
        SplitList = []
        for key in self.results.keys():
            for plot in self.plots[key][1]:
                NewAnalysis = copy.deepcopy(self)
                NewAnalysis.label = plot + ":" + self.plots[key][0]
                NewAnalysis.results = {key : self.results[key]}
                NewAnalysis.plots = {key : [self.plots[key][0],[plot]]}
                SplitList.append(NewAnalysis)
                
        return SplitList
    

class SRanalysis:
    """Class to store signal region-type of analyses with efficiency maps.
    Stores the basic analysis info and contains a method for obtaining the efficiency from the database"""    
    def __init__(self):
        self.label = ""
        self.sqrts = 0
        self.lum = 0
       
    def __str__(self):
        return self.label

    def getEfficiencyFor(self,element):
        """ function to be used to get the efficiency from a database (dummy for now) """      
        return False
