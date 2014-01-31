#!/usr/bin/env python

"""
.. module:: analysis
        :synopsis: Encapsulates all data around one result of one analysis, i.e.
            the association with one plot, one reference cross section result,
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

class ULanalysis(object):
    """Class to store upper limit-type analyses. Stores the conditions and the elements constrained by the analysis
    as well as basic analysis info """    
    def __init__(self):
        self.label = ""
        self.sqrts = 0
        self.lum = 0
        self.run = None
        self.conditions = None    #List of condition strings
        self.constraint = None    #Constraint string
        self.elementsEff = {}     #Dictionary with constrained elements as keys and efficiencies as values
       
    def __str__(self):
        return self.label

    def getEfficiencyFor(self,element):
        """Get (simple) efficiency for element. Equals zero if element is not constrained by the analysis or
        the element multiplicative factor if it is."""            
        for el in self.elementsEff:
            if element.particlesMatch(el): return self.elementsEff[el]
        return 0.       #Return zero, if element is not found
    

class SRanalysis(object):
    """Class to store signal region-type of analyses with efficiency maps.
    Stores the basic analysis info and contains a method for obtaining the efficiency from the database"""    
    def __init__(self):
        self.label = ""
        self.sqrts = 0
        self.lum = 0
        self.run = None
       
    def __str__(self):
        return self.label

    def getEfficiencyFor(self,element):
        """ function to be used to get the efficiency from a database (dummy for now) """
        if not element: return False      
        return False
