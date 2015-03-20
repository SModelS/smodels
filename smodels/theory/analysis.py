"""
.. module:: theory.analysis
   :synopsis: Encapsulates all data types around one result of one analysis,
              i.e. the association with one plot and one 
              reference cross section result.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

from smodels.experiment import limitGetter, smsResults
from smodels.theory.printer import Printer
from smodels.theory.auxiliaryFunctions import _memoize

class ULanalysis(Printer):
    """
    Class to store one upper limit-type analysis.    
    Stores the conditions and the elements constrained by the analysis as well
    as basic analysis info.

    :ivar conditions: List of conditions strings    
    :ivar constraint: Constraint string
    :ivar elementsEff: Dictionary with constrained elements as keys and
       efficiencies as values    
    :ivar label: Analysis label/name
    :ivar sqrts: Analysis center-of-mass energy
    :ivar lum: Analysis luminosity
    
    """
    def __init__(self):
        self.label = ""
        self.sqrts = 0
        self.lum = 0
        self.conditions = None
        self.constraint = None
        self.elementsEff = {}

    def __str__(self):
        return self.label

    def getEfficiencyFor(self, element):
        """
        Get (trivial) efficiency for element.        
        Returns zero if element is not constrained by the analysis or the
        element multiplicative factor if it is.
        
        :returns: 1 if element is in constraint, zero otherwise  
              
        """
        for el in self.elementsEff:
            if element.particlesMatch(el):
                return self.elementsEff[el]
        return 0.

    @_memoize
    def getUpperLimitFor(self, mass):
        """
        Get the experimental upper limit for a specific mass array.
        
        :parameter mass: mass vector for computing the upper limit
        :returns: experimental upper limit for cross-section times BR (float with unit or Unum object)  
            
        """
        
        return limitGetter.getPlotLimit(mass, self)
    
    def getBranchCondition(self):
        """
        Most analyses include assumptions about the masses of the elements
        appearing in their constraints.
        This method returns a string describing this condition
        
        :returns: string describing branch condition (from the branchcondition field)
                  or None if no condition is found
                  
        """
        
        ananame, txname = self.label.split(':')
        
        return smsResults.getBranchCondition(ananame,txname)
    
    def formatData(self,outputLevel):
        """
        Select data preparation method through dynamic binding.
        
        :parameter outputLevel: general control for the output depth to be printed 
                            (0 = no output, 1 = basic output, 2 = detailed output,...
                            
        """
        return Printer.formatULanalysisData(self,outputLevel)    


class EManalysis(Printer):
    """
    Class to store a efficiency map-type of analysis.    
    Stores the basic analysis info and contains a method for obtaining the
    efficiency maps from the database.
    
    :ivar label: Analysis label/name
    :ivar sqrts: Analysis center-of-mass energy
    :ivar lum: Analysis luminosity
    
    """
    def __init__(self):
        self.label = ""
        self.sqrts = 0
        self.lum = 0

    def __str__(self):
        return self.label

    def getEfficiencyFor(self, element):
        """
        Get efficiency for element from the database.        
        Returns zero if a efficiency is not found.
        
        .. warning:: not implemented yet
                
        :returns: efficiency value (float). zero, if element is not found  
              
        """
        if not element:
            return False
        return False

    def getLimitFor(self):
        """
        Get experimental limit for a cross-section.
        
        .. warning:: not implemented yet
        
        :returns: experimental upper limit for cross-section
           (float with unit or Unum object)     
                     
        """
        return False

