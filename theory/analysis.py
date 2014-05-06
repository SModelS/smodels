"""
.. module:: theory.analysis
   :synopsis: Encapsulates all data types around one result of one analysis,
   i.e. the association with one plot and one reference cross section result.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

from smodels.experiment import limitGetter


class ULanalysis(object):
    """
    Class to store upper limit-type analyses.
    
    Stores the conditions and the elements constrained by the analysis as well
    as basic analysis info.
    
    self.conditions -- List of condition strings
    
    self.constraint -- Constraint string
    
    self.elementsEff -- Dictionary with constrained elements as keys and
    efficiencies as values
    
    """
    def __init__(self):
        self.label = ""
        self.sqrts = 0
        self.lum = 0
        self.run = None
        self.conditions = None
        self.constraint = None
        self.elementsEff = {}

    def __str__(self):
        return self.label

    def getEfficiencyFor(self, element):
        """
        Get (simple) efficiency for element.
        
        Returns zero if element is not constrained by the analysis or the
        element multiplicative factor if it is.
        
        :returns: float -- zero, if element is not found
        
        """
        for el in self.elementsEff:
            if element.particlesMatch(el):
                return self.elementsEff[el]
        return 0.

    def getUpperLimitFor(self, mass):
        """
        Get the experimental upper limit for a specific mass array.
        
        :param mass: mass vector for computing the upper limit
        :returns: experimental upper limit for cross-section times BR
        
        """
        return limitGetter.getPlotLimit(mass, self)


class SRanalysis(object):
    """
    Class to store signal region-type of analyses with efficiency maps.
    
    Stores the basic analysis info and contains a method for obtaining the
    efficiency from the database.
    
    """
    def __init__(self):
        self.label = ""
        self.sqrts = 0
        self.lum = 0
        self.run = None

    def __str__(self):
        return self.label

    def getEfficiencyFor(self, element):
        """ TODO: remove before release?
        Get efficiency from a database (dummy for now).
        
        """
        if not element:
            return False
        return False

    def getLimitFor(self):
        """ TODO: remove before release?
        Get experimental limit for a cross-section in a specific signal region.
        (dummy for now)
        
        :returns: experimental upper limit for cross-section in the signal
        region
        
        """
        return False

