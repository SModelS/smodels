"""
.. module:: experiment.analysisObjects
   :synopsis: Holds the analysis objects (Upper Limit type or Efficiency Map type) and
              basic methods.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""

from smodels.theory.printer import Printer
from smodels.theory.particleNames import elementsInStr
from smodels.theory import element
from numpy import sqrt,inf
from scipy import stats,special,integrate


class Analysis:
    """
    Parent class holding the shared analysis info and methods. Base class for the ULanalysis
    and EManalysis classes.
    :ivar label: Analysis label/name
    :ivar sqrts: Analysis center-of-mass energy
    :ivar lum: Analysis luminosity
    """
 
    def __init__(self,globalInfo,data):   
        self.info = globalInfo
        self.data = data
        self.sqrts = globalInfo.getInfo('sqrts')
        self.lum = globalInfo.getInfo('lumi')
        self.label = globalInfo.getInfo('id')


class ULanalysis(Analysis,Printer):
    """
    Class to store one upper limit-type analysis.    
    Stores the conditions and the elements constrained by the analysis as well
    as basic analysis info.

    :ivar conditions: List of conditions strings    
    :ivar constraint: Constraint string
    :ivar elementsEff: Dictionary with constrained elements as keys and
    efficiencies as values    

    """
    def __init__(self,globalInfo,ULdata,txnameInfo):
        Analysis.__init__(self,globalInfo,ULdata)
        Printer.__init__(self)
        self.txInfo = txnameInfo
        self.label += ':' + txnameInfo.txname
        
        self.conditions = txnameInfo.getInfo('fuzzycondition')
        self.constraint = txnameInfo.getInfo('constraint')
        self.elementsEff = self._generateEffs()
        
    def __str__(self):
        return self.label
    
    def _generateEffs(self):
        """
        Generates a trivial efficiency dictionary {element: eff},
        where eff = 1 for all the elements appearing in the constraints or conditions.
        
        :return: dictionary of Element objects and efficiency values ({Element : eff,..})
        """
        
        #First generate all the elements which appear in the constraint or conditions
        elStrings = elementsInStr(self.constraint)
        if self.conditions: elStrings += elementsInStr(self.conditions)
        elementsEff = {}
        elStrings = set(elStrings)
        for elstr in elStrings:
            el = element.Element(elstr)
        elementsEff[el] = 1.
        return elementsEff
        

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

    def getUpperLimitFor(self, mass):
        """
        Get the experimental upper limit for a specific mass array.
        
        :parameter mass: mass vector for computing the upper limit
        :returns: experimental upper limit for cross-section times BR (float with unit or Unum object)      
        """
        
        return self.data.getULFor(mass)

    
    def formatData(self,outputLevel):
        """
        Select data preparation method through dynamic binding.
        :parameter outputLevel: general control for the output depth to be printed 
                            (0 = no output, 1 = basic output, 2 = detailed output,...
        """
        return Printer.formatULanalysisData(self,outputLevel)    


class EManalysis(Analysis,Printer):
    """
    Class to store a efficiency map-type of analysis.    
    Stores the basic analysis info and contains a method for obtaining the
    efficiency maps from the database.
  
    """
    
    def __init__(self,globalInfo,EMdata):
        Analysis.__init__(self,globalInfo,EMdata)
        Printer.__init__(self)
        
    def getEfficiencyFor(self, element):
        """
        Get efficiency for element.
        
        Returns zero if element is not constrained by the analysis.        
        :returns: float -- zero, if element is not found        
        """
        
        return self.data.getEffFor(element)    


    def getUpperLimitFor(self,dummy=None):
        """ Get experimental limit for the signal cross-section*efficiency in the analysis signal region.
                        
        :returns: 95\% C.L. experimental upper limit for the signal cross-section in the signal
                  region        
        """
                       
        if self.data.obsEvents is None or self.data.expBGerror is None or self.data.xpectedBG is None:
            return False
        
        Nobs =  self.data.obsEvents  #Number of observed events
        Nexp = self.data.expectedBG  #Number of expected BG events
        alpha = 0.05 #95% C.L.
        Nmax = 0.5*stats.chi2.isf(alpha,2*(Nobs+1)) - Nexp  #Upper limit on number of signal events
        maxSignalEvents = Nmax  #DOES NOT INCLUDE SYSTEMATIC UNCERTAINTIES
                
        return maxSignalEvents/self.lum
    
    def getPValue(self,signalxsec):
        """
        Computes the p-value using the signal cross-section (signalxsec) and the systematic
        error in the BG (bgsysError = systematic error/expected BG).
        Assumes a Gaussian distribution for the BG systematical error.
        
        :param signalxsec: signal cross-section*efficiency for the signal region with units (Unum object)
        :returns: the p-value (float)
        """
        
        Nbg = float(self.expectedBG)  #Number of expected BG events
        NbgErr = float(self.expBGerror)  #Systematical uncertainty in the BG
        Nobs = self.obsEvents  #Number of observed events
        Nsig = float(signalxsec*self.lum)  #Number of signal events
        
        #Signal + BG prediction:
        Ntot = Nbg+Nsig
        #Normalization        
        n = (1./2.)*(1. + special.erf(Ntot/(sqrt(2.)*NbgErr)))
        #P-value integrand
        def pint(x):            
            pInt = stats.poisson.cdf(Nobs,x)   #poisson.cdf with mean x (=total number of predicted events distributed according to gaussian)
            pInt *= stats.norm.pdf(x,loc=Ntot,scale=NbgErr)  #systematical error weight
            return pInt
        #P-value integral
        p = n*integrate.quad(pint,0,inf)[0]        
        return p
    
    def formatData(self,outputLevel):
        """
        Select data preparation method through dynamic binding.
        :parameter outputLevel: general control for the output depth to be printed 
                            (0 = no output, 1 = basic output, 2 = detailed output,...
        """
#        return Printer.formatEManalysisData(self,outputLevel)       
    
