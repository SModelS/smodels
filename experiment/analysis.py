#!/usr/bin/env python

"""
.. module:: SMSAnalysis
        :synopsis: Encapsulates the basic data around one result of one analysis, i.e.
            the association with one plot, one reference cross section result, etc
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
        
"""
        
from theory.ParticleNames import Reven, PtcDic
from tools.PhysicsUnits import addunit, rmvunit
import SMSInterpolation
import logging, copy


class AnalysisPlot:
    """
    Stores the basic result for a single analysis result (plot): constraint and conditions and plot name
    """        
    def __init__(self):           
        self.result = None ## result (constraint) in string format
        self.condition = None    ## conditions for the respective result (in string format)
        self.label = None  ## plot label for the corresponding result


class Analysis:
    """
    Stores global information about the analysis and the list of AnalysisPlots
    """ 
    def __init__(self):
        self.label = ""
        self.sqrts = 0
        self.lum = 0
        self.listOfPlots = [] ## holds a list of AnalysisPlot objects (usually with a single entry after split)
        self.run = ""        

    def __str__(self):
        return self.label
          
    def split(self):
        """ if the analysis contains more than one AnalysisPlot, splits in a list of simple analyses
        with a single plot each. Returns a list of simple analyses. If the analysis is already
        simple, return itself"""
        
        if len(self.listOfPlots) == 1: return self
        
        singleAnalyses = []
        for analysisPlot in self.listOfPlots:            
            newAnalysis = Analysis()
            newAnalysis.label = analysisPlot.label + ":" + self.label
            newAnalysis.sqrts = self.sqrts
            newAnalysis.lum = self.lum
            newAnalysis.run = self.run
            newAnalysis.listOfPlots = [analysisPlot]
            singleAnalyses.append(newAnalysis)
                
        return singleAnalyses
    
    
    def getPlotLimit(self,inmass,plot=None,complain = False):
        """ Get upper limit on sigma*BR for a specific array of masses from plot in listOfPlots.
        If plot = None, get the limit from first plot in listOfPlots.
        inmass: array of masses
        plot: AnalysisPlot object"""
  

        if plot is None:  plot = self.listOfPlots[0]
        massarray = copy.deepcopy(inmass)
#Skip empty mass arrays:
        if len(massarray) < 1: 
            if complain: print "[LimitGetter.py] length of massarray < 1"
            return False
#Make sure the two branches have equal masses:
        if massarray[0] != massarray[1]:
            if complain: print "[LimitGetter.py] masses differ between branches"
            return False
        masslist = [rmvunit(mass,'GeV') for mass in massarray[0]]
#Run label:
        run = self.run
        if run == "": run = None   #If run has not been defined, use latest run
        CMSlabel = plot.label   #CMS-type label
        ana = self.label  #analysis name
        limit = SMSInterpolation.UpperLimit(ana, CMSlabel, masslist,debug=complain,run=run)
 
        return limit