#!/usr/bin/env python

"""
.. module:: SMSAnalysis
    :synopsis: ...
    
.. moduleauthor:: someone <email@example.com>
    
"""
    
""" all data classes necessary to create a SModelS description of events """

from SMSDataObjects import GTop, EElement
from ParticleNames import Reven, PtcDic
import ClusterTools, TheoryPrediction

class EAnalysis:  
  """ FIXME currently an analysis is a container for the association between
      a topology and n results. Maybe the structure becomes clearer, if the
      association is made only between a topology and one result? """
  def __init__(self):
    self.label = ""
    self.Top = GTop()
    self.sqrts = 0
    self.lum = 0
    self.results = {} ## pairs of (constraint,condition)
    self.plots = {} ## pairs of (constraint,???)
    self.run = ""
    self.masscomp = 0.2

  def __str__(self):
    return self.label

#Given the constraints dictionary, automatically fill the element list with the
#elements corresponding to the strings in the dictionary, skipping repeated ones
  def generateElements(self):
     
    ListOfStrs = []
    vertnumb = self.Top.vertnumb
    vertparts = self.Top.vertparts
#Syntax check:    
    for k in range(len(vertnumb)):
      if len(vertparts[k]) != vertnumb[k]:
        print "[SMSAnalysis.generateElements]: Inconsistent data: ninsertions=%d len(insertions)=%d for ``%s''." % ( vertnumb[k], len(vertparts[k]), self.Top )
        return False
    
#Get all element strings:    
    inelements = self.results.items()
    
    for iii in range(len(inelements)):
      for ii in range(2):  
        con = inelements[iii][ii].replace(" ","")
        while "[" in con:  #String has element        
          st = con[con.find("[[["):con.find("]]]")+3] #Get duplet
          con = con.replace(st,"")  # Remove element duplet
          element=EElement ( st )
          ptclist = element.allParticles()
          ## ptclist = strtoel(st)   # Get particle list
#Syntax checks:
          for ib in range(2):
            for ipt in range(len(ptclist[ib])):
              if len(ptclist[ib][ipt]) != vertparts[ib][ipt]:
                print "[SMSAnalysis.generateElements]: Wrong syntax2"
                return False
              for ptc in ptclist[ib][ipt]:
                if not ptc in Reven.values() and not PtcDic.has_key(ptc):
                  print "GenerateElements: Unknown particle ",ptc
                  return False    
          ListOfStrs.append(st)
     
#Remove repeated elements:
    ListOfStrs = set(ListOfStrs)
#Now add all elements to element list    
    while len(ListOfStrs) > 0:
      NewEl=EElement ( ListOfStrs.pop() )
      #ptclist = strtoel(ListOfStrs.pop()) 
      #NewEl = EElement()
      #NewEl.B = [BElement(),BElement()]      
      #for ib in range(2):
      #  NewEl.B[ib].particles = ptclist[ib]
      self.Top.ElList.append(NewEl)
  
  

#Check if the plots listed in results exist
  def getPlots(self,verbose=True):
    from Experiment import SMSResults

    run = self.run    
    if run == "": run = None  #If run has not been defined, use latest
    for res in self.results.keys():
      if not self.plots.has_key(res):
        if verbose: print "[SMSAnalysis.getPlots] Plot for result",res,"in Analysis",self.label,"not found"
        topo = ""
        ana = []
      else:
        topo = self.plots[res][0]
        analyses = self.plots[res][1]
        
      for ana in analyses:
        if not SMSResults.exists(ana,topo,run):
          if verbose: print "[SMSAnalysis.getPlots] Histogram for ",topo," in ",ana," for run ",run," not found"

  #Loop over all elements in SMSTopList and add the weight to the 
  #matching elements in Analysis.
  def add(self,SMSTopList):
    import copy

    #Get analysis center of mass energy:
    sqrts = self.sqrts
    if type(sqrts) == type(1.) or type(sqrts) == type(1): sqrts = addunit(sqrts,'TeV')
    CMdic = ClusterTools.CMdic

    for itop in range(len(SMSTopList)):
      NewTop = SMSTopList[itop]  
    #Check if topologies match:
      if not NewTop.isEqual ( self.Top): continue
      
    #Loop over (event) element list:
      for iel in range(len(NewTop.ElList)):
    #Loop over analysis elements:
        for jel in range(len(self.Top.ElList)):
          NewEl = copy.deepcopy(NewTop.ElList[iel])
          neweight = NewEl.weight
    #Remove weights which do not match the analysis center of mass energy      
          if sqrts.asNumber() and len(CMdic) > 0:
            for k in neweight.keys():
              if CMdic[k] != sqrts: neweight.pop(k)
            for k in CMdic.keys():
              if CMdic[k] == sqrts and not neweight.has_key(k):
                neweight.update({k : addunit(0.,'fb')})

          OldEl = copy.deepcopy(self.Top.ElList[jel])
          if not NewEl.isSimilar(OldEl,order=False,igmass=True): continue   #Check if particles match
          
    #If particles match, descend to nested mass list and look for match
          added = False
          for imass in range(len(self.Top.ElList[jel].B[0].masses)):
            OldEl.weight = self.Top.ElList[jel].weight[imass]
            for ib in range(len(NewEl.B)):            
              OldEl.B[ib].masses = copy.deepcopy(self.Top.ElList[jel].B[ib].masses[imass])
              
    #Check if elements match (with identical masses) for any branch ordering
            if NewEl.isSimilar(OldEl,order=False):
              self.Top.ElList[jel].weight[imass] = ClusterTools.sumweights([OldEl.weight,neweight])
              added = True
              break   #To avoid double counting only add event to one mass combination
            
          if not added:
    #Check for both branch orderings, but only add one (if matches) to avoid double counting
            if not NewEl.isSimilar(OldEl,order=True,igmass=True):
              NewEl.B[0] = NewTop.ElList[iel].B[1]
              NewEl.B[1] = NewTop.ElList[iel].B[0]
            for ib in range(len(NewEl.B)):
              self.Top.ElList[jel].B[ib].masses.append(NewEl.B[ib].masses)
            self.Top.ElList[jel].weight.append(neweight)

  def evaluateResult(self,res,uselimits = False):
    """ Main method for evaluating the theoretical predictions to the analysis given the cross-section
       times branching ratio to specific final states. It combines equivalent masses using the Cluster tools
       and calls evaluateCluster to compute the theoretical predictions for each cluster. Returns a list of
       dictionaries with the cluster average mass and theoretical values for the result and the condition(s) in
       the analysis."""
    import copy, ClusterTools
    from ClusterTools import DoCluster, MassDist, GoodMass

    output = []
    if not self.plots.has_key(res) or not self.results.has_key(res):
      print "EvalRes: Wrong analysis input"
      return False

  #Get minimum distance parameter (it can be analysis dependent)
    dmin = self.masscomp
  #List of analyses and anlysis itself (necessary to compute distances)
    analyses = self.plots[res]
    ClusterTools.DistAnalyses = [analyses,self]

  #Create a mass list with all masses appearing in the analysis which have similar branch masses:
    Goodmasses = []
    Top = copy.deepcopy(self.Top)
    for iel in range(len(Top.ElList)):
      for imass in range(len(Top.ElList[iel].B[0].masses)):
        mass = [Top.ElList[iel].B[0].masses[imass],Top.ElList[iel].B[1].masses[imass]]
        gmass = GoodMass(mass,MassDist,dmin)
        if gmass:
           Top.ElList[iel].B[0].masses[imass] = gmass[0]
           Top.ElList[iel].B[1].masses[imass] = gmass[1]
           if not gmass in Goodmasses: Goodmasses.append(gmass)

  #Cluster masses:
    MCluster = DoCluster(Goodmasses,MassDist,dmin)

    if MCluster == None or MCluster == False:
      MCluster = []
      for i in range(len(Goodmasses)): MCluster.append(set([i]))
      print "[SMSAnalysis.py] Cluster failed, using unclustered masses"
      
  #Loop over clusters to evaluate constraints and conditions inside each cluster
    for cluster in MCluster:
  #Get masses in cluster
      masscluster = []
      for ic in cluster: masscluster.append(Goodmasses[ic])
  #Shrink Topology elements to cluster:
      NewTop = Top.clusterTopology(masscluster)

  #Now NewTop contains only elements with a common mass (replaced by the average mass)
  #Evaluate result inside cluster
      result = NewTop.evaluateCluster ( res )
  #Evaluate conditions
      conditions = NewTop.evaluateCluster ( self.results[res] )

  #Save cluster result
      mavg = [NewTop.ElList[0].B[0].masses,NewTop.ElList[0].B[1].masses]

  #Check if average mass is inside the cluster (exp. limit for average mass ~ exp. limit for individual masses):
      davg = -1.
      for mass in masscluster:
        davg = max(davg,MassDist(mass,mavg))
      if davg == -1. or davg > dmin:
        print "EvalRes: Wrong clustering"
        continue


      output.append({'mass' : mavg, 'result' : result, 'conditions' : conditions})

    return TheoryPrediction.TheoryPrediction( output )

  def evaluateResults(self, uselimits = False ):
    """ evaluate all the analysis'es results """
    for res in self.results:
      self.evaluateResult( res, uselimits )
      
  def Split(self):
    """ if the analysis contains more than one result or plot, splits in a list of simple analyses
    with a single result/plot. Returns a list of simple analyses. If the analysis is already
    simple, return the a one element list with itself"""
    import copy
    
    
    SplitList = []
    for key in self.results.keys():
      for plot in self.plots[key][1]:
        NewAnalysis = copy.deepcopy(self)
        NewAnalysis.label = plot + ":" + self.plots[key][0]
        NewAnalysis.results = {key : self.results[key]}
        NewAnalysis.plots = {key : [self.plots[key][0],[plot]]}
        SplitList.append(NewAnalysis)
        
    return SplitList
        
        
