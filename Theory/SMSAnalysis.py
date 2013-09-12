#!/usr/bin/env python

"""
.. module:: SMSAnalysis
    :synopsis: Encapsulates all data around one result of one analysis, i.e.
      the association with one plot, one reference cross section result, etc
      FIXME is this true?
    
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
    
"""
    
from SMSDataObjects import ATop, AElement, CTop
from ParticleNames import Reven, PtcDic
from Tools.PhysicsUnits import addunit, rmvunit

class EAnalysis:  
  def __init__(self):
    self.label = ""
    self.Top = ATop()
    self.sqrts = 0
    self.lum = 0
    self.results = {} ## pairs of (constraint,condition)
    self.plots = {} ## pairs of (constraint,plots/analyses)
    self.run = ""
    self.masscomp = 0.2 ## maximum allowed relative difference in upper limits for considering two mass arrays similar
    self.ResultList = [] ## a list of result objects (TheoryPrediction.ClusterOutput)

  def __str__(self):
    return self.label


  def generateElements(self):
    """ Given the results dictionary, fills the analysis topology list with the \
    analysis elements corresponding to the strings in the dictionary, skipping repeated ones \
    Information about theoretical cross-sections can be passed through CrossSection.XSectionInfo (CrossSection.XSecInfoList object).\
    If CrossSection.XSectionInfo is not defined default values will be used and stored in CrossSection.XSectionInfo."""
     
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
          element=AElement(st)
          ptclist = element.getParticleList()

#Syntax checks:
          for ib,branch in enumerate(ptclist):
            for iv,vertex in enumerate(branch):
              if len(vertex) != vertparts[ib][iv]:
                print "[SMSAnalysis.generateElements]: Wrong syntax2"
                return False
              for ptc in vertex:
                if not ptc in Reven.values() and not PtcDic.has_key(ptc):
                  print "GenerateElements: Unknown particle ",ptc
                  return False    
          ListOfStrs.append(st)
     
#Remove repeated elements:
    ListOfStrs = set(ListOfStrs)
#Now add all elements to the element list with zero weight
    while len(ListOfStrs) > 0:
      NewEl=AElement(PartStr=ListOfStrs.pop())
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

  #Loop over all topologies in SMSTopList and add their elements to the matching topology in Analysis.
  def add(self,SMSTopList):
    """ Main method for matching decomposition generated topologies to the analysis topology \
      Add the elements in the topology list to the analysis elements\
      :returns: True if successful
    """
    if SMSTopList==False:
      return False
    for NewTop in SMSTopList:
    #Check if topologies match:
      if not NewTop == self.Top: continue   
    #Loop over (event) element list:
      for NewElement in NewTop.ElList: self.Top.addEventElement(NewElement,self.sqrts)

    return True

  def computeTheoryPredictions(self):
    """ Main method for evaluating the theoretical predictions to the analysis \
      given the cross-section times branching ratio to specific final states. \
      It combines equivalent masses using the Cluster tools and calls \
      evaluateCluster to compute the theoretical predictions for each cluster. \
      Produces self.ResultsList.

      :returns: True if successful, None if list is empty.
    """
    from ClusterTools import DoCluster, GoodMass, MassAvg
    from Experiment import LimitGetter

    dmin = self.masscomp
  #Create a mass list with all masses appearing in the analysis elements which have similar branch masses:
    Goodmasses = []
    for El in self.Top.ElList:
      for massweight in El.MassWeightList:
        gmass = GoodMass(massweight.mass,self.MassDist,dmin)
#        print massweight.mass
#        print gmass,dmin,"\n"
        if gmass:
           massweight.mass = gmass
           if not gmass in Goodmasses: Goodmasses.append(gmass)

  #Cluster masses:
    MCluster = DoCluster(Goodmasses,self.MassDist,dmin,MassAvg)

    if MCluster == None or MCluster == False:
      MCluster = []
      for i in range(len(Goodmasses)): MCluster.append(set([i]))
      print "[SMSAnalysis.py] Cluster failed, using unclustered masses"
      
  #Loop over clusters to evaluate constraints and conditions inside each cluster
    self.ResultList = []  # Clear out results
    for cluster in MCluster:
  #Get masses in cluster
      masscluster = []
      for ic in cluster: masscluster.append(Goodmasses[ic])
  #Create clustered topology from the analysis topology elements belonging to cluster:
      NewTop = CTop(self.Top,masscluster)

  #Now NewTop contains only elements with a common mass (replaced by the average mass)
  #Evaluate theoretical predictions for analysis inside cluster
      ClusterResult = NewTop.evaluateCluster( self.results )
      ClusterResult.mass = NewTop.clustermass
      ClusterResult.explimit = LimitGetter.GetPlotLimit(ClusterResult.mass,self,complain=False)


  #Check if average mass is inside the cluster (exp. limit for average mass ~ exp. limit for individual masses):
      mavg = NewTop.clustermass
      davg = -1.
      for mass in masscluster:
        davg = max(davg,self.MassDist(mass,mavg))
      if davg == -1. or davg > dmin:
        import logging
        log = logging.getLogger(__name__)
        log.warning ( "Upper limit for average mass is not similar to upper limit for constituent masses (davg=%s, dmin=%s). Wont cluster." % ( davg, dmin ) )
        continue

      self.ResultList.append(ClusterResult)

    if not self.ResultList: return None      
    return True

      
  def split(self):
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
  
  
  def MassDist(self,mass1,mass2):
    """ definition of distance between two mass arrays. The function is defined so it uses the analysis
    experimental limits to define distance """
    from  Experiment import LimitGetter

  
  #If masses differ by more than 100%, do not define distance
    if abs(mass1[0][0]-mass2[0][0])/(mass1[0][0]+mass2[0][0]) > 0.5:
      return None

#Get upper bounds for each mass:
    xmass1 = LimitGetter.GetPlotLimit(mass1,self,complain=False)
    xmass2 = LimitGetter.GetPlotLimit(mass2,self,complain=False)

    if type(xmass1) != type(addunit(1.,'pb')) and (xmass1==None or xmass1==False):
      #print "[SMSAnalysis.MassDist] no limit for",self.label,"plot 1, masses=",mass1
      return None
    if type(xmass2) != type(addunit(1.,'pb')) and (xmass2==None or xmass2==False):
      #print "[SMSAnalysis.MassDist] no limit for",self.label,"plot 2, masses=",mass2
      return None

    x1 = rmvunit(xmass1,'fb')
    x2 = rmvunit(xmass2,'fb')
    d = 2.*abs(x1-x2)/(x1+x2)

    if d < 0.: return None   #Skip masses without an upper limit

    return d

        
        
