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
from ClusterTools import DoCluster, GoodMass, MassAvg
from Experiment import LimitGetter, SMSResults
import logging, copy

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
  
  
  def MassDist(self,mass1,mass2):
    """ definition of distance between two mass arrays. The function is defined so it uses the analysis
    experimental limits to define distance """
 

#Get upper bounds for each mass if input are not numbers:
    if type(mass1) != type(1.) or type(mass2) != type(1.):
#If masses differ by more than 100%, do not define distance
      if abs(mass1[0][0]-mass2[0][0])/(mass1[0][0]+mass2[0][0]) > 0.5: return None
      xmass1 = self.MassPosition(mass1)
      xmass2 = self.MassPosition(mass2)
    else:
      xmass1 = mass1
      xmass2 = mass2

    if xmass1 is None or xmass2 is None: return None

    d = 2.*abs(xmass1-xmass2)/(xmass1+xmass2)

    if d < 0.: return None   #Skip masses without an upper limit

    return d


  def MassPosition(self,mass,nounit=True):
    """ gives the mass position in upper limit space, using the analysis experimental limit data.
        If nounit=True, the result is given as number assuming fb units """

    xmass = LimitGetter.GetPlotLimit(mass,self,complain=False)
    if type(xmass) != type(addunit(1.,'pb')): return None
    if nounit: xmass = rmvunit(xmass,'fb')
    return xmass
        
        
