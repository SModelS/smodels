""" all data classes necessary to create a SModelS description of events """

from SMSmethods import GTop, eltostr, Reven, EElement, BElement, PtcDic

class EAnalysis:  
  """ an analysis/topology pair """
  def __init__(self):
    self.label = ""
    self.Top = GTop()
    self.sqrts = 0
    self.lum = 0
    self.results = {}
    self.plots = {}
    self.run = ""
    self.masscomp = 0.2


#Given the constraints dictionary, automatically fill the element list with the
#elements corresponding to the strings in the dictionary, skipping repeated ones
  def GenerateElements(self):
     
    ListOfStrs = []
    vertnumb = self.Top.vertnumb
    vertparts = self.Top.vertparts
#Syntax check:    
    for k in range(len(vertnumb)):
      if len(vertparts[k]) != vertnumb[k]:
        print "GenerateElements: Inconsistent data: ninsertions=%d len(insertions)=%d for ``%s''." % ( vertnumb[k], len(vertparts[k]), self.Top )
        return False
    
#Get all element strings:    
    inelements = self.results.items()
    
    for iii in range(len(inelements)):
      for ii in range(2):  
        con = inelements[iii][ii].replace(" ","")
        while "[" in con:  #String has element        
          st = con[con.find("[[["):con.find("]]]")+3] #Get duplet
          con = con.replace(st,"")  # Remove element duplet
          ptclist = eltostr(st)   # Get particle list
#Syntax checks:
          for ib in range(2):
            for ipt in range(len(ptclist[ib])):
              if len(ptclist[ib][ipt]) != vertparts[ib][ipt]:
                print "GenerateElements: Wrong syntax2"
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
      ptclist = eltostr(ListOfStrs.pop()) 
      NewEl = EElement()
      NewEl.B = [BElement(),BElement()]      
      for ib in range(2):
        NewEl.B[ib].particles = ptclist[ib]
      self.Top.ElList.append(NewEl)
  
  

#Check if the plots listed in results exist
  def GetPlots(self,verbose=True):
    from Experiment import SMSResults

    run = self.run    
    if run == "": run = None  #If run has not been defined, use latest
    for res in self.results.keys():
      if not self.plots.has_key(res):
        if verbose: print "SMSmethods.py: GetPlots: Plot for result",res,"in Analysis",self.label,"not found"
        topo = ""
        ana = []
      else:
        topo = self.plots[res][0]
        analyses = self.plots[res][1]
        
      for ana in analyses:
        if not SMSResults.exists(ana,topo,run):
          if verbose: print "SMSmethods.py: GetPlots: Histogram for ",topo," in ",ana," for run ",run," not found"

  #Loop over all elements in SMSTopList and add the weight to the 
  #matching elements in Analysis.
  def add(self,SMSTopList):
    import SMSglobals, copy

    #Get analysis center of mass energy:
    sqrts = self.sqrts
    if type(sqrts) == type(1.) or type(sqrts) == type(1): sqrts = addunit(sqrts,'TeV')
    CMdic = SMSglobals.CMdic

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
              self.Top.ElList[jel].weight[imass] = sumweights([OldEl.weight,neweight])
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

