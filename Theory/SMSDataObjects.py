""" all data classes necessary to create a SModelS description of events """

from SMSmethods import strtoel, EElement, BElement
import ClusterTools
from ParticleNames import Reven, PtcDic, simParticles
from AuxiliaryFunctions import Ceval

class GTop:
  """ global topology. contains a list of elements, and
    the number of vertices, and .. FIXME andre? """

  def __init__(self):
    self.vertnumb = []
    self.vertparts = []
    self.ElList = []

  def leadingElement ( self ):
    """ often, a topology carries only one element, so
      we have a special accessor for this """
    if len(self.ElList)==0: return None
    return self.ElList[0]

  def elements ( self ):
    return self.ElList

  def __str__(self):
    ret="number of vertices=%s number of vertex particles=%s" % \
        ( self.vertnumb, self.vertparts )
    return ret

  def checkConsistency ( self, verbose=False ):
    """ the number of vertices and insertions per vertex is
      redundant information in a topology, so we can perform
      an internal consistency check """
    for element in self.ElList:
      info=element.getEinfo()
      if self.vertnumb!=info["vertnumb"]:
        if verbose: print "[SMSmethods.py] inconsistent topology!!!"
        return False
      if self.vertparts!=info["vertparts"]:
        if verbose: print "[SMSmethods.py] inconsistent topology!!!"
        return False
    if verbose: print "[SMSmethods.py] topology is consistent."
    return True

  def isEqual ( self, Top2, order=False ):
    """ is this topology equal to Top2?
        Returns true if they have the same number of vertices and particles.
        If order=False and each topology has two branches, ignore branch ordering."""
    if order or len(self.vertnumb) != 2 or len(Top2.vertnumb) != 2:
      if self.vertnumb != Top2.vertnumb: return False
      if self.vertparts != Top2.vertparts: return False
      return True
    else:
      x1 = [self.vertnumb[0],self.vertparts[0]]
      x2 = [self.vertnumb[1],self.vertparts[1]]
      xA = [Top2.vertnumb[0],Top2.vertparts[0]]
      xB = [Top2.vertnumb[1],Top2.vertparts[1]]
      if x1 == xA and x2 == xB: return True
      if x1 == xB and x2 == xA: return True
      return False

#Adds Eelement to ElList
#OBS: NewElement must have the correct branch ordering!
  def addElement(self, NewElement):

#First get global topology info from NewElement:
    Einfo = NewElement.getEinfo()
#Sanity checks:
    if Einfo["vertnumb"] != self.vertnumb or Einfo["vertparts"] != self.vertparts:
      print "[GTop.addElement] wrong element topology"
      return False
#Append element to ElList:
    self.ElList.append(NewElement)
    return True

  def evaluateCluster(self, instr):
    """ Evaluates string expression in instr using the elements and weights
        stored in self. FIXME I dont understand this """
    from Tools.PhysicsUnits import addunit

    outstr = instr.replace(" ","")
  #Get ordered list of elements:
    El = []
    iels = 0
    while "[[[" in outstr:  #String has element
      st = outstr[outstr.find("[[["):outstr.find("]]]")+3] #Get duplet
      ptclist = strtoel(st)   # Get particle list
  #Syntax checks:
      for ib in range(2):
        for ptcL in ptclist[ib]:
          for ptc in ptcL:
            if not ptc in Reven.values() and not PtcDic.has_key(ptc):
              print "EvalRes: Unknown particle ",ptc
              return False
      outstr = outstr.replace(st,"El["+str(iels)+"]")  # Replace element
      El.append(ptclist)   #Store elements
      iels +=1

  #Get list of expressions (separated by commas):
    if ";" in outstr or "Csim" in outstr or "Cgtr" in outstr:
      outstrv = outstr.rsplit(";")
    else:
      outstrv = outstr.rsplit(",")

  #Generate zeroweight entry
    zeroweight = {}
    for wk in self.ElList[0].weight.keys():
      zeroweight[wk]=addunit(0.,'fb')
      # zeroweight.update({wk : addunit(0.,'fb')})

  #Find elements in self corresponding to elements in El and fill Elw with the respective weights:
    Elw = []
    for i in range(len(El)):
      Elw.append(zeroweight)
      for j in range(len(self.ElList)):
        AEl = [self.ElList[j].B[0].particles,self.ElList[j].B[1].particles]
        if simParticles(El[i],AEl,useDict=False):
          Elw[i] = self.ElList[j].weight
          break

  #Evaluate the instr expression (condition or constraint) for each weight entry:
    result = {}
    Els = []
    if len(Elw) > 0:
      for w in Elw[0].keys():
        Els = [weight[w] for weight in Elw]
        eout = [Ceval(x,Els) for x in outstrv]
        if len(eout) == 1: eout = eout[0]
        result[w]=eout ## .update({w : eout})
    else:
      eout = [Ceval(x,Els) for x in outstrv]
      result = eout
    return result

  def massCompressedTopology ( self, mingap ):
    """ if two masses in this topology are degenerate, create
        a compressed copy of this topology """
    import copy
    ETopComp = copy.deepcopy(self)
  #Loop over branches
    for ib in range(len(ETopComp.vertnumb)):
      if ETopComp.vertnumb[ib] < 2: continue
  #Remove all external particles between compressed masses
      for ivert in range(ETopComp.vertnumb[ib]-1):
        massA = ETopComp.ElList[0].B[ib].masses[ivert]
        massB = ETopComp.ElList[0].B[ib].masses[ivert+1]
        if abs(massA-massB) < mingap:
          ETopComp.ElList[0].B[ib].particles[ivert] = []
          ETopComp.vertparts[ib][ivert] = 0

  #Remove all vertices and masses with zero particle emissions:
      while ETopComp.vertparts[ib].count(0) > 1:
        ivert = ETopComp.vertparts[ib].index(0)
        ETopComp.vertnumb[ib] -= 1
        massA = ETopComp.vertparts[ib].pop(ivert)
        massA = ETopComp.ElList[0].B[ib].masses.pop(ivert)
        massA = ETopComp.ElList[0].B[ib].particles.pop(ivert)

    if not ETopComp.isEqual(self):
      return ETopComp
    else:
      return False

  def invisibleCompressedTopology ( self ):
    import copy
    ETopComp = copy.deepcopy(self)
    #Loop over branches
    for ib in range(len(ETopComp.vertnumb)):
      if ETopComp.vertnumb[ib] < 2: continue
      #Remove all external neutrinos
      for ivert in range(ETopComp.vertnumb[ib]):
        if ETopComp.vertparts[ib][ivert] > 0:
          ptcs = ETopComp.ElList[0].B[ib].particles[ivert]
          while ptcs.count('nu') > 0: ptcs.remove('nu')   #Delete neutrinos
          ETopComp.ElList[0].B[ib].particles[ivert] = ptcs
          ETopComp.vertparts[ib][ivert] = len(ptcs)
  #First first non-empty vertex at the end of the branch
      inv  = ETopComp.vertnumb[ib]-1
      while inv > 0 and ETopComp.vertparts[ib][inv-1] == 0: inv -= 1
  #Remove empty vertices at the end of the branch:
      ETopComp.vertnumb[ib] = inv + 1
      ETopComp.vertparts[ib] = self.vertparts[ib][0:inv]
      ETopComp.vertparts[ib].append(0)
      ETopComp.ElList[0].B[ib].particles = self.ElList[0].B[ib].particles[0:inv]
      ETopComp.ElList[0].B[ib].masses = self.ElList[0].B[ib].masses[0:inv+1]


    if not ETopComp.isEqual(self):
      return ETopComp
    else:
      return False

  def clusterTopology ( self, masscluster ):
    """ Given this object and the clustered masses, return a new topology
        with only the elements belonging to the cluster
        (all masses are replaced by the average mass and elements with
        equivalent masses have their weights combined) """
  #Compute average mass in cluster
    mavg = ClusterTools.MassAvg(masscluster,"harmonic")

  #Keep only elements which belong to the cluster
    NewTop = GTop()
    NewTop.vertnumb = self.vertnumb
    NewTop.vertparts = self.vertparts

    for iel in range(len(self.ElList)):
      for imass in range(len(self.ElList[iel].B[0].masses)):
        mass = [self.ElList[iel].B[0].masses[imass],self.ElList[iel].B[1].masses[imass]]
        ptc = [self.ElList[iel].B[0].particles,self.ElList[iel].B[1].particles]
        weight = self.ElList[iel].weight[imass]

  #If mass is in cluster, add element to NewTop:
        if mass in masscluster:
          Elm = EElement()
          Elm.weight = weight
          for ib in range(len(mavg)):
            Elm.B.append(BElement())
            Elm.B[ib].masses = mavg[ib]
            Elm.B[ib].particles = ptc[ib]
          match = False
          for iel2 in range(len(NewTop.ElList)):
            ptcB = [NewTop.ElList[iel2].B[0].particles,NewTop.ElList[iel2].B[1].particles]
            if simParticles(ptcB,ptc,useDict=False):
              match = True
              oldweight = NewTop.ElList[iel2].weight
              NewTop.ElList[iel2].weight = sumweights([oldweight,weight])
              break
          if not match: NewTop.addElement(Elm)
    return NewTop

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
          ptclist = strtoel(st)   # Get particle list
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
      ptclist = strtoel(ListOfStrs.pop()) 
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

  def evaluateResult(self,res,uselimits = False):
    """ Evaluate theoretical predictions for the analysis result and conditions.
        FIXME what does that mean? """
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

    return output

  def evaluateResults(self, uselimits = False ):
    """ evaluate all the analysis'es results """
    for res in self.results:
      self.evaluateResult( res, uselimits )

