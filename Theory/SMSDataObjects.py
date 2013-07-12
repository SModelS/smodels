#!/usr/bin/env python

"""
.. module:: SMSDataObjects
   :synopsis: Our basic data objects: TopologyList, GTop, EElement, BElement
    
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

import ClusterTools
from ParticleNames import Reven, PtcDic, simParticles
from AuxiliaryFunctions import Ceval
from Tools.PhysicsUnits import addunit, rmvunit

class BElement:
  """ A branch-element """

  def __init__( self, S=None ):
    """ A branch-element can be constructed from a string S (e.g. ('[b,b],[W]')"""
    self.masses = []
    self.particles = []
    self.momID = 0
    if type(S)==type(""):
      st = S.replace(" ","")
      while "[" in st or "]" in st:
        ptcs = st[st.find("[")+1:st.find("],[")].split(",")
        for ptc in ptcs:
          if not ptc in Reven.values() and not PtcDic.has_key(ptc):
            print "[BElement] unknown particle:",ptc
            return
        spcts=str(ptcs).replace("'","").replace(" ","")
        st=st.replace(spcts,"",1)
        self.particles.append ( ptcs )


  def isEqual ( ElA, ElB, order=True ):
    if not simParticles(ElA.particles,ElB.particles,useDict=False): return False
    if ElA.masses != ElB.masses: return False
    return True

  def __eq__ ( self, other ):
    return self.isEqual ( other )

  def isSimilar ( self, elB, order=True, igmass=False ):
    """ compare elB with self.
        If particles are similar and all masses equal, returns True,
        otherwise returns False.
        If order = False, test both branch orderings (for an element doublet only)
        If igmass = True, only compare particles """
    if type (elB) != type(self): return False
    if not simParticles(self.particles,elB.particles): return False
    if not igmass and self.masses != elB.masses: return False
    return True

  def __str__ ( self ):
    """ the canonical SModels description of the BElement. """
    st = str(self.particles).replace("'","")
    st = st.replace(" ","")
    return st

  def describe ( self ):
    """ a lengthy description of the BElement """
    from Tools.PhysicsUnits import rmvunit
    ret="particles=%s masses=%s" % \
       ( self.particles, [ rmvunit(x,"GeV") for x in self.masses ] )
    return ret

class EElement:
  """ An Event Element, contains of several branches and weight information """
  def __init__( self, S=None ):
    """ If S != None, an Element is created from a string description """
    self.B = []
    self.weight = []
    if S:
      st = S.replace(" ","").replace("'","")
      st = st[st.find("[[["):st.find("]]]")+3]
      b1=st[2:st.find("]],[[")+1]
      b2=st[st.find("]],[[")+4:st.find("]]]")+1]
      import copy
      self.B.append ( copy.deepcopy( BElement ( b1 ) ) )
      self.B.append ( copy.deepcopy ( BElement ( b2 ) ) )

#Get global topology info from element structure
  def getEinfo( self ):
    vertnumb = []
    vertparts = []
    for el in self.B:
      vertnumb.append(len(el.masses))
      vertparts.append([len(x) for x in el.particles])
      if len(vertparts[len(vertparts)-1]) == vertnumb[len(vertnumb)-1]-1:
        vertparts[len(vertparts)-1].append(0)  #Append 0 for stable LSP
    return {"vertnumb" : vertnumb, "vertparts" : vertparts}

  def allParticles ( self ):
    """ returns all particles from all branches """
    ret=[]
    for b in self.B:
      ret.append ( b.particles )
    return ret

  def isSimilar ( ElA, ElB,order=True,igmass=False ):
    """ Compare two EElements
        If particles are similar and all masses equal, returns True,
        otherwise returns False
        If order = False, test both branch orderings (for an element doublet only)
        If igmass = True, only compare particles """
    if type(ElA) != type(ElB): return False
    El1 = ElA.B
    El2 = ElB.B
    if len(El1) == 2 and not order:
      ptcsA = [El2[0].particles,El2[1].particles]
      massA = [El2[0].masses,El2[1].masses]
      ptcs = [El1[0].particles,El1[1].particles]
      mass = [El1[0].masses,El1[1].masses]
      ptcs_b = [El1[1].particles,El1[0].particles]
      mass_b = [El1[1].masses,El1[0].masses]
      if igmass:
        mass = massA
        mass_b = massA
      if simParticles(ptcsA,ptcs) and mass == massA:
        return True
      elif simParticles(ptcsA,ptcs_b) and mass_b == massA:
        return True
      else:
        return False
    else:
      for i in range(len(El1)):
        if not simParticles(El1[i].particles,El2[i].particles): return False
        if not igmass and El1[i].masses != El2[i].masses: return False
    return True

  def isEqual ( ElA, ElB,order=True):
    """ Compare two EElements
        If all masses and particles are equal, returns True,
        otherwise returns False
        If order = False, test both branch orderings (for an element doublet only) """
    if type(ElA) != type(ElB): return False
    El1 = ElA.B
    El2 = ElB.B
    if len(El1) == 2 and not order:
      ptcsA = [El2[0].particles,El2[1].particles]
      massA = [El2[0].masses,El2[1].masses]
      ptcs = [El1[0].particles,El1[1].particles]
      mass = [El1[0].masses,El1[1].masses]
      ptcs_b = [El1[1].particles,El1[0].particles]
      mass_b = [El1[1].masses,El1[0].masses]
      if simParticles(ptcsA,ptcs,useDict=False) and mass == massA:
        return True
      elif simParticles(ptcsA,ptcs_b,useDict=False) and mass_b == massA:
        return True
      else:
        return False
    else:
      for i in range(len(El1)):
        if not simParticles(El1[i].particles,El2[i].particles,useDict=False): return False
        if El1[i].masses != El2[i].masses: return False

    return True

  def __eq__ ( self, other ):
    return self.isEqual ( other )

  def __str__ ( self ):
    """ returns the canonical name of the element, e.g. [[jet],[jet]] """
    ret="["
    for i in self.B:
      ret+=str(i)+","
    if len(ret)>1:
      ret=ret[:-1]
    ret+="]"
    return ret

  def describe ( self ):
    """ returns a lengthy description of the event element """
    ret="Branch #1={{"+str(self.B[0])+"}}, Branch #2={{"+str(self.B[1])+"}}"
    return ret

class AElement:
  """ An Analysis Element, contains a string with the particle list, a dictionary with mass arrays and the respective weights and the\
   analysis-dependent weight format """

  def __init__( self, PartStr=None, zeroweight = None ):
    """ If PartStr != None, the Element is created with particle string \
        If zeroweight != None, fill the weight format with the analysis-dependent weight"""
    self.ParticleStr = ""
    self.MassWeightList = []
    self.WeightFormat = {}
    if PartStr: self.ParticleStr = PartStr
    if zeroweight: self.WeightFormat = zeroweight

  def getParticleList(self):
    """ Converts the particle string in self.ParticleStr to a list of particles"""

    if self.ParticleStr == "": return None
    particles = [[],[]]
    S = self.ParticleStr
    st = S.replace(" ","").replace("'","")
    st = st[st.find("[[["):st.find("]]]")+3]
    branches=[st[2:st.find("]],[[")+1],st[st.find("]],[[")+4:st.find("]]]")+1]]
    for ib,branch in enumerate(branches):
      st = branch
      while "[" in st or "]" in st:
        ptcs = st[st.find("[")+1:st.find("],[")].split(",")
        for ptc in ptcs:
          if not ptc in Reven.values() and not PtcDic.has_key(ptc):
            print "[AElement] unknown particle:",ptc
            return False
        spcts=str(ptcs).replace("'","").replace(" ","")
        st=st.replace(spcts,"",1)
        particles[ib].append(ptcs)

    return particles

  def getEinfo( self ):
    """ Get global topology info from particle string """
    vertnumb = []
    vertparts = []
    ptcs = self.getParticleList()
    for branch in ptcs:
      vertnumb.append(len(branch))
      vertparts.append([len(v) for v in branch])
      if len(vertparts[len(vertparts)-1]) == vertnumb[len(vertnumb)-1]-1:
        vertparts[len(vertparts)-1].append(0)  #Append 0 for stable LSP
    return {"vertnumb" : vertnumb, "vertparts" : vertparts}


class CElement:
  """ A Cluster Element, contains a simple string with the particle list and its weight """

  def __init__( self, PartStr=None, weight = None ):
    """ If PartStr != None, the Element is created with particle string \
        If zeroweight != None, fill the weight format with the analysis-dependent weight"""
    self.ParticleStr = ""
    self.Weight = {}
    if PartStr: self.ParticleStr = PartStr
    if weight: self.Weight = weight

  def getParticleList(self):
    """ Converts the particle string in self.ParticleStr to a list of particles"""

    if self.ParticleStr == "": return None
    particles = [[],[]]
    S = self.ParticleStr
    st = S.replace(" ","").replace("'","")
    st = st[st.find("[[["):st.find("]]]")+3]
    branches=[st[2:st.find("]],[[")+1],st[st.find("]],[[")+4:st.find("]]]")+1]]
    for ib,branch in enumerate(branches):
      st = branch
      while "[" in st or "]" in st:
        ptcs = st[st.find("[")+1:st.find("],[")].split(",")
        for ptc in ptcs:
          if not ptc in Reven.values() and not PtcDic.has_key(ptc):
            print "[AElement] unknown particle:",ptc
            return False
        spcts=str(ptcs).replace("'","").replace(" ","")
        st=st.replace(spcts,"",1)
        particles[ib].append(ptcs)

    return particles

  def getEinfo( self ):
    """ Get global topology info from particle string """
    vertnumb = []
    vertparts = []
    ptcs = self.getParticleList()
    for branch in ptcs:
      vertnumb.append(len(branch))
      vertparts.append([len(v) for v in branch])
      if len(vertparts[len(vertparts)-1]) == vertnumb[len(vertnumb)-1]-1:
        vertparts[len(vertparts)-1].append(0)  #Append 0 for stable LSP
    return {"vertnumb" : vertnumb, "vertparts" : vertparts}


class MassWeight:
  """ A simple object holding a mass array and its respective weight """
  def __init__( self, mass=None, weight = None ):
    self.mass = []
    self.weight = {}
    if mass: self.mass = mass
    if weight: self.weight = weight

class CTop:
  """ cluster global topology. contains a list of cluster elements (CElements),
    the number of vertices and the number of particle insertions in each vertex and the cluster mass"""

  def __init__(self,Top=None,masscluster=None):
    """ if Top != None and masscluster != None creates a cluster topology
        from just the elements in Top.ElList belonging to the masscluster
        (all masses are replaced by the average mass and elements with
        equivalent masses have their weights combined) """

    self.vertnumb = [] ## number of vertices per branch, e.g. [2,2]
    self.vertparts = [] ##  number of insertions per vertex per branch, e.g. [[1,0],[1,0]] for T2
    self.clustermass = None
    self.ElList = [] # list of AElements

    if Top and masscluster:
      if type(Top) != type(ATop()):
        print "[SMSDataObjects.CTop] Input must be an analysis topology (ATop) !"
        return False
      if type(masscluster) != type([]):
        print "[SMSDataObjects.CTop] Input masscluster must be a list !"
        return False
#Compute average mass in cluster
      self.clustermass = ClusterTools.MassAvg(masscluster,"harmonic")
      self.vertnumb = Top.vertnumb
      self.vertparts = Top.vertparts
#Add elements to cluster
      for El in Top.ElList: self.addAnalysisElement(El,masscluster)



  def addAnalysisElement(self,NewElement,masscluster):
    """  Adds an analysis element (AElement) to the corresponding cluster elements (CElements) in ElList \
         Both the particles and branch orderings must be identical!"""
    import copy

    if type(NewElement) != type(AElement()):
      print "[SMSDataObjects.py] wrong input! Must be an AElement object"
      return False

#Consistency checks:
    for OldEl in self.ElList:
      if type(OldElement) != type(CElement()):
        print "[SMSDataObjects.py] wrong input! Elements in ATop must be an AElement object"
        return False
      if OldEl.getEinfo() != NewElement.getEinfo(): return False
    if self.clustermass != ClusterTools.MassAvg(masscluster,"harmonic"):
        print "[SMSDataObjects.py] wrong masscluster input!"
        return False

    newparticles = NewElement.getParticleList()

    for massweight in NewElement.MassWeightList:
      mass = massweight.mass
      weight = massweight.weight
      if not mass in masscluster: continue
  #If mass is in cluster, add element to Cluster Topology
      match = False
      for OldEl in self.ElList:
        oldparticles = OldElement.getParticleList()
        if simParticles(oldparticles,newparticles,useDict=False):
          match = True
          oldweight = OldEl.Weight
          OldEl.Weight = ClusterTools.sumweights([oldweight,weight])
          break

      if not match:
        NewEl = CElement(NewElement.ParticleStr,weight)        
        self.ElList.append(NewEl)
    return True

  def evaluateCluster(self,results):
    """ Evaluates the constraints and conditions in results using the
      respective theoretical cross section predictions for each element in the cluster topology. 

      :type results: a dictionary with the constraints and the conditions \
         of the experiment results

      :returns: an XSecPredictionForCluster object
    """
    from Tools.PhysicsUnits import addunit
    from AuxiliaryFunctions import getelements, eltonum
    import TheoryPrediction
  
  #To store the result:
    ClusterResult = TheoryPrediction.XSecPredictionForCluster()
   
  #Get constraints and conditions:
    consts = results.keys()
    if len(consts) > 1:
      print "evaluateCluster: Analysis contains more than one entry"
      return False    
    conds = results.values()[0]
    if ";" in conds or "Csim" in conds or "Cgtr" in conds:
      conds = conds.rsplit(";")
    else:
      conds = conds.rsplit(",")
      
  #Get a list of all elements appearing in results:
    allEl = set(getelements(consts) + getelements(conds))
    
  #Generate zeroweight and list of numerical elements with zero weights
    zeroweight = {}
    for wk in self.ElList[0].Weight.keys(): zeroweight[wk]=addunit(0.,'fb')
    nEll = [zeroweight]*len(allEl)  
  #Build a dictionary to map the relevant elements to its respective numerical element and the element to its weight
    thdic = {}
    eldic = {}
    iel = 0
    for el in allEl:
      ptcsA = CElement(el).getParticleList()
      nel = "nEl["+str(iel)+"]"
      thdic[el] = nel
      for El in self.ElList:
        ptcsB = El.getParticleList()
        if simParticles(ptcsA,ptcsB,useDict=False): nEll[iel] = El.Weight
      iel += 1

  #Replace string elements by their respective numerical element (nEl):
    consts_num = {}
    conds_num = {}
    for const in consts: consts_num[const] = eltonum(const,thdic)
    for cond in conds: conds_num[cond] = eltonum(cond,thdic)
    
  #Loop over weights and
    const_res = {}
    cond_res = {}
  #Evaluate each constraint
    for ckey in consts_num.keys():
      const = consts_num[ckey]
      res = {}
      for weight in zeroweight.keys():
        nEl = []
        for el in nEll: nEl.append(el[weight])  #Select weight
        res[weight] = Ceval(const,nEl)
      const_res[ckey] = res
  #Evaluate each condition   
      for ckey in conds_num.keys():
        cond = conds_num[ckey]
        res = {}
        for weight in zeroweight.keys():
          nEl = []
          for el in nEll: nEl.append(el[weight])  #Select weight
          res[weight] = Ceval(cond,nEl)
        cond_res[ckey] = res   
      
    ClusterResult.conditions_dic = cond_res
    ClusterResult.result_dic = const_res

    return ClusterResult




class ATop:
  """ analysis global topology. contains a list of analysis elements (AElements),
    the number of vertices and the number of particle insertions in each vertex"""

  def __init__(self):
    self.vertnumb = [] ## number of vertices per branch, e.g. [2,2]
    self.vertparts = [] ##  number of insertions per vertex per branch, e.g. [[1,0],[1,0]] for T2
    self.ElList = [] # list of AElements

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

  def __eq__ ( self, other ):
    return self.isEqual ( other )

  def leadingElement ( self ):
    """ often, a topology carries only one element, so
      we have a special accessor for this """
    if len(self.ElList)==0: return None
    return self.ElList[0]

  def addEventElement(self,NewElement):
    """  Adds an event element (EElement) to the corresponding analysis elements (AElements) in ElList\
         The event element and analysis elements DO NOT need to have the same branch ordering"""

    import copy

    if type(NewElement) != type(EElement()):
      print "[SMSDataObjects.py] wrong input! Must be an EElement object"
      return False

    newparticles_a = [NewElement.B[0].particles,NewElement.B[1].particles]
    newparticles_b = [NewElement.B[1].particles,NewElement.B[0].particles]

    for OldElement in self.ElList:
      if type(OldElement) != type(AElement()):
        print "[SMSDataObjects.py] wrong input! Elements in ATop must be an AElement object"
        return False

      oldparticles = OldElement.getParticleList()
#Check if particles match
      if not simParticles(newparticles_a,oldparticles) and not simParticles(newparticles_b,oldparticles): return False
#Format the new weight to the analysis-dependent format (remove weights which do not match the analysis format and add zero to missing weights)
      neweight = copy.deepcopy(NewElement.weight)
      for key in OldElement.WeightFormat.keys()+neweight.keys():
        if not neweight.has_key(key): neweight[key] = addunit(0.,'fb')
        if not OldElement.WeightFormat.has_key(key): neweight.pop(key)
    
#Check if masses match
      added = False
      OldEl = EElement(OldElement.ParticleStr)  #Create temporary EElement for easy comparison
      for massweight in OldElement.MassWeightList:
        OldEl.B[0].masses = massweight.mass[0]
        OldEl.B[1].masses = massweight.mass[1]
        OldEl.weight = massweight.weight
        if NewElement.isSimilar(OldEl,order=False):
          massweight.weight = ClusterTools.sumweights([OldEl.weight,neweight])
          added = True
          break   #To avoid double counting only add the event weight to one mass combination

#If no identical mass was found, add entry to mass-weight dictionary
      if not added:
        newmass = [NewElement.B[0].masses,NewElement.B[1].masses]
        if not NewElement.isSimilar(OldEl,order=True,igmass=True): newmass = [newmass[1],newmass[0]] #Check for correct branch ordering
        newmassweight = MassWeight(newmass,neweight)
        OldElement.MassWeightList.append(newmassweight)


    return True




class GTop:
  """ global topology. contains a list of event elements (EElements),
    the number of vertices and the number of particle insertions in each vertex"""

  def __init__(self):
    self.vertnumb = [] ## number of vertices per branch, e.g. [2,2]
    self.vertparts = [] ##  number of insertions per vertex per branch, e.g. [[1,0],[1,0]] for T2
    self.ElList = [] # list of EElements

  def leadingElement ( self ):
    """ often, a topology carries only one element, so
      we have a special accessor for this """
    if len(self.ElList)==0: return None
    return self.ElList[0]

  def elements ( self ):
    return self.ElList

  def describe(self):
    """ a lengthy description """
    ret="number of vertices=%s number of vertex particles=%s number of elements=%d" % \
        ( self.vertnumb, self.vertparts, len(self.ElList) )
    return ret

  def __str__(self):
    return str(self.vertparts).replace(" ","")

  def checkConsistency ( self, verbose=False ):
    """ the number of vertices and insertions per vertex is
      redundant information in a topology, so we can perform
      an internal consistency check """
    for element in self.ElList:
      info=element.getEinfo()
      if self.vertnumb!=info["vertnumb"]:
        if verbose: print "[SMSDataObjects.py] inconsistent topology!!!"
        return False
      if self.vertparts!=info["vertparts"]:
        if verbose: print "[SMSDataObjects.py] inconsistent topology!!!"
        return False
    if verbose: print "[SMSDataObjects.py] topology is consistent."
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

  def __eq__ ( self, other ):
    return self.isEqual ( other )

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


  def massCompressedTopology ( self, mingap ):
    """ if two masses in this topology are degenerate, create
        a compressed copy of this topology """
    import copy
    from Tools.PhysicsUnits import rmvunit
    mingap=rmvunit ( mingap, "GeV" )
    ETopComp = copy.deepcopy(self)
  #Loop over branches
    for ib in range(len(ETopComp.vertnumb)):
      if ETopComp.vertnumb[ib] < 2: continue
  #Remove all external particles between compressed masses
      for ivert in range(ETopComp.vertnumb[ib]-1):
        massA = rmvunit(ETopComp.ElList[0].B[ib].masses[ivert],"GeV")
        massB = rmvunit(ETopComp.ElList[0].B[ib].masses[ivert+1],"GeV")
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


class TopologyList:
  """ implements a list of topologies, knows how to correctly add 
      a topology """
  def __init__ ( self, topos=[] ):
    """ if topos are given, we add all of them sequentially """
    self.topos=[]
    for topo in topos:
      self.add ( topo )

  def __len__ ( self ): 
    return len(self.topos)

  def __getitem__ ( self, n ):
    return self.topos[n]

  def addList ( self, List ):
    for topo in List: self.add ( topo )

  def __str__ ( self ):
    s="TopologyList:\n" 
    for topo in self.topos:
      s+=str(topo)+"\n"
    return s

  def describe( self ):
    s="TopologyList:\n" 
    for topo in self.topos:
      s+=str(topo)+"\n"
    return s

  def add ( self, topo ):
    """ Check if elements in topo matches an entry in self.topos. If it does,
    add weight.  If the same topology exists, but not the same element, add
    element.  If neither element nor topology exist, add the new topology and
    all its elements 

    :type topo: GTop
    
    """
    import copy
    for (inew,element) in enumerate(topo.ElList): ## range(len(topo.ElList)):
      if len(element.B)<2:
        print "[SMSDataObjects.TopologyList] error: assumed at least two branches"
        continue
      NewEl_a = element
      NewEl_b = copy.deepcopy(NewEl_a)
      NewEl_b.B[1] = element.B[0]
      NewEl_b.B[0] = element.B[1]   #Check both orderings
      equaltops = -1
      equalels = -1
      i = -1
      while (equaltops < 0 or equalels < 0) and i < len(self.topos)-1:
        i += 1
        if topo.isEqual(self.topos[i],order=False):  #First look for matching topology
          equaltops = i
        else: continue

        for j in range(len(self.topos[i].ElList)):  #Search for matching element
          OldEl = self.topos[i].ElList[j]
          if OldEl.isEqual(NewEl_a):
            equalels = j
            NewEl = NewEl_a
            break
          elif OldEl.isEqual(NewEl_b):
            equalels = j
            NewEl = NewEl_b
            break


  #If element exists, add weight:
      if equalels >= 0:
        if len(OldEl.weight) != len(NewEl.weight):
          print "Wrong number of weights"
        else:
          w1 = OldEl.weight
          w2 = NewEl.weight
          self.topos[equaltops].ElList[equalels].weight = ClusterTools.sumweights([w1,w2])



  #When combining elements, keep the smallest set of PDG mother IDs (not used in the analysis, only relevant to set a standard):
          if min(abs(NewEl.B[0].momID),abs(NewEl.B[1].momID)) < min(abs(OldEl.B[0].momID),abs(OldEl.B[1].momID)):
            for ib in range(2): self.topos[equaltops].ElList[equalels].B[ib].momID = NewEl.B[ib].momID


  #If topology and/or element does not exist, add:
    if equaltops == -1:
      self.topos.append(topo)
    elif equalels == -1:
      if topo.isEqual(self.topos[equaltops],order=True):
        NewEl = NewEl_a
      else:
        NewEl = NewEl_b
      if not self.topos[equaltops].addElement(NewEl):
        print "Error adding element"
        print '\n'

