#!/usr/bin/env python

"""
.. module:: SMSDataObjects
    :synopsis: our basic data objects: TopologyList, GTop, EElement, BElement
    
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>, Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

import ClusterTools
from ParticleNames import Reven, PtcDic, simParticles
from AuxiliaryFunctions import Ceval

class BElement:
  """ A branch-element """

  def __init__(self, S=None ):
    """ A branch-element can be constructed from a string S """
    self.masses = []
    self.particles = []
    self.momID = 0
    if type(S)==type(""):
      st = S.replace(" ","")
      # print "[SMSDataObjects.py] Construct a BElement from a string: st=",st
      while "[" in st or "]" in st:
        ptcs = st[st.find("[")+1:st.find("],[")].split(",")
        for ptc in ptcs:
          if not ptc in Reven.values() and not PtcDic.has_key(ptc):
            print "[BElement] unknown particle:",ptc
            return
        spcts=str(ptcs).replace("'","").replace(" ","")
        st=st.replace(spcts,"",1)
        self.particles.append ( ptcs )
      # print "[SMSDataObjects] ptcs=",ptcs

  def isEqual ( ElA, elB, order=True ):
    if simParticles(ElA.particles,ElB.particles,useDict=False): return False
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
    ret="particles=%s masses=%s" % \
       ( self.particles, [ rmvunit(x,"GeV") for x in self.masses ] )
    return ret

class EElement:
  """ An Event Element, contains of several branches and weight information """
  def __init__(self, S=None ):
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
  def getEinfo(self):
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

  def isEqual ( ElA, ElB,order=True,igmass=False ):
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


class GTop:
  """ global topology. contains a list of elements,
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
    ret="number of vertices=%s number of vertex particles=%s" % \
        ( self.vertnumb, self.vertparts )
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

  def evaluateCluster(self, instr):
    """ Given a mathematical expression in string format involving elements, it evaluates the expression
    using the respective theoretical predictions for each element."""
    from Tools.PhysicsUnits import addunit

    outstr = instr.replace(" ","")
  #Get ordered list of elements:
    El = []
    iels = 0
    while "[[[" in outstr:  #String has element
      st = outstr[outstr.find("[[["):outstr.find("]]]")+3] #Get duplet
      element=EElement ( st )
      # ptclist = strtoel(st)   # Get particle list
      ptclist = element.allParticles()
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
              NewTop.ElList[iel2].weight = ClusterTools.sumweights([oldweight,weight])
              break
          if not match: NewTop.addElement(Elm)
    return NewTop

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

