import sys
sys.path.append ( "../" )

from Tools.PhysicsUnits import addunit, rmvunit
import copy
from ParticleNames import Rodd, Reven, PtcDic, ptype, simParticles
from ClusterTools import DoCluster, GoodMass, MassAvg, sumweights, MassDist,\
       ClusterDist

class BElement:
  """ A branch-element """

  def __init__(self, S=None ):
    """ A branch-element can be constructed from a string S """
    self.masses = []
    self.particles = []
    self.momID = 0
    if type(S)==type(""):
      st = S.replace(" ","")
      # print "[SMSmethods.py] Construct a BElement from a string: st=",st
      while "[" in st or "]" in st:
        ptcs = st[st.find("[")+1:st.find("],[")].split(",")
        for ptc in ptcs:
          if not ptc in Reven.values() and not PtcDic.has_key(ptc):
            print "[BElement] unknown particle:",ptc
            return
        spcts=str(ptcs).replace("'","").replace(" ","")
        st=st.replace(spcts,"",1)
        self.particles.append ( ptcs )
      # print "[SMSmethods] ptcs=",ptcs

  def isEqual ( ElA, elB, order=True ):
    if simParticles(ElA.particles,ElB.particles,useDict=False): return False
    if ElA.masses != ElB.masses: return False
    return True

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
    """ Compare two BElements or Eelements.
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


def AddToList(SMSTop,SMSTopList):
  """ Check if elements in SMSTop matches an entry in SMSTopList. If it does,
  add weight.  If the same topology exists, but not the same element, add
  element.  If neither element nor topology exist, add the new topology and
  all its elements """

  for inew in range(len(SMSTop.ElList)):
    NewEl_a = SMSTop.ElList[inew]
    NewEl_b = copy.deepcopy(NewEl_a)
    NewEl_b.B[1] = SMSTop.ElList[inew].B[0]
    NewEl_b.B[0] = SMSTop.ElList[inew].B[1]   #Check both orderings
    equaltops = -1
    equalels = -1
    i = -1
    while (equaltops < 0 or equalels < 0) and i < len(SMSTopList)-1:
      i += 1
      if SMSTop.isEqual(SMSTopList[i],order=False):  #First look for matching topology
        equaltops = i
      else: continue

      for j in range(len(SMSTopList[i].ElList)):  #Search for matching element
        OldEl = SMSTopList[i].ElList[j]
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
        SMSTopList[equaltops].ElList[equalels].weight = sumweights([w1,w2])



#When combining elements, keep the smallest set of PDG mother IDs (not used in the analysis, only relevant to set a standard):
        if min(abs(NewEl.B[0].momID),abs(NewEl.B[1].momID)) < min(abs(OldEl.B[0].momID),abs(OldEl.B[1].momID)):
          for ib in range(2): SMSTopList[equaltops].ElList[equalels].B[ib].momID = NewEl.B[ib].momID


#If topology and/or element does not exist, add:
  if equaltops == -1:
    SMSTopList.append(SMSTop)
  elif equalels == -1:
    if SMSTop.isEqual(SMSTopList[equaltops],order=True):
      NewEl = NewEl_a
    else:
      NewEl = NewEl_b
    if not SMSTopList[equaltops].addElement(NewEl):
      print "Error adding element"
      print '\n'


def strtoel(invar):
  """ Converts an SMS description string to a nested particle list """
  st = invar.replace(" ","")
  st = st[st.find("[[["):st.find("]]]")+3]
  st_B = []
  st_B.append(st[2:st.find("]],[[")+1])
  st_B.append(st[st.find("]],[[")+4:st.find("]]]")+1])

  ptclist = [[],[]]
  for ib in range(2):
    while "[" in st_B[ib] or "]" in st_B[ib]:
      ptcs = st_B[ib][st_B[ib].find("[")+1:st_B[ib].find("],[")].split(",")

#Syntax check:
      for ptc in ptcs:
        if not ptc in Reven.values() and not PtcDic.has_key(ptc):
          print "[strtoel]: Unknown particle:",ptc
          return False

      ptclist[ib].append(ptcs)
      sptcs = str(ptcs).replace("'","")
      sptcs = str(sptcs).replace(" ","")
      st_B[ib] = st_B[ib].replace(sptcs,"",1)

  #print "ptclist=",ptclist
  return ptclist

