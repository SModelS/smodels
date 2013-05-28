from ParticleNames import Rodd, Reven, PtcDic, ptype, simParticles
import copy
from ClusterTools import sumweights

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
