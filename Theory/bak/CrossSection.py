#!/usr/bin/env python

"""
.. module:: CrossSection
    :synopsis: A class that encapsulates the result of the computation of the reference cross section.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from Tools.PhysicsUnits import addunit, rmvunit
import copy
import logging

UseXSecs = None

class XSectionInfo:
 """A simple class to store information about the cross-section (center of mass, order and label)"""
 def __init__ (self):
   self.sqrts = None
   self.order = None
   self.label = None

 def __eq__(self,other):
   if type(other) != type(XSectionInfo()): return False
   if other.sqrts != self.sqrts: return False
   if other.order != self.order: return False


class XSection:
  """A simple class to store the information about a single cross-section (value, paritcle ids, center of mass, order and label)
     order = 0 (LO), 1 (NLO) or 2 (NLL)."""
  def __init__ (self):
    self.info = XSectionInfo()
    self.value = None
    self.pid = None

  def __mul__(self,other):
    newXsec = copy.deepcopy(self)
    if type(other) == type(1.):
      newXsec.value = newXsec.value*other
    else:
      print "[XSection.mul]: Xsections can only be multiplied by floats"
      return False
    return newXsec
    
  def __add__(self,other):
    if type(other) == type(XSection()):
      if self.info == other.info:
        res = copy.deepcopy(self)
        res.value += other.value
        return res
    print "[XSection.add]: Trying to add",type(other),"to a XSection objetc"
    return False

  def __eq__(self,other):
    if type(other) != type(XSection()): return False
    if other.info != self.info: return False
    if other.value != self.value: return False
    if other.pid != self.pid: return False
    return True

  def __str__ (self):
    """cross-section information in string format"""
    st = 'label: '+self.label+', sqrts: '+str(rmvunit(self.sqrts,'TeV'))+' TeV, order: '+str(self.order)+', value:'+str(self.value)
    return st

     
class XSectionList:
  """A simple class to store the list of cross-sections to be used"""    
  
  def __init__ (self):
    """Creates a list of XSection objects from the input string with None cross-section values."""
    self.XSections=[]

  def __mul__(self,other):
    newList = copy.deepcopy(self)
    for ixsec,xsec in enumerate(newList.XSections): newList.XSections[ixsec] = xsec*other
    return newList

  def __str__(self):
    return str(self.getDictionary())


  def getXsecsFor(self,item):
    """ Returns a list of XSection objects for item (label, pid, sqrts) """
    xsecList = XSectionList()
    for xsec in self.XSections:
      if type(item) == type(xsec.info.label) and item == xsec.info.label:
        xsecList.XSections.append(xsec)
      elif type(item) == type(xsec.info.sqrts) and item == xsec.info.sqrts:  
        xsecList.XSections.append(xsec)
      elif type(item) == type(xsec.pid) and item == xsec.pid:  
        xsecList.XSections.append(xsec)
      elif type(item) == type(1) and (item in xsec.pid):
        xsecList.XSections.append(xsec)

    return xsecList

  def delete(self,XSec):
    """Deletes the cross-section entry from the list"""
    for ixsec,xsec in enumerate(self.XSections):
      if xsec == XSec: self.XSections.pop(ixsec)

  def getInfo(self):
    """Gets the basic info about the cross-sections appearing in the list (order,value and label). Returns a list of XSectionInfo objects """
    allInfo = []
    for xsec in self.XSections:
      info = xsec.info
      if not info in allInfo: allInfo.append(info)
    return allInfo

  def getLabels(self):
    """Gets all the labels appearing in the list."""
    allLabels = []
    allInfo = self.getInfo()
    for info in allInfo: allLabels.append(info.label)
    return list(set(allLabels)) 

  def getPIDpairs(self):
    """Gets all the particle ID pairs appearing in the list."""
    allPidPairs = []
    for xsec in self.XSections: allPidPairs.append(xsec.pid)
    return list(set(allPidPairs)) 

  def getPIDs(self):
    """Gets all the particle IDs appearing in the list."""
    allPids = []
    for xsec in self.XSections: allPids.extend(xsec.pid)
    return list(set(allPids)) 


  def getMaxXsec(self):
    """Gets the maximum cross-section value appearing in the list."""
    maxxsec= addunit(0.,'fb')
    for xsec in self.XSections:
      if xsec.value > maxxsec: maxxsec = xsec.value
    return maxxsec

  def makeUniformLabels(self):
    """ Fills the cross-section list with zero cross-section entries, so all particles have cross-section entries for all labels """

    global UseXSecs

    allPids = self.getPIDpairs()
    allInfo = self.getInfo()
    if not UseXSecs is None: allInfo = UseXSecs
    for pid in allPids:
      Xsecs = self.getXsecsFor(pid)
      hasLabels = [info.label for info in Xsecs.getInfo()]
      for info in allInfo:
        if not info.label in hasLabels:
           newentry = XSection()
           newentry.value = addunit(0.,'fb')
           newentry.pid = pid
           newentry.info = info
           self.XSections.append(newentry)


  def getDictionary(self,groupBy="pids"):
    """ Converts the list of XSection objects to a nested dictionary. First level keys are the particles IDs (if groupBy=pids) or labels
    (if groupBy=labels) and values are the cross-section labels or particle IDs and the cross-section value."""

    XsecDic = {}
    self.makeUniformLabels()  #Make sure all particles have entries for all labels

    if groupBy == "pids":
      allPids = self.getPIDpairs()
      for pid in allPids:
        XsecDic[pid] = {}
        Xsecs = self.getXsecsFor(pid)
        for xsec in Xsecs.XSections: XsecDic[pid][xsec.info.label] = xsec.value        

    elif groupBy == "labels":
      allLabels = self.getLabels()
      for label in allLabels:
        XsecDic[label] = {}
        Xsecs = self.getXsecsFor(label)
        for xsec in Xsecs.XSections:
          XsecDic[label][xsec.pid] = xsec.value

    return XsecDic
    

def getXsecFromFile(slhafile):
  """ obtain dictionary cross-sections from input SLHA file. """

  global UseXSecs

  XsecsInFile = XSectionList()  #To store information about all cross-sections in the SLHA file
  slha = open(slhafile, 'r')
  lines = slha.readlines()
  xsecblock = False
  for l in lines:
    if l.startswith("#") or len(l)<2: continue
    if 'XSECTION' in l:
      xsecblock = True
      sqrtS =  eval(l.split()[1])/1000.    #Values in the SLHA file are in GeV
      pids = (eval(l.split()[5]),eval(l.split()[6]))
      continue
    if not xsecblock: continue  #ignores other entries
    cs_order = eval(l.split()[1])
    cs = addunit(eval(l.split()[6]),'fb')
    wlabel = str(int(sqrtS))+' TeV'
    if cs_order == 0: wlabel += ' (LO)'
    elif cs_order == 1: wlabel += ' (NLO)'
    elif cs_order == 2: wlabel += ' (NLL)'
    else:
      print '[SLHADecomposer] unknown QCD order in XSECTION line', l
      return False
    xsec = XSection()
    xsec.info.sqrts = addunit(sqrtS,'TeV')
    xsec.info.order = cs_order
    xsec.info.label = wlabel
    xsec.value = cs
    xsec.pid = pids
    XsecsInFile.XSections.append(xsec)

  slha.close()

  if UseXSecs is None:
    UseXSecs = XsecsInFile.getInfo()
    log = logging.getLogger(__name__)
    log.warning ( "Cross-section information not found. Using values from SLHA file" )
  else:
    for xsec in XsecsInFile.XSections:
      if not xsec.info in UseXSecs: XsecsInFile.XSections.delete(xsec)   #Remove entries which do not match the previously defined cross-sections

  return XsecsInFile
