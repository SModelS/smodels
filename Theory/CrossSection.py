#!/usr/bin/env python

"""
.. module:: CrossSection
    :synopsis: A class that encapsulates the result of the computation of the reference cross section.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

class CrossSection:
  """ A wrapper around this complicated result dictionary, to 
      make it easier to use this dictionary """

  def __init__ ( self, data ):
    self.data=data

  # make it behave much like a dictionary
  def __len__ ( self ): return len(self.data)
  def __getitem__ ( self, i ): return self.data[i]
  def items ( self ): return self.data.items()
  def __str__ ( self ): return str(self.data)
  def keys ( self ): return self.data.keys()  

  def weights ( self ): return self.data["Wdic"]
  def crossSections ( self ): return self.data["Xsecdic"]
  def crossSectionsInfo ( self ): return self.data["XsecList"]

  def getCrossSection ( self, pidmom1, pidmom2, order="NLL", sqrts=8 ):
    """ returns production cross section for a given production mode (pidmom1, pidmom2)
        at 7 or 8 TeV in chosen order LO or NLL """
        
    from Tools.PhysicsUnits import addunit        
    if pidmom1 > pidmom2:
      pidmom1, pidmom2 = pidmom2, pidmom1

    if type(order) == type('str'):
      if order.lower() == "nll": order = 2
      elif order.lower() == "nlo": order = 1
      elif order.lower() == "lo": order = 0
    if type(order) != type(1):
      print '[CrossSection]: Unknown cross section order:',order
      return None

    if type(sqrts) == type(1) or type(sqrts) == type(1.):
      sqrts = addunit(float(sqrts),'TeV')
    if type(sqrts) != type(addunit(1.,'TeV')):
      print '[CrossSection]: Unknown sqrts value:',sqrts
      return None

#Get label: 
    for xsec in self.crossSectionsInfo.xsecs:
      if xsec.order != order: continue
      if xsec.sqrts != sqrts: continue
      k = xsec.label
      break
      
    if not k or not self.crossSections().has_key(k):
        print '[CrossSection]: No cross sections found for',k
        return None
    allxsecs=self.crossSections()[k]
    if not allxsecs:
      print '[CrossSection]: No cross sections for ',k
      return None
    if allxsecs.has_key((pidmom1,pidmom2)):
      return allxsecs[(pidmom1,pidmom2)]
    else:
      print '[CrossSection]: Cross Sections only available for %s' %str(allxsecs.keys())
      return None
      
    return None

  def getSumOfCrossSections ( self, pidmoms, order="NLL", sqrts=8 ):
    """takes a list of pids and returns the integrated production cross section
       for all combinations of pids from the given list at 7 or 8TeV 
       with order = LO or NLL"""
    from Tools.PhysicsUnits import addunit, rmvunit
    
    if type(order) == type('str'):
      if order.lower() == "nll": order = 2
      elif order.lower() == "nlo": order = 1
      elif order.lower() == "lo": order = 0
    if type(order) != type(1):
      print '[CrossSection]: Unknown cross section order:',order
      return None

    if type(sqrts) == type(1) or type(sqrts) == type(1.):
      sqrts = addunit(float(sqrts),'TeV')
    if type(sqrts) != type(addunit(1.,'TeV')):
      print '[CrossSection]: Unknown sqrts value:',sqrts
      return None
  
  #Get label: 
    for xsec in self.crossSectionsInfo.xsecs:
      if xsec.order != order: continue
      if xsec.sqrts != sqrts: continue
      k = xsec.label
      break
  
    
    if not k or not self.crossSections().has_key(k):
      print '[CrossSection]: No cross sections for ',k
      return None
    allxsecs=self.crossSections()[k]
    if not allxsecs:
      print '[CrossSection]: No cross sections for',k
      return None
    Sum=0
    for (key,value) in allxsecs.items():
      if abs(key[0]) in pidmoms and abs(key[1]) in pidmoms:
        value=rmvunit(value, 'fb')
#        print 'k:', key, 'v:',value
        Sum+=value
    return Sum

  def crossSectionLightSquarks ( self, order="NLL", sqrts=8 ):
    """returns the integrated production cross section of all light squarks productions
       at sqrts and order = order"""
    squarks=[1000001,1000002,1000003,1000004,2000001,2000002,2000003,2000004]
    return self.getSumOfCrossSections ( squarks, order, sqrts )

  def lhefile ( self, sqrts=8, order=0 ):
    from Tools.PhysicsUnits import addunit
    
    if type(sqrts) == type(1) or type(sqrts) == type(1.):
      sqrts = addunit(float(sqrts),'TeV')
    if type(sqrts) != type(addunit(1.,'TeV')):
      print '[CrossSection]: Unknown sqrts value:',sqrts
      return None
  
    #Get label:
    for xsec in self.crossSectionsInfo().xsecs:
      if xsec.order != order: continue
      if xsec.sqrts != sqrts: continue
      k = xsec.label
      break
    
    if not k or not self.data["lhefiles"].has_key(k):
      print "[CrossSection.py] lhefile for",sqrts,"does not exist."
      return None
    else:    
      return self.data["lhefiles"][k]

class SingleXSecInfo:
  """A simple class to store the information about a single cross-section (center of mass, order and label)
     order = 0 (LO), 1 (NLO) or 2 (NLL). If Str != None, creates from string (format = 7 tev NLO, 8 TeV (NLL),..)"""
  def __init__ (self, Str=None ):
    from Tools.PhysicsUnits import addunit
    self.sqrts=None
    self.order=None
    self.label=None
    if Str:
      orders = ['lo','nlo','nll']
      inputerror = False
      inStr = Str.replace(' ','').lower()
      for iorder,order in enumerate(orders):
        if order in inStr: self.order = iorder
      if self.order is None or not 'tev' in inStr: inputerror = True
      sqrt=inStr[:inStr.index('tev')]
      iS = 1
      sqrtS = 0.
      while iS and iS <= len(sqrt):
        try:
          sqrtS = eval(sqrt[-iS:])
          iS += 1
        except:
          iS = 0
      if type(sqrtS) != type(1.) and type(sqrtS) != type(1): inputerror = True
      if inputerror:
        print "[XSecComputer.py] Unknown input cross-section format, ignoring",inStr
        return
      else:
        self.sqrts = addunit(abs(float(sqrtS)),'TeV')
        self.label = Str.lstrip().rstrip()

  def __eq__(self,other):
    if type(other) != type(SingleXSecInfo()): return False
    if other.sqrts != self.sqrts: return False
    if other.order != self.order: return False
    return True

  def __str__ ( self ):
    """cross-section information in string format"""
    from Tools.PhysicsUnits import rmvunit
    st = 'label: '+self.label+', sqrts: '+str(rmvunit(self.sqrts,'TeV'))+' TeV, order: '+str(self.order)
    return st

      
class XSecInfoList:
  """A simple class to store the list of cross-sections to be used"""    
  
  def __init__ (self, Str='7 TeV (LO), 8 TeV (LO), 7 TeV (NLL), 8 TeV (NLL)'):
    """Creates a list of SingleXSecInfo objects from the input string. All cross-sections must be in TeV!"""
    self.xsecs=[]
    if Str:
      inStrs = Str.split(',')
    else:
      inStrs = []
    for inStr in inStrs:
      try:
        Xsec = SingleXSecInfo(inStr)
        if type(Xsec) == type(SingleXSecInfo()): self.xsecs.append(Xsec)
      except:
        pass        
          
  def __getitem__ ( self, label=None ):
    """Return the SingleXSecInfo object with the respective label"""
    if not label: return None
    for xsec in self.xsecs:
      if xsec.label == label: return xsec
    return None
          
  def sort(self,reverse=False):
    """Sort cross-sections by center of mass energy"""
    self.xsecs = sorted(self.xsecs, key=lambda xsec: xsec.sqrts, reverse=reverse)
    return True

  def getSqrts(self):
    """Return a list with all sqrts values appearinf in xsecs"""
    allsqrts = []
    for xsec in self.xsecs:
      if not xsec.sqrts in allsqrts: allsqrts.append(xsec.sqrts)
    return allsqrts
