#!/usr/bin/python

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

  def weights ( self ): return self.data["Wdic"]
  def crossSections ( self ): return self.data["Xsecdic"]

  def getCrossSection ( self, pidmom1, pidmom2, order="NLL", sqrts=8 ):
    """ returns production cross section for a given production mode (pidmom1, pidmom2)
        at 7 or 8 TeV in chosen order LO or NLL """
    if pidmom1 > pidmom2:
      pidmom1, pidmom2 = pidmom2, pidmom1
    k='%d TeV (%s)' %(sqrts, order)
    if not self.crossSections().has_key(k):
      print '[CrossSection]: No cross sections for %s at %d TeV.' %(order,sqrts)
      return None
    allxsecs=self.crossSections()[k]
    if not allxsecs:
      print '[CrossSection]: No cross sections for %s at %d TeV.' %(order,sqrts)
      return None
    if allxsecs.has_key((pidmom1,pidmom2)):
      return allxsecs[(pidmom1,pidmom2)]
    else:
      print '[CrossSection]: Cross Sections only available for %s' %str(allxsecs.keys())
      return None

  def getSumOfCrossSections ( self, pidmoms, order="NLL", sqrts=8 ):
    """ takes a list of pids and returns the integrated production cross section
        for all combinations of pids from the given list at 7 or 8TeV 
        with order = LO or NLL """
    from Tools.PhysicsUnits import rmvunit
    k='%d TeV (%s)' %(sqrts, order)
    if not self.crossSections().has_key(k):
      print '[CrossSection]: No cross sections for %s at %d TeV.' %(order,sqrts)
      return None
    allxsecs=self.crossSections()[k]
    if not allxsecs:
      print '[CrossSection]: No cross sections for %s at %d TeV.' %(order,sqrts)
      return None
    Sum=0
    for (key,value) in allxsecs.items():
      if abs(key[0]) in pidmoms and abs(key[1]) in pidmoms:
        value=rmvunit(value, 'fb')
#        print 'k:', key, 'v:',value
        Sum+=value
    return Sum

  def crossSectionLightSquarks ( self, order="NLL", sqrts=8 ):
    """ returns the integrated production cross section of all light squarks productions
        at 7 or 8TeV with order = LO or NLL """
    squarks=[1000001,1000002,1000003,1000004,2000001,2000002,2000003,2000004]
    return self.getSumOfCrossSections ( squarks, order, sqrts )

  def lhefile ( self, sqrts ):
    from Tools.PhysicsUnits import rmvunit
    sqrts=rmvunit(sqrts,"TeV")
    if sqrts==7: return self.data["lhe7file"]
    if sqrts==8: return self.data["lhe8file"]
    print "[CrossSection.py] lhefile for",sqrts,"does not exist."
    return None
