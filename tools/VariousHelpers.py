#!/usr/bin/env python

"""
.. module:: VariousHelpers
    :synopsis: Various helper classes.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

def getMaxLum(List):
  """ Goes through all analyses in the list and returns the maximum luminosity.
  If one of the luminosities has not been define, return None. """
  from PhysicsUnits import addunit
  maxlum = addunit(0.,'fb-1')
  for Ana in List:
    if type(Ana.lum) == type(1) or type(Ana.lum) == type(1.): Ana.lum = addunit(Ana.lum,'fb-1')
    if Ana.lum.asNumber():
      maxlum = max(maxlum,Ana.lum)
    else:
      return None
  return maxlum

def nCPUs():
  """ obtain the number of CPU cores on the machine, for several
      platforms and python version """
  try:
    import multiprocessing
    return multiprocessing.cpu_count()
  except Exception,e:
    pass
  try:
    import psutil
    return psutil.NUM_CPUS
  except Exception,e:
    pass
  try:
    import os
    res = int(os.sysconf('SC_NPROCESSORS_ONLN'))
    if res>0: return res
  except Exception,e:
    pass
  return None

##import logging 
##FORMAT = '%(levelname)s in %(module)s.%(funcName)s(): %(message)s'
##import sys
#### workaround for python2.4
##if sys.version_info[0]<3 and sys.version_info[1]<5:
##  FORMAT = '%(levelname)s in %(module)s.s(): %(message)s'
### http://docs.python.org/2/library/logging.html#logging.Formatter
##logging.basicConfig(level=logging.INFO,format=FORMAT)
