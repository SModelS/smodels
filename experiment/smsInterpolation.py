#!/usr/bin/env python

"""
.. module:: smsInterpolation
    :synopsis: smsInterpolation is called by smsResults.getSmartUpperLimit.
        UpperLimit takes arbitrary input masses and checks if there is a
        corresponding upper limit for the given analysis and topology.
        The upper limit is returned in 'pb'.
        If several histograms with different x-values are available,
        an interpolation is performed.

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""

import smsResults, smsHelpers
import numpy as np
from scipy.interpolate import griddata
from tools.physicsUnits import rmvunit, addunit

def getxval(mx, my, mz,mass=False):
  """Calculate x-value for one point.
      mx=Mother-mass, my=LSP-mass, mz is information on the intermediate mass as given in the axes-information.
      If mass=True is selected: return intermediate mass instead of x-value."""
  mx=rmvunit(mx,"GeV")
  my=rmvunit(my,"GeV")
  mz=rmvunit(mz,"GeV")
  if mz.find('x')==-1 and mz.find('C')==-1 and mz.find('y')==-1:
    xfac = float(mz)/100
    if mass: return xfac*mx+(1-xfac)*my
    return float(mz)/100
  if mz.find('x')>-1: tx = float(mz[mz.find('x')+1:mz.find('x')+4])
  else: tx=None
#  if mz.find('y')>-1: ty = float(mz[mz.find('y')+1])
#  else: ty=None
  if mz.find('C')>-1: c = float(mz.split("C")[1])
  else: c=None
  z = 0.
  if tx: z += tx*my/100
#  if ty: z += ty*mx
  if c: z+=c
  if mass: return z
  # print "[smsInterpolation.py] z=",z,"mx=",mx,"my=",my
  xval = (z-my)/(mx-my)
  return xval


def doGridData(ana,topo,masses,d,debug=True,run=None):
  """ make np.array and use scipy.griddata function for ana, topo """
  masslist=[]
  ullist=[]

  for ds in d:
    if not ds['mz']:
      if debug: print "[smsInterpolation] error: No information on intermediate mass availabel for %s/%s." % ( ana, topo )
      return None

    D=None
    L=None
    M1=None

    if ds['mz'][0].find('D')>-1: D=float(ds['mz'][0].split('=')[1])
    elif ds['mz'][0].find('LSP')>-1: L=float(ds['mz'][0].split("LSP")[1])
    elif ds['mz'][0].find('M1')>-1: M1=float(ds['mz'][0].split("M1")[1])

    ul_dict=smsHelpers.getUpperLimitDictionary(ana,getHistName(topo,ds['mz'][0]),run)

    for x in ul_dict:
      for y in ul_dict[x]:

        if D:
          massv=[0.,0.,0.]
          massv[getaxis('x',ds['axes'])]=x
          massv[getaxis('y',ds['axes'])]=y
          if massv[getaxis('x',ds['mz'][0])]==0.: massv[getaxis('x',ds['mz'][0])]=massv[getaxis('y',ds['mz'][0])]+D
          if massv[getaxis('y',ds['mz'][0])]==0.: massv[getaxis('y',ds['mz'][0])]=massv[getaxis('x',ds['mz'][0])]-D

        elif L:
          massv=[0.,0.,L]
          massv[getaxis('x',ds['axes'])]=x
          massv[getaxis('y',ds['axes'])]=y

        elif M1:
          massv=[M1,0., 0.]
          massv[getaxis('x',ds['axes'])]=x
          massv[getaxis('y',ds['axes'])]=y


        else: massv=[x,getxval(x,y,ds['mz'][0],mass=True),y]

        masslist.append(massv)
        ullist.append(ul_dict[x][y])

  p=np.array(masslist)
  v=np.array(ullist)

  mx=rmvunit(masses[0],"GeV")
  my=rmvunit(masses[1],"GeV")
  mz=rmvunit(masses[2],"GeV")

  r=griddata(p,v,(mx,my,mz),method="linear")

  if np.isnan(r):
    if debug: print "[smsInterpolation] error: masses out of range for %s/%s (no extrapolation)" % ( ana, topo )
    return None

  return addunit(float(r),'pb')


def getAxis(w,a):
  """For w=x,y, find according index in the masses-list, using the axes-information a=(mx-my)."""
  ml = []
  ml.append(a.find('M1'))
  ml.append(a.find('M2'))
  ml.append(a.find('M0'))
  if w=='x':
    return getIndex(ml,second=True)
  if w=='y':
    return getIndex(ml)
  

def compareMasses(masses, d):
  """Check if input masses are comparable to masses in the histogram corresponding to the information given in axes-dictionary d."""
  from tools.physicsUnits import rmvunit
  try: #check if histogram axes are M1, M0, return 1 if x-value of histogram is comparable to x value for given masses, 0 if not
    x1=getxval(masses[0],masses[-1],d['mz'][0])
    x2=float(rmvunit(masses[1],"GeV")-rmvunit(masses[-1],"GeV"))/(rmvunit(masses[0],"GeV")-rmvunit(masses[-1],"GeV"))
    if abs(x1-x2)/(x1+x2)<0.1:
      return True
    else: return None
  except:
    if d['mz'][0].find('LSP')>-1: #check if histogram for fixed LSP mass, return 1 if my is comparable to LSP mass of the histogram, 0 if not
      mlsp = float(d['mz'][0].split("LSP")[1])
      if abs(mlsp-rmvunit(masses[-1],"GeV"))/(mlsp+rmvunit(masses[-1],"GeV"))<0.1:
        return True
      else: return None
    elif d['mz'][0].find('D')>-1: #check for fixed deltaM
      ml = []
      ml.append(d['mz'][0].find('M1'))
      ml.append(d['mz'][0].find('M2'))
      ml.append(d['mz'][0].find('M0'))
      deltam = float(d['mz'][0].split('=')[1])
      deltain = rmvunit(masses[getIndex(ml,second=True)],"GeV")-rmvunit(masses[getIndex(ml)],"GeV")
      if deltain<0: return None
      if abs(deltain-deltam)/(deltain+deltam)<0.1: return True
      else: return None
    elif d['mz'][0].find('M1')>-1: #check for fixed m_mother
      mmother=float(d['mz'][0].split("M1")[1])
      if abs(mmother-rmvunit(masses[0],"GeV"))/(mmother+rmvunit(masses[0],"GeV"))<0.1:
        return True
      else: return None

def getIndex(ls, second=False):
  """Find index of list element with maximum value. If the last element is the maximum, return -1.
      If second=True, find second largest list element."""
  ind = np.argmax(ls)
  if second:
    v=0
    for i in range(len(ls)):
      if not i==ind:
        if ls[i]>=v:
          v=ls[i]
          ind2=i
    if ind2==len(ls)-1: return -1
    return ind2
  if ind==len(ls)-1: return -1
  return ind
    

def getHistName(topo, mz):
  """Build histogram name for given topology and mz information (mz as given in the axes-information)."""
  if mz == None:
    return topo
  elif 'D' in mz:
    return topo+'D'+mz.split('=')[1]
  else: return topo+mz

def upperLimit(ana, topo, masses,debug=True,run=None):
  """Returns upper limit for ana, topo, for given masses. masses: list of masses, with (mother, intermediate(s), LSP).
      For intermediate masses: if possible do interpolation over upper limits for different x-values.
      If interpolation is not possible: check if masses are comparable to the assumptions in the histogram."""
  d = smsResults.getaxes(ana, topo)
  if not run: run=smsHelpers.getRun(ana)
  if not d:
    if debug: print "[smsInterpolation] error: %s/%s not found." % ( ana, topo )
    return None
  if len(masses)==2 and not d[0]['mz']:
    return smsResults.getUpperLimit(ana, topo, masses[getAxis('x',d[0]['axes'])], masses[getAxis('y',d[0]['axes'])], interpolate=True)
  if len(masses)==2 and d[0]['mz']:
    if debug: print "[smsInterpolation] error: Need intermediate mass input for %s/%s." % ( ana, topo )
    return None
  if len(masses)>2 and not d[0]['mz']:
    if debug: print "[smsInterpolation] error: No intermediate mass in %s/%s." % ( ana, topo )
    return None
  if len(masses)>3 or len(d[0]['mz'])>1:
    if debug: print "[smsInterpolation] error: More than one intermediate mass in %s/%s. Cannot find upper limit for topologies with more than one intermediate mass." % ( ana, topo )
    return None
  if len(masses)>2 and len(d) == 1:
    if compareMasses(masses,d[0]):
      return smsResults.getUpperLimit(ana, getHistName(topo,d[0]['mz'][0]),masses[getAxis('x',d[0]['axes'])],masses[getAxis('y',d[0]['axes'])],interpolate=True)
    if debug: print "[smsInterpolation] error: Only one histogram available for %s/%s, cannot interpolate for intermediate mass." % ( ana, topo )
    return None
  return doGridData(ana,topo,masses,d,debug,run)
