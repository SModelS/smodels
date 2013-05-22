#SMSInterpolation is called by SMSResults.getSmartUpperLimit.
#UpperLimit takes arbitrary input masses and checks if there is a
#corresponding upper limit for the given analysis and topology.
#The upper limit is returned in 'pb'.
#If several histograms with different x-values are available,
#an interpolation is performed.

import SMSResults, numpy, SMSHelpers, copy
from SMSHelpers import rmvunit, addunit

def gethistname(topo, mz):
  """Build histogram name for given topology and mz information (mz as given in the axes-information)."""
  if mz == None or mz == "050":
    return topo
  elif 'D' in mz:
    return topo+'D'+mz.split('=')[1]
  else: return topo+mz

def getxval(mx, my, mz,mass=False):
  """Calculate x-value for one point.
      mx=Mother-mass, my=LSP-mass, mz is information on the intermediate mass as given in the axes-information.
      If mass=True is selected: return intermediate mass instead of x-value."""
  if mz.find('x')==-1 and mz.find('C')==-1 and mz.find('y')==-1:
    xfac = float(mz)/100
    if mass: return xfac*mx+(1-xfac)*my
    return float(mz)/100
  if mz.find('x')>-1: tx = float(mz[mz.find('x')+1:mz.find('x')+4])
  else: tx=None
#  if mz.find('y')>-1: ty = float(mz[mz.find('y')+1])
#  else: ty=None
  if mz.find('C')>-1: c = float(mz[mz.find('C')+1:mz.find('C')+4])
  else: c=None
  z = 0.
  if tx: z += tx*my/100
#  if ty: z += ty*mx
  if c: z+=c
  if mass: return z
  xval = (z-my)/(mx-my)
  return xval

def getindex(ls, second=False):
  """Find index of list element with maximum value. If the last element is the maximum, return -1.
      If second=True, find second largest list element."""
  ind = numpy.argmax(ls)
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
    
def getaxis(w,a):
  """For w=x,y, find according index in the masses-list, using the axes-information a=(mx-my)."""
  ml = []
  ml.append(a.find('M1'))
  ml.append(a.find('M2'))
  ml.append(a.find('M0'))
  if w=='x':
    return getindex(ml,second=True)
  if w=='y':
    return getindex(ml)
  
def compareM(masses, d):
  """Check if input masses are comparable to masses in the histogram corresponding to the information given in axes-dictionary d."""
  try: #check if histogram axes are M1, M0, return 1 if x-value of histogram is comparable to x value for given masses, 0 if not
    x1=getxval(masses[0],masses[-1],d['mz'][0])
    x2=float(masses[1]-masses[-1])/(masses[0]-masses[-1])
    if abs(x1-x2)/(x1+x2)<0.1:
      return True
    else: return None
  except:
    if d['mz'][0].find('LSP')>-1: #check if histogram for fixed LSP mass, return 1 if my is comparable to LSP mass of the histogram, 0 if not
      mlsp = float(d['mz'][0][d['mz'][0].find('P')+1:d['mz'][0].find('P')+4])
      if abs(mlsp-masses[-1])/(mlsp+masses[-1])<0.1:
        return True
      else: return None
    elif d['mz'][0].find('D')>-1: #check for fixed deltaM
      ml = []
      ml.append(d['mz'][0].find('M1'))
      ml.append(d['mz'][0].find('M2'))
      ml.append(d['mz'][0].find('M0'))
      deltam = float(d['mz'][0].split('=')[1])
      deltain = float(masses[getindex(ml,second=True)]-masses[getindex(ml)])
      if deltain<0: return None
      if abs(deltain-deltam)/(deltain+deltam)<0.1: return True
      else: return None


def UpperLimit(ana, topo, masses,debug=True):
  """Returns upper limit for ana, topo, for given masses. masses: list of masses, with (mother, intermediate(s), LSP).
      For intermediate masses: if possible do interpolation over upper limits for different x-values.
      If interpolation is not possible: check if masses are comparable to the assumptions in the histogram."""
  d = SMSResults.getaxes(ana, topo)
  if not d:
    if debug: print "[SMSInterpolation] error: %s/%s not found." % ( ana, topo )
    return None
  if len(masses)==2 and not d[0]['mz']:
    return SMSResults.getUpperLimit(ana, topo, masses[getaxis('x',d[0]['axes'])], masses[getaxis('y',d[0]['axes'])])
  if len(masses)==2 and d[0]['mz']:
    if debug: print "[SMSInterpolation] error: Need intermediate mass input for %s/%s." % ( ana, topo )
    return None
  if len(masses)>2 and not d[0]['mz']:
    if debug: print "[SMSInterpolation] error: No intermediate mass in %s/%s." % ( ana, topo )
    return None
  if len(masses)>3 or len(d[0]['mz'])>1:
    if debug: print "[SMSInterpolation] error: More than one intermediate mass in %s/%s. Cannot find upper limit for topologies with more than one intermediate mass." % ( ana, topo )
    return None
  if len(masses)>2 and len(d) == 1:
    if compareM(masses,d[0]):
      return SMSResults.getUpperLimit(ana, gethistname(topo,d[0]['mz'][0]),masses[getaxis('x',d[0]['axes'])],masses[getaxis('y',d[0]['axes'])])
    if debug: print "[SMSInterpolation] error: Only one histogram available for %s/%s, cannot interpolate for intermediate mass." % ( ana, topo )
    return None
  xsecs = []
  xvals = []
  d_it=copy.deepcopy(d)
  for ds in d_it:
    if not ds['mz']:
      if debug: print "[SMSInterpolation] error: No information on intermediate mass availabel for %s/%s." % ( ana, topo )
      return None
    if 'LSP' in ds['mz'][0] or 'D' in ds['mz'][0]:
      if compareM(masses,ds):
        if SMSResults.getUpperLimit(ana, gethistname(topo,ds['mz'][0]),masses[getaxis('x',ds['axes'])],masses[getaxis('y',ds['axes'])]):
          return SMSResults.getUpperLimit(ana, gethistname(topo,ds['mz'][0]),masses[getaxis('x',ds['axes'])],masses[getaxis('y',ds['axes'])])
      d.remove(ds)
      continue
    xs=rmvunit(SMSResults.getUpperLimit(ana, gethistname(topo,ds['mz'][0]),masses[getaxis('x',ds['axes'])],masses[getaxis('y',ds['axes'])]),'pb')
    if xs:
      xsecs.append(xs)
      xvals.append(getxval(masses[0],masses[-1],ds['mz'][0]))
    else:
      d.remove(ds)
  if len(xsecs)<2:
    if d:
      if compareM(masses,d[0]): return SMSResults.getUpperLimit(ana, gethistname(topo,d[0]['mz'][0]),masses[getaxis('x',d[0]['axes'])],masses[getaxis('y',d[0]['axes'])])
    if debug: print "[SMSInterpolation] error: Available histograms for %s/%s could not be combined, cannot interpolate" % ( ana, topo )
    return None
  p = numpy.polyfit(xvals,xsecs,len(xsecs)-1)
  X = float(masses[1]-masses[-1])/float(masses[0]-masses[-1])
  XSec = float(numpy.polyval(p,X))
  if X>max(xvals) or X<min(xvals):
    if debug: print "[SMSInterpolation] error: x value for %s/%s out of range, no extrapolation" % ( ana, topo )
    return None
#  dix= max(xvals)-min(xvals)
#  diy=max(xsecs)-min(xsecs)
#  if (X>max(xvals) and X-max(xvals)>dix) or (X<min(xvals) and min(xvals)-X>dix):
#    print "x value for intermediate mass out of extrapolation range"
#    return 0
#  if (XSec>max(xsecs) and XSec-max(xsecs)>diy) or (XSec<min(xsecs) and min(xsecs)-XSec>diy):
#    print "xsec value out of extrapolation range"
#    return 0
#  if XSec<0:
#    print "Negative crosssection, cannot extrapolate"
#    return 0
  return addunit(XSec,'pb')
