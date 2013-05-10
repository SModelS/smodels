#!/usr/bin/python

import SMSResults, numpy, SMSHelpers
from SMSHelpers import rmvunit, addunit

def gethistname(topo, mz):
	if mz == None or mz == "050":
		return topo
	else: return topo+mz

def getxval(mx, my, mz):
	if mz.find('x')==-1 and mz.find('C')==-1 and mz.find('y')==-1:
		return float(mz)/100
	if mz.find('x')>-1: tx = float(mz[mz.find('x')+1])
	else: tx=None
	if mz.find('y')>-1: ty = float(mz[mz.find('y')+1])
	else: ty=None
	if mz.find('C')>-1: c = float(mz[mz.find('C')+1:mz.find('C')+4])
	else: c=None
	z = 0.
	if tx: z += tx*my
	if ty: z += ty*mx
	if c: z+=c
	xval = (z-my)/(mx-my)
	return xval

def getindex(ls, second=False):
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
	ml = []
	ml.append(a.find('M1'))
	ml.append(a.find('M2'))
	ml.append(a.find('M0'))
	if w=='x':
		return getindex(ml,second=True)
	if w=='y':
		return getindex(ml)
	
def compareM(masses, d):
	try: #check if histogram axes are M1, M0, return 1 if x-value of histogram is comparable to x value for given masses, 0 if not
		x1=getxval(masses[0],masses[-1],d['mz'][0])
		x2=float(masses[1]-masses[-1])/(masses[0]-masses[-1])
		print x1,x2
		if abs(x1-x2)/(x1+x2)<0.1:
			return True
		else: return 0
	except:
		if d['mz'][0].find('LSP')>-1: #check if histogram for fixed LSP mass, return 1 if my is comparable to LSP mass of the histogram, 0 if not
			mlsp = float(d['mz'][0][d['mz'][0].find('P')+1:d['mz'][0].find('P')+4])
			if abs(mlsp-masses[-1])/(mlsp+masses[-1])<0.1:
				return True
			else: return 0
		elif d['mz'][0].find('D')>-1: #check for fixed deltaM
			ml = []
			ml.append(d['mz'][0].find('M1'))
			ml.append(d['mz'][0].find('M2'))
			ml.append(d['mz'][0].find('M0'))
			deltam = d['mz'][0].split('=')[1]
			deltain = masses[getindex(ml,second=True)]-masses[getindex(ml)]
			if deltain<0: return 0
			if abs(deltain-deltam)/(deltain+deltam)<0.1: return True
			else: return 0


def UpperLimit(ana, topo, masses): #masses: list of masses, with (mother, intermediate(s), LSP)
	d = SMSResults.getaxes(ana, topo)
	print d
	if len(masses)==2 and not d[0]['mz']:
		return SMSResults.getUpperLimit(ana, topo, masses[getaxis('x',d[0]['axes'])], masses[getaxis('y',d[0]['axes'])])
	if len(masses)==2 and d[0]['mz']:
		print "Need intermediate mass input for this topology."
		return 0
	if len(masses)>2 and not d[0]['mz']:
		print "No intermediate mass for this Topology."
		return 0
	if len(masses)>2 and len(d) == 1:
		if compareM(masses,d[0]):
			return SMSResults.getUpperLimit(ana, gethistname(topo,d[0]['mz'][0]),masses[getaxis('x',d[0]['axes'])],masses[getaxis('y',d[0]['axes'])])
		print "Only one histogram available, cannot interpolate for intermediate mass."
		return 0
	xsecs = []
	xvals = []
	for ds in d:
		if not ds['mz']:
			print "No information on intermediate mass availabel."
			return 0
		if 'LSP' in ds['mz'][0] or 'D' in ds['mz'][0]:
			if compareM(masses,ds):
				return SMSResults.getUpperLimit(ana, gethistname(topo,ds['mz'][0]),masses[getaxis('x',ds['axes'])],masses[getaxis('y',ds['axes'])])
			else:
				d.remove(ds)
				continue
		xsecs.append(rmvunit(SMSResults.getUpperLimit(ana, gethistname(topo,ds['mz'][0]),masses[getaxis('x',ds['axes'])],masses[getaxis('y',ds['axes'])]),'pb'))
		xvals.append(getxval(masses[0],masses[-1],ds['mz'][0]))
	if len(xsecs)<2:
		if compareM(masses,d[0]): return SMSResults.getUpperLimit(ana, gethistname(topo,d[0]['mz'][0]),masses[getaxis('x',d[0]['axes'])],masses[getaxis('y',d[0]['axes'])])
		print "Available histograms could not be combined, cannot interpolate"
		return 0
	print xsecs, xvals
	p = numpy.polyfit(xvals,xsecs,len(xsecs)-1)
	X = float(masses[1]-masses[-1])/float(masses[0]-masses[-1])
	XSec = float(numpy.polyval(p,X))
	print X
	if X>max(xvals) or X<min(xvals):
		print "x value out of range, no extrapolation"
		return 0
#	dix= max(xvals)-min(xvals)
#	diy=max(xsecs)-min(xsecs)
#	if (X>max(xvals) and X-max(xvals)>dix) or (X<min(xvals) and min(xvals)-X>dix):
#		print "x value for intermediate mass out of extrapolation range"
#		return 0
#	if (XSec>max(xsecs) and XSec-max(xsecs)>diy) or (XSec<min(xsecs) and min(xsecs)-XSec>diy):
#		print "xsec value out of extrapolation range"
#		return 0
#	if XSec<0:
#		print "Negative crosssection, cannot extrapolate"
#		return 0
	return addunit(XSec,'pb')
