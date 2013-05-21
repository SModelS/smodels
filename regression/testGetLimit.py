#!/usr/bin/python

import sys
sys.path.append ( "../" )

import ROOT, SMSmethods, mySMSgetlimit, SMSResults, SMSInterpolation
from SMSHelpers import addunit, rmvunit

def recreateHist(ana,topo,mz=None, run=''):
	""" recreate ROOT TH2F histogram of a given analysis and topology (only topologies with no intermediate masses) """
	if SMSResults.hasDictionary(ana):
		print "Cannot recreate histogram from dictionary"
		return None	
#	xmin=rmvunit(SMSResults.getLowX(ana,topo),'GeV')
#	ymin=rmvunit(SMSResults.getLowY(ana,topo),'GeV')
#	xmax=rmvunit(SMSResults.getUpX(ana,topo),'GeV')
#	ymax=rmvunit(SMSResults.getUpY(ana,topo),'GeV')
#	bwx=rmvunit(SMSResults.getBinWidthX(ana,topo),'GeV')
#	bwy=rmvunit(SMSResults.getBinWidthY(ana,topo),'GeV')
#	bx=int((xmax-xmin)/bwx)
#	by=int((ymax-ymin)/bwy)
#	print bx,by
	xmin=0.
	xmax=1000.
	bwx=20.
	ymin=0.
	ymax=1000.
	bwy=20.
	bx=int((xmax-xmin)/bwx)
        by=int((ymax-ymin)/bwy)
	h = ROOT.TH2F('h','h',bx,xmin,xmax,by,ymin,ymax)
	x=xmin+bwx/2
	y=ymin+bwy/2
	a=SMSmethods.EAnalysis()
	while x<xmax:
		while y<ymax:
			massv = [x]
			if mz: massv.append(SMSInterpolation.getxval(x,y,mz,mass=True))
			massv.append(y)
			v=rmvunit(mySMSgetlimit.GetPlotLimit([massv,massv],[topo,[ana]],a)[0][1],'pb')
			if v: h.Fill(x,y,v)
			y+=bwy
		y=ymin+bwy/2
		x+=bwx
	return h
