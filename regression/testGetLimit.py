#!/usr/bin/python

import sys
sys.path.append ( "../" )

import ROOT, SMSmethods, SMSgetlimit, SMSResults
from SMSHelpers import addunit, rmvunit

def recreateHist(ana,topo, run=''):
	""" recreate ROOT TH2F histogram of a given analysis and topology (only topologies with no intermediate masses) """
	if SMSResults.hasDictionary(ana):
		print "Cannot recreate histogram from dictionary"
		return None
	xmin=rmvunit(SMSResults.getLowX(ana,topo),'GeV')
	ymin=rmvunit(SMSResults.getLowY(ana,topo),'GeV')
	xmax=rmvunit(SMSResults.getUpX(ana,topo),'GeV')
	ymax=rmvunit(SMSResults.getUpY(ana,topo),'GeV')
	bwx=rmvunit(SMSResults.getBinWidthX(ana,topo),'GeV')
	bwy=rmvunit(SMSResults.getBinWidthY(ana,topo),'GeV')
	bx=int((xmax-xmin)/bwx)
	by=int((ymax-ymin)/bwy)
	print bx,by
	h = ROOT.TH2F('h','h',bx,xmin,xmax,by,ymin,ymax)
	x=xmin+bwx/2
	y=ymin+bwy/2
	a=SMSmethods.EAnalysis()
	while x<xmax:
		while y<ymax:
			v=rmvunit(SMSgetlimit.GetPlotLimit([[x,y],[x,y]],[topo,[ana]],a)[0][1],'pb')
			if v: h.Fill(x,y,v)
			y+=bwy
		y=ymin+bwy/2
		x+=bwx
	return h
