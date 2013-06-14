#!/usr/bin/env python

"""
.. module:: TxNames
    :synopsis: missing

.. moduleauthor:: missing <email@example.com>

"""

def getTx(element):
	"""takes EElement and returns Tx-Name"""
	import copy
	E=copy.deepcopy(element)
	Branch1=E.B[0]
	Branch2=E.B[1]
	d=E.getEinfo()
	v1=d["vertnumb"][0]
	v2=d["vertnumb"][1]
	p1=d["vertparts"][0]
	p2=d["vertparts"][1]
#Sort branches
	if v1 < v2 or (v1==v2  and sum(p1)<sum(p2)):
		Branch1, Branch2 = Branch2, Branch1
		v1, v2 = v2, v1
		p1, p2 = p2, p1
#read branch infortmation
	b1=[ v1,p1,Branch1.particles ]
	b2=[ v2,p2,Branch2.particles ]
#find SMS
	if b1[0]==1 or b2[0]==1:
		return "bad input"
	if b1[1]==[2,0] and b2[1]==[2,0]:
		return getT1(ptsCount(b1[2],b2[2]))
	if (b1[1]==[2,2,0] or b1[1]==[2,1,0]) and b2[1]==[2,0]:
		return getT3(ptsCount(b1[2],b2[2]), b2)
	if (b1[1]==[2,2,0] and b2[1]==[2,2,0]) or (b1[1]==[2,1,0] and b2[1]==[2,1,0]):
       		return getT5(ptsCount(b1[2],b2[2]))
	if b1[1]==[1,0] and b2[1]==[1,0]:
                return getT2(ptsCount(b1[2],b2[2]))
	if b1[1]==[1,1,0] and b2[1]==[1,0]:
		return getT4(ptsCount(b1[2],b2[2]))
        if b1[1]==[1,1,0] and b2[1]==[1,1,0]:
                return getT6(ptsCount(b1[2],b2[2]))
	return "?"

def ptsCount(b1,b2):
	pts=[]
	b1.extend(b2)
	for b in b1:
		pts.extend(b)
	return pts

def getT1(pts):
	if pts.count('b')==4:
		return "T1bbbb"
	elif pts.count('t+')+pts.count('t-')==4:
		return "T1tttt"
	elif pts.count('b')==2 and pts.count('t+')+pts.count('t-')==2:
		return "T1tbtb"
	else: return "T1"

def getT3(pts, b2):
        if pts.count('W+')+pts.count('W-')==1 and b2[2][0].count('b')==2:
                return "T3Wb"
        elif pts.count('W+')+pts.count('W-')==1:
                return "T3W"
	elif pts.count('Z')==1:
		return "T3Z"
	elif (pts.count('e+')==1 and pts.count('e-')==1) or (pts.count('mu+')==1 and pts.count('mu-')==1):
		return "T3lh"
	elif (pts.count('e+')+pts.count('e-')==1 or pts.count('mu+')+pts.count('mu-')==1) and pts.count('nu')==1:
		return "T3lnu"
        else: return "T3?"

def getT5(pts):
        if pts.count('photon')==2:
                return "T5gg"
        elif pts.count('W+')+pts.count('W-')==2:
                return "T5WW"
	elif pts.count('W+')+pts.count('W-')==1 and pts.count('Z')==1:
		return "T5WZ"
	elif pts.count('Z')==2:
		return "T5ZZ"
	elif pts.count('l+')+pts.count('l-')==2 and pts.count('nu')==2:
		return "T5lnu"
	elif pts.count('t+')+pts.count('t-')==4:
		return "T5tttt"
        else: return "T5?"

def getT2(pts):
        if pts.count('b')==2:
                return "T2bb"
	elif pts.count('b')==1 and pts.count('W+')+pts.count('W-')==1:
		return "T2bW" #gibt es so nicht, sollte T6bbww geben?
        elif pts.count('t+')==1 and pts.count('t-')==1:
                return "T2tt"
        elif pts.count('e+')+pts.count('e-')==2 or pts.count('mu+')+pts.count('mu-')==2:
                return "TSlepSlep"
        elif pts.count('W+')+pts.count('W-')==1 and pts.count('Z')==1:
                return "TChiWZ"
        elif pts.count('W+')==1 and pts.count('W-')==1:
                return "TChiWW"
        elif pts.count('Z')==2:
                return "TChiZZ"

	else: return "T2"

def getT4(pts):
	return "T4"

def getT6(pts):
        if pts.count('b')==2 and pts.count('Z')==2:
                return "T6bbZZ"
	elif pts.count('l+')+pts.count('l-')==3 and pts.count('nu')==1:
                return "TChiNuSlep" #same as TChiChipmSnuSlep?
        elif pts.count('l+')+pts.count('l-')==2 and pts.count('nu')==2:
                return "TChipmSnuSlep"
	elif pts.count('t')==2 and pts.count('W')==2:
		return "T6ttWW"
	elif pts.count('b')==2 and pts.count('W')==2:
		return "T6bbWW"
	elif pts.count('t')==2 and pts.count('Z')==2:
		return "T6ttZZ"
        else: return "T6?"

