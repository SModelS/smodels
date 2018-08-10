#!/usr/bin/env python

"""
.. module:: set efficiencies to zero when they do not satisfy the M_cut
"""


import sys
import os,glob
import argparse
import types

mcuts = []
for f in glob.glob('effmap*.dat'):
    ff = open(f,'r')
    data = ff.readlines()
    ff.close()
    mcuts.append(eval(f[f.find('mrec')+4:f.find('.dat')])/0.6)

mcuts = sorted(list(set(mcuts)))


for f in glob.glob('effmap*.dat'):

    ff = open(f,'r')
    data = ff.readlines()
    ff.close()
    mcut = eval(f[f.find('mrec')+4:f.find('.dat')])/0.6
    header = [x for x in data[0].split() if x and x.strip() != '#']
    mhscp_index = header.index('mhscp')
    newdata = [data[0]]
    for l in data[1:]:
        pt = [eval(x) for x in l.split()]
        mhscp = pt[mhscp_index]
        eff = pt[-1]
        if mhscp < mcut and eff:
            print 'Efficiency should be zero! (mhscp < 0.6*mreco)'
            print l
        if mcuts.index(mcut)+1 < len(mcuts):
            mcut_next = mcuts[mcuts.index(mcut)+1]
        else:
            mcut_next = None

        if mcut_next and mhscp > mcut_next and eff:
#            print 'Setting\n',l,'\nto zero'
            continue
        newdata.append(l)

    ff = open(f.replace('.dat','_clean.txt'),'w')
    for l in newdata:
        ff.write(l)
			
sys.exit()
