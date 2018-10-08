#!/usr/bin/env python

"""
.. module:: convertToTxt
   :synopsis: used to convert the original HSCP efficiency maps to .txt tables

"""
import sys
import os


maps = ['./HSCPc%i00/maps.py'%i for i in range(4)]

#Dictionary to convert from original database to the arxiv model labeling convention
modelsDict = {'M1' : 'T0', 'M2' : 'T0M', 'M3' : 'T2','M4' : 'T2M', 'M5' : 'T6', 'M6' : 'T5M', 'M7' : 'T5', 'M8' : 'T1'}
modelsDict = dict([[val,key] for key,val in modelsDict.items()])

for ifile,mfile in enumerate(maps):    
    execfile(mfile)
    
    for oldtx in Dict:
        newtx = oldtx.replace(oldtx.replace('HSCP',''),modelsDict[oldtx.replace('HSCP','')])
        data = Dict[oldtx]
        ypts = None
        zpts = None
        if newtx == 'HSCPM1' or newtx == 'HSCPM2':
            xpts = [pt[0][0][0] for pt in data if pt[1] > 0.]
        elif newtx == 'HSCPM3' or newtx == 'HSCPM4':
            xpts = [pt[0][0][0] for pt in data if pt[1] > 0.]
            ypts = [pt[0][0][1] for pt in data if pt[1] > 0.]
        elif newtx == 'HSCPM5':
            xpts = [pt[0][0][0] for pt in data if pt[1] > 0.]
            ypts = [pt[0][0][1] for pt in data if pt[1] > 0.]
            zpts = [pt[0][0][2] for pt in data if pt[1] > 0.]
        elif newtx == 'HSCPM7' or newtx == 'HSCPM6':
            xpts = [pt[0][1][0] for pt in data if pt[1] > 0.]
            ypts = [pt[0][1][1] for pt in data if pt[1] > 0.]
            zpts = [pt[0][1][2] for pt in data if pt[1] > 0.]
        elif newtx == 'HSCPM8':
            xpts = [pt[0][0][0] for pt in data if pt[1] > 0.]
            ypts = [pt[0][0][1] for pt in data if pt[1] > 0.]
        elif newtx == 'HSCPM2':
            xpts = [pt[0][0] for pt in data if pt[1] > 0.]
            ypts = [pt[0][1] for pt in data if pt[1] > 0.]
            zpts = [pt[0][2] for pt in data if pt[1] > 0.]
        effpts = [pt[1] for pt in data if pt[1] > 0.]
            
        newMap = open('eff_%s_c%i00.txt'%(newtx,ifile),'w')
        newMap.write('## %s efficiencies\n' %newtx)
        header = "## x"
        if ypts: header += '\ty'
        if zpts: header += '\tz'
        header += '\t effi\n'
        newMap.write(header)
        for i,x in enumerate(xpts):
            line = str(x)
            if ypts: line += '\t%s' %str(ypts[i])
            if zpts: line += '\t%s' %str(zpts[i])
            line += "\t%s" %str(effpts[i])
            newMap.write(line+'\n')
        newMap.close()
                 
    
        
        
    
    
