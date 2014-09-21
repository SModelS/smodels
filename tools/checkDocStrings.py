#!/usr/bin/env python

"""
.. module:: tools.checkDocStrings.py
   :synopsis:  Goes over all .py files in the parent folder (and subfolders)
   and check for mis-matched code-docstrings or docstrings that are too short.
"""

import fnmatch
import os,glob


minDocStringLength = 20

matches = []

for f in glob.glob("../*"):
    for root, dirnames, filenames in os.walk(f):
        for filename in fnmatch.filter(filenames, '*.py'):
            matches.append(os.path.join(root, filename))
        
tagLinesWith = ['def ','class ']
ignoreLinesWith = ['__init__(self','__eq__(self','__ne__(self','__str__(self',
                   '__iter__(self','__getitem__(self','__mul__(self',
                   '__rmul__(self','__add__(self','__setitem__(self',
                   '__len__(self','__repr__(self','__cmp__(self']

info = {}
badBlocks = {}
for fname in sorted(matches):
    if 'checkDocStrings' in fname: continue #Ignore this file
    info[fname] = {'docstrings' : [], 'codesnippet' : []}
    badBlocks[fname] = {'short docstring' : [], 'missing docstring' : []}
    fcheck = open(fname,'r')
    fdata = fcheck.read().lstrip()
    if fdata[0] == '"""': idoc = 'evenBlocks'
    else: idoc = 'oddBlocks'
    fdata = fdata.split('"""')
    if len(fdata) < 2:
        info[fname]['docstrings'] = []
    else:
        for iblock,block in enumerate(fdata):
            pblock = 'odd'
            if iblock%2 == 0: pblock = 'even'
            if pblock in idoc: #Text block is a docstring
                info[fname]['docstrings'].append(block)
                if len(block.replace(' ',"")) < minDocStringLength: badBlocks[fname]['short docstring'].append(True)
                else: badBlocks[fname]['short docstring'].append(False)
            else: #Text block is regular code
                codesnippet = ""
                nTag = 0;
                nIgTag = 0;
                for l in block.split('\n'):
                    ignoreLine = False
                    for tag in ignoreLinesWith:
                        if tag in l: ignoreLine = True
                    if ignoreLine:
                        nIgTag += 1
                        continue
                    nTag += sum([l.count(tag) for tag in tagLinesWith])          
                if nTag != 1 and nIgTag != 1:
                    badBlocks[fname]['missing docstring'].append(True)
                    doPrint = True   #Print blocks with more than one (class or method)/docstring
                    codesnippet += block
                    info[fname]['codesnippet'].append(codesnippet)
                else:
                    badBlocks[fname]['missing docstring'].append(False)
                    lines = block.split('\n')
                    tagFound = False              
                    while not tagFound and len(lines) > 0:
                        line = lines.pop(-1)
                        if not line.replace('\n', '').strip(): continue  #ignore empty lines                      
                        codesnippet = line + "\n" + codesnippet                    
                        for tag in tagLinesWith:
                            if tag in line: tagFound = True
                    if not codesnippet: codesnippet = "NOT FOUND"
                    info[fname]['codesnippet'].append(codesnippet)

nmissing = 0
nshort = 0
badFiles = []
for fname in sorted(info.keys()):
    finfo = info[fname]    
    print '-----------File:',fname
    if not finfo['docstrings']:
        print 'No doc strings'
        continue
    for idoc,doc in enumerate(finfo['docstrings']):
        if not badBlocks[fname]['short docstring'][idoc] and not badBlocks[fname]['missing docstring'][idoc]: continue
        elif badBlocks[fname]['missing docstring'][idoc]:
            if "module::" in doc: continue #Ignore module descriptions
            if 'if__name__=="__main__":' in finfo['codesnippet'][idoc].replace(" ",''): continue #Ignore main descriptions
            nmissing += 1
            if not fname in badFiles: badFiles.append(fname)     
            print '---Code snippet (missing doc string):',finfo['codesnippet'][idoc]
            print '---Doc string:',doc            
        elif badBlocks[fname]['short docstring'][idoc]:
            nshort += 1
            if not fname in badFiles: badFiles.append(fname)
            print '---Code snippet:',finfo['codesnippet'][idoc]
            print '---Doc string (too short):',doc
            
print '\n\n\n'
print 'Number of files = ',len(info)
print 'Number of bad files = ',len(badFiles)
print 'Number of docstrings = ',sum([len(info[fname]['docstrings']) for fname in info.keys()])
print 'Number of missing docstrings = ',nmissing
print 'Number of short docstrings = ',nshort
print '\n Bad files:\n'
for f in sorted(badFiles): print f
            
