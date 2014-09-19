#!/usr/bin/env python

import fnmatch
import os,glob,sys

matches = []

for f in glob.glob("*"):
    for root, dirnames, filenames in os.walk(f):
        for filename in fnmatch.filter(filenames, '*.py'):
            matches.append(os.path.join(root, filename))
        
codesnippetTags = ['def','class']
info = {}
for fname in matches:
    info[fname] = {'docstrings' : [], 'codesnippet' : []}
    fcheck = open(fname,'r')
    fdata = fcheck.read().lstrip()
    if fdata[0] == '"""': idoc = 'evenBlocks'
    else: idoc = 'oddBlocks'
    fdata = fdata.split('"""')
    if len(fdata) < 2:
        info[fname]['docstrings'] = None
    else:
        for iblock,block in enumerate(fdata):
            pblock = 'odd'
            if iblock%2 == 0: pblock = 'even'
            if pblock in idoc: #Text block is a docstring
                info[fname]['docstrings'].append(block)
            else: #Text block is regular code
                  codesnippet = ""
                  lines = block.split('\n')
                  tagFound = False                  
                  while not tagFound and len(lines) > 0:
                    line = lines.pop(-1)
                    if not line.replace('\n', '').strip(): continue  #ignore empty lines                      
                    codesnippet = line + "\n" + codesnippet                    
                    for tag in codesnippetTags:
                        if tag in line: tagFound = True
                  if not codesnippet: codesnippet = "NOT FOUND"
                  info[fname]['codesnippet'].append(codesnippet)

for fname,finfo in info.items():    
    print '-----------File:',fname
    if not finfo['docstrings']:
      print 'No doc strings'
      continue
    for idoc,doc in enumerate(finfo['docstrings']):
        print '---Code snippet:',finfo['codesnippet'][idoc]
        print '---Doc string:',doc
