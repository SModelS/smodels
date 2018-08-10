#!/usr/bin/python

f = open('names.txt','r')
lines = f.readlines()
data = {}
defoults = ['link','file','commeny']
for l in lines:
    if l[:1] == 'T':
        topo = l.split()[0]
        data[topo] = {}
        continue
    if not l.split()==[] and  not 'Fig' in l and not 'fig' in l:
        obj = l.split()[0]
        print 
        fileName = l.split()[1] 
        if not fileName in defoults: 
            data[topo][obj] = fileName
print data

for topo in data:
    if 'limit' in data[topo]:
        print "histo"