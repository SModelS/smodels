#!/usr/bin/python

f=open("T6bbWWx200.txt")
lines=f.readlines()
f.close()

g=open("T6bbWWx200_corr.txt","write")
for line in lines:
    tokens=line.split()
    g.write ("%s %s %s\n" % ( tokens[1], tokens[0], tokens[2] ) )
g.close()
