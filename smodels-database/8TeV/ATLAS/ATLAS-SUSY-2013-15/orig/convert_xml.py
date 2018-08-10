#!/usr/bin/python

f=open("exclusion_T6bbWWoffM1300.txt")
lines=f.readlines()
f.close()

for line in lines:
    if not line[:6]=="<point":
        continue
    dx=line.find("dx")
    line=line[dx+4:].replace("dy", " ").replace("' />","").replace("='","").replace("'","")
    print line.replace("\n","").replace("  "," ")
