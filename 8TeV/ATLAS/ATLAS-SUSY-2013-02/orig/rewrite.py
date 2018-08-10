#!/usr/bin/python

a=open("limit_T5WWLSP060.txt")
lines=a.readlines()
a.close()

mchi0=60.

f=open("limit_T5WWLSP060.corr.txt","w")
for line in lines:
    values=map ( float, line.split () )
    print values
    mgluino=values[0]
    x=values[1]
    ul=values[2]
    mchip=x*(mgluino-mchi0)+mchi0
    f.write ( "%f %f %f\n" % ( mgluino, mchip, ul ) )

