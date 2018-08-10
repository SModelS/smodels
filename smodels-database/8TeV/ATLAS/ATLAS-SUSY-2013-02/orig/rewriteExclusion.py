#!/usr/bin/python

a=open("exclusion_T5WWLSP060.txt")
lines=a.readlines()
a.close()

mchi0=60.

f=open("exclusion_T5WWLSP060.corr.txt","w")
for line in lines:
    values=map ( float, line.split () )
    print values
    mgluino=values[0]
    x=values[1]
    mchip=x*(mgluino-mchi0)+mchi0
    f.write ( "%f %f \n" % ( mgluino, mchip  ) )

