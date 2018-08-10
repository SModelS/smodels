#!/usr/bin/python

a=open("T6WWLSP060.data")
lines=a.readlines()
a.close()

mchi0=60.

f=open("T6WWLSP060corr.data","w")
for line in lines:
    if line.find("#")==0:
        continue
    values=map ( float, line.split () )
    mgluino=values[0]
    x=values[1]
    ul=values[2]
    mchip=x*(mgluino-mchi0)+mchi0
    print mgluino,mchip
    f.write ( "%f %f %f\n" % ( mgluino, mchip, ul ) )

f.close()
