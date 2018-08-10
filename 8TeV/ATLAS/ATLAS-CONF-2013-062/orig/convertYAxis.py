#!/usr/bin/python

for fig in ["Fig14a", "Fig14b"]:

  orig = open("%sOrig.txt" %fig)
  corr = open("%s.txt" %fig, 'w')

  for l in orig.readlines():
    lo = l.split()
    newM = float(lo[1])*float(lo[0])+(1-float(lo[1]))*60
    ln = lo[0]+"  "+str(newM)+"   "+lo[2]+"\n"
    corr.write(ln)

orig.close()
corr.close()
