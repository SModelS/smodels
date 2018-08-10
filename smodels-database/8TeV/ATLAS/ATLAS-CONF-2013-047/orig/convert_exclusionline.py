#!/usr/bin/python

orig = open("T6WWLSP060_orig_exc.dat")
corr = open("T6WWLSP060_exc.dat","write")

for l in orig.readlines():
    lo = l.split()
    if lo:
      #print lo[0]
      #print len(lo)
      try:
        first = float(lo[0])
        #if len(lo) < 3:
          #continue
        newM = float(lo[1])*float(lo[0])+(1-float(lo[1]))*60
        ln = lo[0]+"  "+str(newM)+"\n"
        corr.write(ln)
      except (ValueError):
        #print "lo[0] no mx"
        pass

orig.close()
corr.close()

