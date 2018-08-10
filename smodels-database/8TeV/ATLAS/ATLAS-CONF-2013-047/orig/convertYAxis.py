#!/usr/bin/python

for fig in ["T5WWLSP060", "T6WWLSP060"]:

  orig = open("%sOrig.dat" %fig)
  corr = open("%s.dat" %fig, 'w')

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
	ln = lo[0]+"  "+str(newM)+"   "+lo[2]+"\n"
	corr.write(ln)
      except (ValueError):
	#print "lo[0] no mx"
	pass
              
    

orig.close()
corr.close()

