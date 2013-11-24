#!/usr/bin/python


infile = open('../tests.dat','r')
files = infile.readlines()
files = [ f.replace('\n','').split()[3] for f in files]
infile.close()
infile = open('bak.dat','r')
files2 = infile.readlines()
files2 = [ f.replace('\n','') for f in files2]
infile.close()

for f in files:
  fmatch = None
  for f2 in files2:
    try:
      if f2.replace('bak/','').split()[1]  == f:
        fmatch = f2.split()[0]
        break
    except: pass

  print f,fmatch

  newfile = open(f,'a')
  oldfile = open(fmatch,'r')
  olddata = oldfile.read()
  metadata = olddata[olddata.find('#END'):]  
  newfile.write(metadata)
  newfile.close()
  

  
