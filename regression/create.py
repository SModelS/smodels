#!/usr/bin/python

import os, sys

def create ( Nr ):
  File="test%d.py" % Nr
  print "%d -> %s" % ( Nr, File ),
  os.system ( "python %s > %d.log" % ( File, Nr ) )
  oldstdout=sys.stdout
  sys.stdout=open("/dev/null","w")
  exec("import test%d" % Nr)
  sys.stdout.close()
  sys.stdout=oldstdout
  doc=eval("test%d.__doc__" % Nr)
  print doc

if len(sys.argv)>1:
  for i in sys.argv[1:]:
    create ( int(i) )
  sys.exit(0)

Files=[]
for File in os.listdir("."):
  if File[-3:]!=".py": continue
  if File[:4]!="test": continue
  Nr=int(File[4:-3])
  Files.append(Nr)
Files.sort()
for Nr in Files:
  create ( Nr )
