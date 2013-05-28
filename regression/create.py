#!/usr/bin/python

import os, sys

def create ( Nr ):
  File="test%d.py" % Nr
  print "%d -> %s" % ( Nr, File )
  os.system ( "python %s | tee %d.log" % ( File, Nr ) )

if len(sys.argv)>1:
  for i in sys.argv[1:]:
    create ( int(i) )
  sys.exit(0)

for File in os.listdir("."):
  if File[-3:]!=".py": continue
  if File[:4]!="test": continue
  Nr=int(File[4:-3])
  create ( Nr )
