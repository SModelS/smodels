#!/usr/bin/python

import os, commands, sys
from TestTools import ok

def run ( Nr ):
  File="test%d.py" % Nr
  print "%s" % ( File ),
  user=os.getlogin()
  logfile="/tmp/%s%d.log" % ( user, Nr )
  os.system ( "python %s > %s" % ( File, logfile ) )
  cmd="diff %d.log %s" % ( Nr, logfile )
  out=commands.getoutput( cmd )
  ret=ok ( "", out, False )
  os.unlink ( logfile )
  print ret

if len(sys.argv)>1:
  for i in sys.argv[1:]:
    run ( int(i) )
  sys.exit(0)

Files=[]
for File in os.listdir("."):
  if File[-3:]!=".py": continue
  if File[:4]!="test": continue
  Nr=int(File[4:-3])
  Files.append(Nr)
Files.sort()
for Nr in Files:
  run ( Nr )
