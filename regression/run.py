#!/usr/bin/python

import os, commands, sys
from TestTools import ok

def run ( Nr ):
  File="test%d.py" % Nr
  print "%d:" % Nr,
  user=os.getlogin()
  logfile="/tmp/%s%d.log" % ( user, Nr )
  oldstdout=sys.stdout
  sys.stdout=open(logfile,"w")
  exec("import test%d" % Nr)
  sys.stdout.close()
  sys.stdout=oldstdout
  doc=eval("test%d.__doc__" % Nr)
  # os.system ( "python %s > %s" % ( File, logfile ) )
  cmd="diff %d.log %s" % ( Nr, logfile )
  out=commands.getoutput( cmd )
  ret=ok ( "", out, False )
  print doc,ret
  if ret.find("failed")>-1:
    print "try this:\n%s" % cmd
  else:
    os.unlink ( logfile )

if len(sys.argv)>1:
  for i in sys.argv[1:]:
    run ( int(i) )
  sys.exit(0)

Files=[]
for File in os.listdir("."):
  if File[-3:]!=".py": continue
  if File[:4]!="test": continue
  try:
    Nr=int(File[4:-3])
    Files.append(Nr)
  except Exception,e:
    pass
Files.sort()
for Nr in Files:
  run ( Nr )
