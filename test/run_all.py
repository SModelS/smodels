#!/usr/bin/python

""" Systematically run all scripts and check if we get an error """

import os, sys
from TestTools import ok

Dir="../bin/"

scripts=os.listdir ( Dir )

print "Start to systematically check all scripts:\n"

for script in scripts:
  if script[-3:]!=".py": continue
  if script=="run_all.py": continue
  if script=="set_path.py": continue
  print "%s: " % ( script ),
  cmd= "python " + Dir + script + " > /dev/null"
  ret=os.system ( cmd )
  print ok ( 0, ret )

print "\nFinally check SMSmain.py: ",
sys.stdout.flush()
os.chdir("..")
tmpfile="/tmp/compare"+os.getlogin()
cmd= "python SMSmain.py > %s" % tmpfile
ret=os.system ( cmd )
print ok ( 0, ret )

os.chdir("test/")
import commands
out=commands.getoutput("diff SMSmain.log %s" % tmpfile)
print "Diff of the logs: %s" % ok ( "", out, False )
