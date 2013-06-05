#!/usr/bin/python

import os, sys

if len(sys.argv)<2:
  print "Usage:",sys.argv[0],"<destination_directory>"
  sys.exit(0)

dest=sys.argv[1]
print "Installing the database to %s:" % dest

DB="/afs/hephy.at/user/w/walten/public/sms"

Dirs=[ "2011", "2012", "RPV7", "RPV8", "ATLAS8TeV","8TeV" ]

for Dir in Dirs:
  cmd="cp -r %s/%s %s" % (DB, Dir, dest)
  print cmd
  os.system ( cmd )
