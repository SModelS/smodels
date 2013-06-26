#!/usr/bin/python

import os, sys

def usage():
  print "Usage:",sys.argv[0]," [-h] [-scp] <destination_directory>"
  print "        -scp: use scp to smodels instead of local cp"
  print "        -h: show this help"
  sys.exit(0)

if len(sys.argv)<2:
  usage()

useScp=False

for i in sys.argv[1:]:
  if i=="-scp": useScp=True
  if i=="-h": usage()

dest=sys.argv[-1]
print "Installing the database to %s:" % dest

if not os.path.exists ( dest ):
  os.mkdir ( dest )

DB="/afs/hephy.at/user/w/walten/public/sms"

Dirs=[ "2011", "2012", "RPV7", "RPV8", "ATLAS8TeV","8TeV" ]

if useScp:
  for Dir in Dirs:
    cmd="scp -r smodels.hephy.at:%s/%s %s " % (DB, Dir, dest)
    print cmd
    os.system ( cmd )
else:
  cmd="cp -r"
  for Dir in Dirs:
    cmd+=" %s/%s " % (DB, Dir)
  cmd+= dest
  print cmd
  os.system ( cmd )
