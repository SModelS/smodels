#!/usr/bin/python

""" a simple script that downloads the results database to a target directory, either locally via 'cp' (needs an afs installation on this machine), or via scp to smodels """

import os, sys

def usage():
  print "Usage:",sys.argv[0]," [-h] [-r] [-scp] <destination_directory>"
  print "        -scp: use scp to smodels instead of local cp"
  print "        -h: show this help"
  print "        -r: remove old local database, if exists"
  sys.exit(0)

if len(sys.argv)<2:
  usage()

useScp=False
force=False

for i in sys.argv[1:]:
  if i=="-scp": useScp=True
  if i=="-h": usage()
  if i=="-r": force=True

dest=sys.argv[-1]
print "Installing the database to %s:" % dest

if not os.path.exists ( dest ):
  os.mkdir ( dest )

DB="/afs/hephy.at/user/w/walten/public/sms"

Dirs=[ "2011", "2012", "RPV7", "RPV8", "ATLAS8TeV","8TeV" ]

if useScp:
  for Dir in Dirs:
    Target="%s/%s" % (dest, Dir)
    if os.path.exists ( Target ):
      print "Warning:",Target,"exists already."
      if force:
        print "Requested removal of",Target
        os.system ( "rm -rf %s" % Target )
    cmd="scp -r smodels.hephy.at:%s/%s %s " % (DB, Dir, dest)
    print cmd
    os.system ( cmd )
else:
  cmd="cp -r"
  for Dir in Dirs:
    Target="%s/%s" % (dest, Dir)
    print "Dir",Target
    if os.path.exists ( Target ):
      print "Warning:",Target,"exists already."
      if force:
        print "Requested removal of",Target
        os.system ( "rm -rf %s" % Target )
    cmd+=" %s/%s " % ( DB, Dir )
  cmd+= dest
  print cmd
  os.system ( cmd )
