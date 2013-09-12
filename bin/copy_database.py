#!/usr/bin/python

""" a simple script that downloads the results database to a target directory, either locally via 'cp' (needs an afs installation on this machine), or via scp to smodels """

import os, sys

def usage():
  print "Usage:",sys.argv[0]," [-h] [-r] [-scp] <destination_directory>"
  print "        -scp: use scp to smodels instead of local cp (from afs)"
  print "        -h: show this help"
  print "        -r: remove old local database, if exists"
  sys.exit(0)

def localCopy ( dest, Dirs, force ):
  cmd="cp -r"
  for Dir in Dirs:
    Target="%s/%s" % (dest, Dir)
    print "Dir",Target
    if os.path.exists ( Target ):
      print "Warning:",Target,"exists already."
      if force:
        print "Requested removal of",Target
        os.system ( "rm -rf %s" % Target )
  #  if not os.path.exists ( Target ):
  #    os.mkdir ( Target )
    cmd+=" %s/%s " % ( DB, Dir )
  cmd+= dest
  print cmd
  os.system ( cmd )
  stripDatabase ( dest, Dirs )

def stripAnalysis ( path ):
  Files=os.listdir ( path )
  print "stripping",path,Files
  for F in Files:
    if not F in [ "sms.py", "sms.root", "info.txt" ]:
    ## if F in [ "orig", "old", "convert.py", "draw.py", "Standardizer.py", "convert.py~", "#Standardizer.py#", "info.txt_", "info.old", "Standardizer.pyc", "#convert.py#", "results" ]:
      cmd="rm -rf %s/%s" % ( path, F )
      print cmd
      os.system ( cmd )

def stripDatabase ( dest, Dirs ):
  for Dir in Dirs:
    path=dest + "/" + Dir
    anas= os.listdir ( path )
    for ana in anas:
      if ana.lower() in [ "todo", "readme", "old", "bad" ]: continue
      if ana[0]==".": continue
      if ana[-3:] in [ ".py", ".sh" ]: continue
      stripAnalysis ( path+"/"+ana )

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
  localCopy ( dest, Dirs, force )
