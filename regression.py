#!/usr/bin/python

green='\033[0;32m'
red='\033[0;31m'
reset='\033[;0m'

import lheDecomposer, sys

Dir="regression"

def verify ( Dict, SMS ):
  if not Dict.has_key ( SMS ): return False
  if Dict[SMS]!=1: return False
  return True

def ok ( A,B ):
  if A==B: return "%sok.%s" % ( green, reset )
  return "%sfailed. [%s]%s" % ( red, B, reset )

checks={}

checkfile=open("%s/checks.txt" % Dir )
checklines=checkfile.readlines()
checkfile.close()

reader="LHEReader_cc.so"

check=[]

debug=False

for a in sys.argv[1:]:
  if a[:2]=="-r":
    reader=a[2:]
    continue
  if a[:2]=="-d":
    debug=True
    continue
  check.append(a)

for line in checklines:
  tokens=line.split(":")
  if not len(tokens)==2:
    print "[regression.py] cannot parse",line
    sys.exit(0)
  File=tokens[0]
  Pass=(len(check)==0)
  for c in check:
    if File.find(c)>-1:
      Pass=True
  if not Pass: continue
  SMSes=tokens[1].replace("\n","").split()
  if len(SMSes)==0:
    print "[regression.py] cannot parse",line
    sys.exit(0)
  if len(SMSes)>1:
    print "[regression.py] cannot yet deal with more than one event"

  o=lheDecomposer.readSequential ( "%s/%s" % ( Dir, File ), readerso=reader, debug=debug)
  if len(o)!=len(SMSes) and len(SMSes)!=1:
    print "[regression.py] lhe file contains %d events, but the check.txt file specifies only %d SMSes." % (len(o),len(SMSes))
    sys.exit(0)
  # print "Check: %s %s %d %s" % (File, SMSes[0],len(SMSes),o)
  for (ctr,foundSMS) in enumerate(o):
    compareAgainst=SMSes[0]
    if len(SMSes)>1: compareAgainst=SMSes[ctr]
    Check= ( foundSMS == compareAgainst )
    checks[File+"_"+str(ctr)]=Check
  print "Checking %20s %7s: %s" % (File,"["+foundSMS+"]",ok(foundSMS,compareAgainst))

passed=sum(checks.values())
All=len(checks)
color=red
if passed==All: color=green
print "%s%s%s/%d tests passed.%s" % ( green, passed, color, All, reset )
