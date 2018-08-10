#!/usr/bin/env python

import os

for f in os.listdir("full"):
  read = None
  of = open(f,"w")
  for l in open("full/"+f):
    if "xdesc" in l:
      read = True
      continue
    if not read or not l.strip(): continue
    ol = l.split()
    os = ol[0]+"  "+ol[3]+"  "+ol[6]+"\n"
    of.write(os)
  of.close()
