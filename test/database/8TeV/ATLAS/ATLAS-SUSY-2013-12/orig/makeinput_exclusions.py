#!/usr/bin/env python

import os

for f in os.listdir("full_data/"):
  if not "exc" in f: continue
  read = None
  of = open(f,"w")
  for l in open("full_data/"+f):
    if "xdesc" in l:
      read = True
      continue
    if not read or not l.strip(): continue
    ol = l.split()
    os = ol[0]+"  "+ol[3]+"\n"
    of.write(os)
  of.close()
