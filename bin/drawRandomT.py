#!/usr/bin/python

import sys, os, random
sys.path.append("/usr/local/smodels/")
from Theory import LHEReader, TopologyBuilder
from Tools import FeynmanGraphs

Dir="/usr/local/smodels/lhe/"

Files=os.listdir( Dir )

File=""
while File[-4:]!=".lhe":
  File=random.choice ( Files )

# File="T1lnu_1.lhe"

filename=Dir+File
T=File.replace(".lhe","")
while T.find("_")!=-1:
  T=T[:T.find("_")]

print 
print "Today's Random Topology is ``%s'':" % T
print

reader = LHEReader.LHEReader( filename )
Event = reader.next()
SMSTop = TopologyBuilder.fromEvent(Event, {} )
FeynmanGraphs.asciidraw ( SMSTop[0].leadingElement() )
