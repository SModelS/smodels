#!/usr/bin/python

import time
import sys
sys.path.insert(0,"../")
from runSModelS import main
from smodels.installation import installDirectory

filename = "%s/inputFiles/slha/gluino_squarks.slha" % \
            (installDirectory() )
parameterFile="%s/test/timingParameters.ini" %installDirectory(), 
out = "%s/test/unitTestOutput.txt" % installDirectory()

suppressStdout = False

if suppressStdout:
    a=sys.stdout
    sys.stdout = open ( "stdout.log", "w" )
t0=time.time()
main( filename, parameterFile, "/tmp", None, 0, True )
t1=time.time()

if suppressStdout:
    sys.stdout = a
print "%.1f secs." % ( t1-t0 )

print "on wnouc, old version: 71 secs."                                               
print "on wnouc, without persistent triangulation: 45 secs."                                               
