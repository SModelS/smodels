#!/usr/bin/env python

import sys, set_path
from Theory import XSecComputer

w=XSecComputer.compute( 1000, "../slha/andrePT1.slha" )
print "w=",w
