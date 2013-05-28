#!/usr/bin/env python

""" Test 1: check the output of the cross section computer """

import sys, set_path
from Theory import XSecComputer

w=XSecComputer.compute( 1000, "../slha/andrePT1.slha" )
print "w=",w
