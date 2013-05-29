#!/usr/bin/env python

""" check the output of the cross section computer """

import set_path
from Theory import XSecComputer

w=XSecComputer.compute( 1000, "../slha/andrePT1.slha" )
print "w=",w.weights()
