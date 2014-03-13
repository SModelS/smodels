#!/usr/bin/env python

""" check the output of the cross section computer """

import set_path
from tools import VariousHelpers
from theory import XSecComputer

w=XSecComputer.compute( 1000, "../slha/andrePT1.slha" )
print "w=",w.weights()
