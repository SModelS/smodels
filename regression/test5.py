#!/usr/bin/env python

""" check the unique naming facility for slha files """

import set_path
from Tools.VariousHelpers import logging
from Theory import SLHATools
slhafile="../slha/andrePT1.slha"
print SLHATools.uniqueName ( slhafile )
