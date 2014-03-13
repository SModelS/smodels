#!/usr/bin/env python

""" check the unique naming facility for slha files """

import set_path
from tools.VariousHelpers import logging
from theory import SLHATools
slhafile="../slha/andrePT1.slha"
print SLHATools.uniqueName ( slhafile )
