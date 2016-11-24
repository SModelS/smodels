#!/usr/bin/env python

"""
.. module:: lheChecks
   :synopsis: Check LHE file format.

.. moduleauthor:: Ursula Laa <ursula.laa@lpsc.in2p3.fr>

"""

from __future__ import print_function
from smodels.tools.ioObjects import LheStatus

def main(args):   
    status = LheStatus(args.filename)
    print(status.status)
