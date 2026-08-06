#!/usr/bin/env python3

"""
.. module:: lheChecks
   :synopsis: Check LHE file format.

.. moduleauthor:: Ursula Laa <ursula.laa@lpsc.in2p3.fr>

"""

from smodels.tools.ioObjects import LheStatus

def main(args):   
    status = LheStatus(args.filename)
    print(status.status)
