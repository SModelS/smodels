#!/usr/bin/env python

"""
.. module:: lheChecks
   :synopsis: Check LHE file format.

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""

from __future__ import print_function
from smodels.tools.ioObjects import LheStatus

def main(args):   
    status = LheStatus(args.filename)
    print(status.status)
