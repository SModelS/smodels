#!/usr/bin/env python3

"""
.. module:: slhaChecks
   :synopsis: Check SLHA file for integrity.

.. moduleauthor:: Ursula Laa <ursula.laa@lpsc.in2p3.fr>
.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>

"""

from __future__ import print_function
from smodels.tools.ioObjects import SlhaStatus

def main(args):   
    status = SlhaStatus( args.filename, 
                 findIllegalDecays=args.illegal, checkXsec=args.xsec, 
                 findMissingDecayBlocks=args.decayBlocks) # pylint: disable-msg=C0103

    print(status.status)
