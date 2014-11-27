#!/usr/bin/env python

"""
.. module:: slhaChecks
   :synopsis: Check SLHA file for integrity.

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>

"""

from __future__ import print_function
from smodels.tools.ioObjects import SlhaStatus
from smodels.tools.physicsUnits import fb

def main(args):   
    status = SlhaStatus(args.filename, maxDisplacement=args.displacement,
                        sigmacut=args.sigmacut*fb, checkLSP=args.lsp,
                        findIllegalDecays=args.illegal, checkXsec=args.xsec,
                        findLonglived=args.longlived, findMissingDecayBlocks=args.decayBlocks) # pylint: disable-msg=C0103
    print(status.status)
