#!/usr/bin/env python

"""
.. module:: slhaChecks
   :synopsis: Check SLHA file for integrity.

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>

"""

from __future__ import print_function
import argparse
from smodels.tools.ioObjects import SlhaStatus
from smodels.tools.physicsUnits import fb, GeV

if __name__ == "__main__":
    argparser = argparse.ArgumentParser() # pylint: disable-msg=C0103
    argparser.add_argument('-f', '--filename',
                           help = 'filename of input slha')
    argparser.add_argument('-xS', '--xsec',
                           help = 'check if file contains xsection block',
                           action = 'store_false')
    argparser.add_argument('-lsp', '--lsp',
                           help = 'check if lsp is neutral and colorless',
                           action = 'store_false')
    argparser.add_argument('-longlived', '--longlived',
                           help = 'check for stable charged particles and visible displaced vertices',
                           action = 'store_false')
    argparser.add_argument('-maxDisp', '--displacement',
                           help = 'give maximum displacement of secondary vertex in m',
                           default = .001)
    argparser.add_argument('-sigmacut','--sigmacut',
                           help = 'give sigmacut in fb',
                           default = .01)
    argparser.add_argument('-illegal','--illegal',
                           help= 'check if all decays are kinematically allowed', 
                           action = 'store_true')
    args = argparser.parse_args() # pylint: disable-msg=C0103
    status = SlhaStatus(args.filename, maxDisplacement=args.displacement,
                        sigmacut=args.sigmacut*fb, checkLSP=args.lsp,
                        findIllegalDecays=args.illegal, checkXsec=args.xsec,
                        findLonglived=args.longlived) # pylint: disable-msg=C0103
    print(status.status)
