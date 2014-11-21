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
from smodels.theory.ioObjects import SlhaStatus 

if __name__ == "__main__":
    argparser = argparse.ArgumentParser() # pylint: disable-msg=C0103
    argparser.add_argument('-f', '--filename',
                           help = 'filename of input slha')
    argparser.add_argument('-maxfl', '--flightlength',
                           help = 'maximum c*tau in m',
                           default = 3.)
    argparser.add_argument('-mD', '--decays',
                           help = 'find missing decay blocks',
                           action = 'store_false')
    argparser.add_argument('-fE', '--empty',
                           help = 'find empty decay blocks',
                           action = 'store_false')
    argparser.add_argument('-xS', '--xsec',
                           help = 'check if file contains xsection block',
                           action = 'store_false')
    argparser.add_argument('-lsp', '--lsp',
                           help = 'check if lsp is neutral and colorless',
                           action = 'store_false')
    argparser.add_argument('-ctau', '--ctau',
                           help = 'check if nlsp has prompt decay',
                           action = 'store_false')
    argparser.add_argument('-maxDisp', '--displacement',
                           help = 'give maximum displacement of secondary vertex in m',
                           default = .001)
    argparser.add_argument('-sigmacut','--sigmacut',
                           help = 'give sigmacut in fb',
                           default = .01)
    argparser.add_argument('-fD', '--displaced',
                           help = 'find displaced vertices',
                           action = 'store_false')
    argparser.add_argument('-mg', '--massgap',
                           help= 'give massgap for mass compression in GeV',
                           default = 5.)
    args = argparser.parse_args() # pylint: disable-msg=C0103
    status = SlhaStatus(args.filename, args.flightlength, args.displacement,
                        args.sigmacut, args.decays, args.displaced,
                        args.massgap, args.empty, args.xsec, args.lsp,
                        args.ctau) # pylint: disable-msg=C0103
    print(status.status)
