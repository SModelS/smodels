#!/usr/bin/env python

"""
.. module:: smodels
   :synopsis: Command line program for SModelS tasks.

.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>

"""

import argparse
from smodels.tools import xsecComputer


def main():
    parser = argparse.ArgumentParser(description="SModelS command line tool.")
    subparsers = parser.add_subparsers(dest='subparser_name')
    
    xseccomputer = subparsers.add_parser('xseccomputer', description="Compute the cross section of an SLHA file.")
    xseccomputer.add_argument('file', type=str, nargs=1, help="SLHA file to compute cross section for")
    xseccomputer.add_argument('-s', '--sqrts', nargs='+', action='append', help="sqrt(s) TeV. Can supply more than one value.", type=int, default=[])
    xseccomputer.add_argument('-e', '--nevents', type=int, default=10000, help="number of events to be simulated.")
    xseccomputer.add_argument('-p', '--tofile', action='store_true', help="write cross sections to file")
    # xseccomputer.add_argument('-S', '--slha', action='store_true', help="input file is an SLHA file")
    # xseccomputer.add_argument('-L', '--lhe', action='store_true', help="input file is an LHE file")
    xseccomputer.add_argument('-k', '--keep', action='store_true', help="do not unlink temporary directory")
    xseccomputer.add_argument('-n', '--NLO', action='store_true', help="compute at the NLO level (default is LO)")
    xseccomputer.add_argument('-N', '--NLL', help="compute at the NLL level (takes precedence over NLL, default is LO)", action='store_true')
    xseccomputer.add_argument('-O', '--LOfromSLHA',help="use LO cross-sections from file to compute the NLO or NLL cross-sections", action='store_true')
    
    
    args = parser.parse_args()
    
    if args.subparser_name == 'xseccomputer':
        xsecComputer.main(args)


if __name__ == '__main__':
    main()
    
