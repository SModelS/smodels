#!/usr/bin/env python

"""
.. module:: smodels
   :synopsis: Command line program for SModelS tasks.

.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>

"""

import argparse
from smodels.tools import xsecComputer
from smodels.tools import slhaChecks


def main():
    parser = argparse.ArgumentParser(description="SModelS command line tool.")
    subparsers = parser.add_subparsers(dest='subparser_name')
    
    xseccomputer = subparsers.add_parser('xseccomputer', description="Compute MSSM cross sections for a SLHA file.")        
    xseccomputer.add_argument('-s', '--sqrts', nargs='+', action='append', help="sqrt(s) TeV. Can supply more than one value.", type=int, default=[])
    xseccomputer.add_argument('-e', '--nevents', type=int, default=10000, help="number of events to be simulated.")
    xseccomputer.add_argument('-p', '--tofile', action='store_true', help="write cross sections to file")
    xseccomputer.add_argument('-k', '--keep', action='store_true', help="do not unlink temporary directory")
    xseccomputer.add_argument('-n', '--NLO', action='store_true', help="compute at the NLO level (default is LO)")
    xseccomputer.add_argument('-N', '--NLL', help="compute at the NLL level (takes precedence over NLL, default is LO)", action='store_true')
    xseccomputer.add_argument('-O', '--LOfromSLHA',help="use LO cross-sections from file to compute the NLO or NLL cross-sections", action='store_true')
    xseccomputer.add_argument('-f', '--filename', help="SLHA file to compute cross sections for", required=True)
      
    slhachecker = subparsers.add_parser('slhachecker', description="Perform several checks on a SLHA file.")    
    slhachecker.add_argument('-xS', '--xsec', help = 'check if file contains xsection block', action = 'store_false')
    slhachecker.add_argument('-lsp', '--lsp', help = 'check if lsp is neutral and colorless', action = 'store_false')
    slhachecker.add_argument('-longlived', '--longlived', help = 'check for stable charged particles and visible displaced vertices', action = 'store_false')
    slhachecker.add_argument('-m', '--displacement', help = 'give maximum displacement of secondary vertex in m', default = .001, type = float)
    slhachecker.add_argument('-s','--sigmacut', help = 'give sigmacut in fb', default = .01, type = float)
    slhachecker.add_argument('-illegal','--illegal', help= 'check if all decays are kinematically allowed', action = 'store_true')
    slhachecker.add_argument('-f', '--filename', help = 'name of input SLHA file', required=True)
    
    
    args = parser.parse_args()
    
    if args.subparser_name == 'xseccomputer':
        xsecComputer.main(args)
    if args.subparser_name == 'slhachecker':    
        slhaChecks.main(args)

if __name__ == '__main__':
    main()
    
