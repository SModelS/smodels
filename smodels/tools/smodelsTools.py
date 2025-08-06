#!/usr/bin/env python3

"""
.. module:: smodelsTools
   :synopsis: Command line program for SModelS tools.

.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function

def main():
    import argparse
    parser = argparse.ArgumentParser(description="SModelS-tools command line tool.")

    parser.add_argument('-v','--verbose', help='verbosity level. '
                        'accepted values are: debug, info, warning, error. [info]',
                                    default = "info", type = str )

    subparsers = parser.add_subparsers(dest='subparser_name')

    subparsers.add_parser('installation', description="Print installation setup and exit.")
    subparsers.add_parser('fixpermissions', description="Fix file permissions for xseccomputer.")
    xseccomputer = subparsers.add_parser('xseccomputer', description="Compute MSSM cross sections for a SLHA file.")
    xseccomputer.add_argument('-f', '--filename', required=True,
            help="SLHA file to compute cross sections for. "
            "If a directory is given, cross sections for all files in the directory are computed." )
    xseccomputer.add_argument('-s', '--sqrts', nargs='+', action='append',
        help="LHC center-of-mass energy in TeV for computing the "
        "cross sections. Can be more than one value, e.g., -s 8 13 for both "
        "8 TeV and 13 TeV cross sections. [13]",
        type=float, default=[])
    xseccomputer.add_argument('-6', '--pythia6', action='store_true',
        help="use pythia6 for LO cross sections")
    xseccomputer.add_argument('-8', '--pythia8', action='store_true',
        help="use pythia8 for LO cross sections (default)")
    xseccomputer.add_argument('-e', '--nevents', type=int, default=10000,
        help="number of events to be simulated [10000].")
    xseccomputer.add_argument('-v', '--verbosity', type=str, default="info",
        help="Verbosity (debug, info, warning, error) [info]")
    xseccomputer.add_argument('-c', '--ncpus', type=int, default=1,
        help="number of CPU cores to be used simultaneously. −1 "
        "means ‘all’. Used only when cross sections are computed for multiple "
        "SLHA files. [1]")
    xseccomputer.add_argument('-p', '--tofile', action='store_true',
        help="write cross sections to file (only highest order)")
    xseccomputer.add_argument('-P', '--alltofile', action='store_true',
        help="write all cross sections to file, including lower orders")
    xseccomputer.add_argument('-n', '--NLO', action='store_true',
        help="compute at the NLO level (default is LO)")
    xseccomputer.add_argument('-N', '--NLL', help="compute at the NLO+NLL level (takes precedence over NLO, default is LO)", action='store_true')
    xseccomputer.add_argument('-O', '--LOfromSLHA', help="use LO cross sections from file to compute the NLO or NLL cross sections", action='store_true')
    xseccomputer.add_argument('-S', '--ssmultipliers', type=str, default=None,
        help="Signal strength multipliers, provided as dictionary of pids")
    xseccomputer.add_argument('-q', '--query', action='store_true',
        help="only query if there are cross sections in the file")
    xseccomputer.add_argument('-C', '--colors', action='store_true',
        help="colored terminal output" )
    xseccomputer.add_argument('-k', '--keep', action='store_true',
        help="do not unlink temporary directory")
    xseccomputer.add_argument( '--noautocompile', action='store_true',
        help="turn off automatic compilation" )

    xsecresummino = subparsers.add_parser('xsecresummino', description="Compute gaugino and slepton cross sections via resummino for a SLHA file.")
    xsecresummino.add_argument('-f', '--filename', required=True,
            help="SLHA file to compute cross sections for. "
            "If a directory is given, cross sections for all files in the directory are computed." )
    xsecresummino.add_argument('-s', '--sqrts', nargs='+', action='append',
        help="LHC center-of-mass energy in TeV for computing the "
        "cross sections. Can be more than one value, e.g., -s 8 13 for both "
        "8 TeV and 13 TeV cross sections. [13]",
        type=float, default=[])
    xsecresummino.add_argument('-part', '--particles', nargs='+', action='append',
        help="list of daughter particles (given as PDG "
        "codes) to compute cross sections for. All valid combinations from "
        "the list will be considered. If no list of particles is given, the channels "
        "info from the resummino.py configuration file is used instead.",
        type=int, default=[])
    xsecresummino.add_argument('-v', '--verbosity', type=str, default="info",
        help="verbosity (debug, info, warning, error). [info]")
    xsecresummino.add_argument('-c', '--ncpus', type=int, default=1,
        help="number of CPU cores to be used simultaneously. -1 means 'all'. Used only when cross sections are computed for multiple SLHA files. [1]")
    xsecresummino.add_argument('-C', '--conf', type=str, default='default',
        help="path to resummino.py configuration file. [smodels/etc/resummino.py]")
    xsecresummino.add_argument('-x', '--xseclimit', type=float, default=None,
        help="cross section limit in pb. If the LO cross section is "
        "below this value, no higher orders will be calculated. The default "
        "is 0.00001, set in the smodels/etc/resummino.py file." )
    xsecresummino.add_argument('-p', '--tofile', action='store_true',
        help="write cross sections to file (only highest order)")
    xsecresummino.add_argument('-P', '--alltofile', action='store_true',
        help="write all cross sections to file, including lower orders")
    xsecresummino.add_argument('-n', '--NLO', action='store_true',
        help="compute at the NLO level (default is LO)")
    xsecresummino.add_argument('-N', '--NLL', help="compute at the NLO+NLL level (takes precedence over NLO, default is LO)", action='store_true')
    # xsecresummino.add_argument('-S', '--ssmultipliers', type=str, default=None,
    #     help="signal strength multipliers, provided as dictionary of pids")
    xsecresummino.add_argument('-k', '--keep', action='store_true',
        help="do not unlink temporary directory")
    xsecresummino.add_argument( '--noautocompile', action='store_true',
        help="turn off automatic compilation" )

    slhachecker = subparsers.add_parser('slhachecker', description="Perform several checks on a SLHA file.")
    slhachecker.add_argument('-xS', '--xsec', help='turn off the check for xsection blocks', action='store_false')
    slhachecker.add_argument('-illegal', '--illegal', help='turn on check for kinematically forbidden decays', action='store_true')
    slhachecker.add_argument('-dB', '--decayBlocks', help='turn off the check for missing decay blocks', action='store_false')
    slhachecker.add_argument('-f', '--filename', help='name of input SLHA file', required=True)

    lhechecker = subparsers.add_parser('lhechecker', description="Check if the input file has LHE format.")
    lhechecker.add_argument('-f', '--filename', help='name of input LHE file', required=True)

    dbBrowser = subparsers.add_parser('database-browser', description="Interface for browsing the Database.")
    dbBrowser.add_argument('-p', '--path_to_database', help='path to SModelS database', required=True)
    dbBrowser.add_argument('-t', '--text', help='load text database, dont even search for binary database file', action='store_true')

    iPlots = subparsers.add_parser('interactive-plots', description="Produces a set of interactive plots for visualizing results from a scan.")
    iPlots.add_argument('-p', '--parameters', help='path to the parameters file [./smodels/etc/iplots_parameters.py]', default = './smodels/etc/iplots_parameters.py')

    iPlots.add_argument('-m', '--modelFile', help='path to the model.py file', default = None)

    iPlots.add_argument('-f', '--smodelsFolder', help='path to the smodels folder or tarball (.tar.gz) with the SModelS python output files.',
                            required=True)
    iPlots.add_argument('-s', '--slhaFolder', help='path to the SLHA folder or tarball (.tar.gz) with the SLHA input files.',
                            required=True)
    iPlots.add_argument('-o', '--outputFolder',
                  help='path to the output folder, where the plots will be stored. [./iplots]',
                  default = "./iplots" )

    iPlots.add_argument('-N', '--npoints', type=int, default=-1,
        help="How many (randomly selected) points will be included in the plot. If -1 all points will be read and included (default = -1).")

    iPlots.add_argument('-v', '--verbosity', type=str, default="info",
        help="Verbosity (debug, info, warning, error) [info]")

    toolbox = subparsers.add_parser( 'toolbox', description=
								                     "Facility to control external dependencies")
    toolbox.add_argument('-c', '--colors', help='turn on terminal colors',
                           action='store_true')
    toolbox.add_argument('-l', '--long', help='long output lines',
                           action='store_true')
    toolbox.add_argument('-m', '--make', help='compile packages if needed',
                           action='store_true')
    environ_info = subparsers.add_parser( 'environ-info', description=
								                     "Simple tool to print out environment info")
    environ_info.add_argument('-c', '--colors', help='turn on terminal colors',
                           action='store_true')
    args = parser.parse_args()

    from smodels.base import smodelsLogging
    smodelsLogging.setLogLevel(args.verbose)

    from smodels import installation
    if args.subparser_name == 'fixpermissions':
        installation.fixpermissions()

    if args.subparser_name == "environ-info":
        from smodels.base.runtime import printEnvironmentInfo
        printEnvironmentInfo( vars(args) )

    if args.subparser_name == 'installation':
        import sys, os
        print(installation.banner())
        print("SModelS version:", installation.version())
        print("Installation directory:",installation.installDirectory())
        path = os.path.abspath(os.path.realpath(__file__))
        print("This executable:",path )
        sys.exit()
    if args.subparser_name == 'toolbox':
        from smodels.tools import toolBox
        toolBox.main ( args )
    if args.subparser_name == 'xseccomputer':
        from smodels.tools import xsecComputer
        xsecComputer.main(args)
    if args.subparser_name == 'xsecresummino':
        from smodels.tools import xsecResummino
        xsecResummino.main(args)
    if args.subparser_name == 'slhachecker':
        from smodels.tools import slhaChecks
        slhaChecks.main(args)
    if args.subparser_name == 'lhechecker':
        from smodels.tools import lheChecks
        lheChecks.main(args)
    if args.subparser_name == 'database-browser':
        from smodels.tools import databaseBrowser
        databaseBrowser.main(args)
    if args.subparser_name == 'interactive-plots':
        from smodels.tools import interactivePlots
        interactivePlots.main(args)

if __name__ == '__main__':
    main()


