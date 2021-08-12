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
                        'accepted values are: debug, info, warning, error.',
                                    default = "info", type = str )

    subparsers = parser.add_subparsers(dest='subparser_name')

    subparsers.add_parser('installation', description="Print installation setup and exit.")
    subparsers.add_parser('fixpermissions', description="Fix file permissions for xseccomputer.")
    xseccomputer = subparsers.add_parser('xseccomputer', description="Compute MSSM cross sections for a SLHA file.")
    xseccomputer.add_argument('-s', '--sqrts', nargs='+', action='append',
        help="sqrt(s) TeV. Can supply more than one value (as a space separated list). Default is both 8 and 13.",
        type=int, default=[])
    xseccomputer.add_argument('-e', '--nevents', type=int, default=10000,
        help="number of events to be simulated [10000].")
    xseccomputer.add_argument('-v', '--verbosity', type=str, default="info",
        help="Verbosity (debug, info, warning, error)")
    xseccomputer.add_argument('-c', '--ncpus', type=int, default=-1,
        help="number of cores to be used simultaneously. -1 means 'all'. ")
    xseccomputer.add_argument('-p', '--tofile', action='store_true',
        help="write cross sections to file (only highest order)")
    xseccomputer.add_argument('-P', '--alltofile', action='store_true',
        help="write all cross sections to file, including lower orders")
    xseccomputer.add_argument('-q', '--query', action='store_true',
        help="only query if there are cross sections in the file")
    xseccomputer.add_argument('-C', '--colors', action='store_true',
        help="colored terminal output" )
    xseccomputer.add_argument( '--noautocompile', action='store_true',
        help="turn off automatic compilation" )
    xseccomputer.add_argument('-k', '--keep', action='store_true',
        help="do not unlink temporary directory")
    xseccomputer.add_argument('-6', '--pythia6', action='store_true',
        help="use pythia6 for LO cross sections")
    xseccomputer.add_argument('-8', '--pythia8', action='store_true',
        help="use pythia8 for LO cross sections (default)")
    xseccomputer.add_argument('-n', '--NLO', action='store_true',
        help="compute at the NLO level (default is LO)")
    xseccomputer.add_argument('-N', '--NLL', help="compute at the NLO+NLL level (takes precedence over NLO, default is LO)", action='store_true')
    xseccomputer.add_argument('-O', '--LOfromSLHA', help="use LO cross sections from file to compute the NLO or NLL cross sections", action='store_true')
    xseccomputer.add_argument('-S', '--ssmultipliers', type=str, default=None,
        help="Signal strength multipliers, provided as dictionary of pids")
    xseccomputer.add_argument('-f', '--filename', required=True,
            help="SLHA file to compute cross sections for. "
            "If a directory is given, compute cross sections for all files in directory." )

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
    iPlots.add_argument('-p', '--parameters', help='path to the parameters file [./iplots_parameters.py]', default = './iplots_parameters.py')

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
        help="Verbosity (debug, info, warning, error)")

    proxydb = subparsers.add_parser( 'proxydb', description=
								                     "create proxy databases for network use")
    proxydb.add_argument ( '-i', '--inputfile',
            help='input database file [db31.pcl]',
            type=str, default="db31.pcl"  )
    proxydb.add_argument ( '-s', '--servername',
            help='the server name [None]',
            type=str, default=None  )
    proxydb.add_argument ( '-v', '--verbose',
            help='verbosity, for server and client [info]',
            type=str, default="info" )
    proxydb.add_argument ( '-R', '--rundir',
            help='the rundir [./]',
            type=str, default="./"  )
    proxydb.add_argument ( '-p', '--serverport',
            help='the server port [31770]',
            type=int, default=None  )
    proxydb.add_argument ( '-l', '--symlink',
            help='set a symlink from outputfile to default.pcl',
            action='store_true'  )
    proxydb.add_argument ( '-o', '--outputfile',
            help='input database file [proxy31.pcl]',
            type=str, default="proxy31.pcl"  )
    proxydb.add_argument('-r', '--run', help='start the server upon creation',
                           action='store_true')

    toolbox = subparsers.add_parser( 'toolbox', description=
								                     "Facility to control external dependencies")
    toolbox.add_argument('-c', '--colors', help='turn on terminal colors',
                           action='store_true')
    toolbox.add_argument('-l', '--long', help='long output lines',
                           action='store_true')
    toolbox.add_argument('-m', '--make', help='compile packages if needed',
                           action='store_true')
    args = parser.parse_args()

    from smodels.tools import smodelsLogging
    smodelsLogging.setLogLevel(args.verbose)

    from smodels import installation
    if args.subparser_name == 'fixpermissions':
        installation.fixpermissions()

    if args.subparser_name == 'installation':
        import sys, os
        print(installation.banner())
        print("SModelS version:", installation.version())
        print("Installation directory:",installation.installDirectory())
        path = os.path.abspath(os.path.realpath(__file__))
        print("This binary:",path )
        sys.exit()
    if args.subparser_name == 'toolbox':
        from smodels.tools import toolBox
        toolBox.main ( args )
    if args.subparser_name == 'proxydb':
        from smodels.tools import proxyDBCreator
        proxyDBCreator.main ( args )
    if args.subparser_name == 'xseccomputer':
        from smodels.tools import xsecComputer
        xsecComputer.main(args)
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

