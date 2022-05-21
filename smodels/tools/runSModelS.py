#!/usr/bin/env python3

"""
.. module:: runSModelS
   :synopsis: Main code for running SModelS.

"""

from __future__ import print_function
import os
from smodels.installation import installDirectory, version
from smodels.tools import modelTester
from smodels.tools import crashReport
from smodels.tools import smodelsLogging
from smodels.tools import runtime
from smodels import particlesLoader
from imp import reload

def main():
    import argparse
    """ Set default input and output files """
    parameterFile = "%s/smodels/etc/parameters_default.ini" % installDirectory()
    outputDir = "./results/"

    """ Get the name of input SLHA file and parameter file """
    ap = argparse.ArgumentParser( description=
            "Run SModelS over SLHA/LHE input files." )
    ap.add_argument('-f', '--filename',
            help='name of SLHA or LHE input file or a directory path (required argument). '
            'If a directory is given, loop over all files in the directory', required=True)
    ap.add_argument('-p', '--parameterFile',
            help='name of parameter file, where most options are defined (optional argument). If not set, use '
            'all parameters from smodels/etc/parameters_default.ini',
            default=parameterFile)
    ap.add_argument('-o', '--outputDir',
            help='name of output directory (optional argument). The default folder is: ' +
            outputDir, default=outputDir)
    ap.add_argument('-d', '--development', help='if set, SModelS will run in development mode and exit if any errors are found.',
            action='store_true')
    ap.add_argument('-t', '--force_txt', help='force loading the text database',
            action='store_true')
    ap.add_argument('-C', '--colors', help='colored output',
            action='store_true')
    ap.add_argument('-V', '--version', action='version', version = version() )
    ap.add_argument('-c', '--run-crashreport',
            help='parse crash report file and use its contents for a SModelS run. '
                 "Supply the crash file simply via '--filename myfile.crash'",
            action='store_true')
    ap.add_argument('-v','--verbose', help='sets the verbosity level (debug, info, warning, error). Default value is info.',
            default = "info", type = str )
    ap.add_argument('-T', '--timeout',
            help='define a limit on the running time (in secs).'
                 'If not set, run without a time limit. If a directory is given as input, '
                 'the timeout will be  applied for each individual file.',
            default = 0, type = int)


    args = ap.parse_args()
    if args.colors:
        from smodels.tools.colors import colors
        colors.on = True


    db=None
    if args.force_txt: db=True
    smodelsLogging.setLogLevel ( args.verbose )

    if args.run_crashreport:
        args.filename, args.parameterFile = crashReport.readCrashReportFile(
                args.filename)
        run(args.filename, args.parameterFile, args.outputDir,
               db, args.timeout, development=True )

    else:
        run(args.filename, args.parameterFile, args.outputDir,
              db, args.timeout, args.development)

def run( inFile, parameterFile, outputDir, db, timeout, development ):
    """
    Provides a command line interface to basic SModelS functionalities.

    :param inFile: input file name (either a SLHA or LHE file)
                   or directory name (path to directory containing input files)
    :param parameterFile: File containing the input parameters (default =
                          smodels/etc/parameters_default.ini)
    :param outputDir: Output directory to write a summary of results to
    :param db: supply a smodels.experiment.databaseObj.Database object, so
            the database doesn't have to be loaded anymore. Will
            render a few parameters in the parameter file irrelevant.
            If None, load the database as described in parameterFile,
            If True, force loading the text database.
    :param timeout: set a timeout for one model point (0 means no timeout)
    :param development: turn on development mode (e.g. no crash report)

    """

    """ Read and check parameter file, exit parameterFile does not exist """
    parser = modelTester.getParameters(parameterFile)

    """ Determine particles module from ini file, if necessary """
    if parser.has_option("particles","model"):
        runtime.modelFile = parser.get( "particles", "model" )
        reload(particlesLoader)

    """ Check database location and load database, exit if not found """
    database, databaseVersion = modelTester.loadDatabase(parser, db)

    """ Get list of input files to be tested """
    fileList, inDir = modelTester.getAllInputFiles(inFile)

    """ Create output directory if missing """
    if not os.path.isdir(outputDir): os.mkdir(outputDir)

    """ Load analysis database, print list """
    listOfExpRes = modelTester.loadDatabaseResults(parser, database)

    """ Test all input points """
    modelTester.testPoints(fileList, inDir, outputDir, parser, databaseVersion,
                 listOfExpRes, timeout, development, parameterFile )

if __name__ == "__main__":
    main()
