   smodelsTools.py xsecresummino [-h] [-s SQRTS [SQRTS ...]] [-j RESUMMINOFILE]                                    [-v VERBOSITY] [-c NCPUS] [-p] [-P]                                     [-n]                                    [-N]   -f FILENAME

options:
  -h, --help            show this help message and exit
  -s SQRTS, --sqrts SQRTS
                        sqrt(s) TeV. Can supply more than one value (as a                        space separated list). Default is both 8 and 13.
  -j RESUMMINOFILE, --json RESUMMINOFILE
                        File which contain all the relevant informations for the cross section calculation with resummino (daughter particles, mode, etc.), default is smodels/etc/resummino.json
  -v VERBOSITY, --verbosity VERBOSITY
                        Verbosity (debug, info, warning, error)
  -c NCPUS, --ncpus NCPUS
                        number of cores to be used simultaneously. -1 means                        'all'.
  -p, --tofile          write cross sections to file (only highest order)
  -P, --alltofile       write all cross sections to file, including lower
                        orders
  -n, --NLO             compute at the NLO level (default is LO)
  -N, --NLL             compute at the NLO+NLL level (takes precedence over
                        NLO, default is LO)
  -f FILENAME, --filename FILENAME
                        SLHA file to compute cross sections for. If a                        directory is given, compute cross sections for all                        files in directory.