#!/usr/bin/env python

import sys, os, commands, argparse, logging
from ConfigParser import SafeConfigParser
import setPath
from smodels.tools.physicsUnits import rmvunit, addunit
from smodels.tools import slhaChecks, ioObjects, xsecComputer
from smodels.experiment import smsHelpers, smsAnalysisFactory, smsResults
from smodels.theory import slhaDecomposer
from smodels.theory.theoryPrediction import theoryPredictionFor,  _getElementsFrom
from smodels.theory import crossSection

log = logging.getLogger(__name__)

#get the name of input slha file (and parameter file)
argparser = argparse.ArgumentParser()
argparser.add_argument('-f', '--filename', help = 'name of SLHA or LHE input file, necessary input', required = True)
argparser.add_argument('-p', '--parameterfile', help = 'name of parameter file, optional argument, default is: ./parameters.ini', default = './parameters.ini')
argparser.add_argument('-o', '--outputfile', help = 'name of output file, optional argument, default is: ./summary.txt', default = 'summary.txt')
args = argparser.parse_args() # pylint: disable-msg=C0103

slhafile = args.filename #get input filename

#use ConfigParser to read input parameters
parser = SafeConfigParser()
#first, read defaults
parser.read('etc/parameters_default.ini')
#then, set parameters given in input file
parser.read(args.parameterfile)

if os.path.exists(args.outputfile): #remove old output file
    log.warning("Removing old output file in "+args.outputfile)
    os.remove(args.outputfile)

#lists to store results
bestresult = []
outputarray = []

#get status, warnings from slhaChecks

#slhatstat == 0, the file is not checked
#slhastat == 1, the check is ok
#slhastat == -1, physical problem in scenario, e.g. charged LSP,
#slhastat == -2  formal problems in the file, e.g. missing decay blocks

slhaStatus = slhaChecks.SlhaStatus(slhafile, sigmacut = parser.getfloat("parameters","sigmacut"), maxDisplacement = .001, checkXsec = not parser.getboolean("options","addMissingXsecs"), massgap = parser.getfloat("parameters","minmassgap"), maxcond = parser.getfloat("parameters","maxcond"))
slhastat, warnings = slhaStatus.status

if slhastat == -1 or slhastat == -3:
    status = -2
    outputStatus = ioObjects.OutputStatus(status, slhastat, warnings)
    outputStatus.printout("file",args.outputfile)
    slhaStatus.printout("file",args.outputfile)
    sys.exit()

# << there should be some automated mechanism to check this thing?>>
#writeXsecs = True #set this true by default, switch to false only in case of lhe decomposition

#set cross section computation bools according to input parameters, slha
'''if not ioPar.doSLHAdec:
    commands.getoutput("./tools/xsecComputer.py %s -S -k -e %i -s %d " %(slhafile,ioPar.nevts, ioPar.sqrts))
else:
    if ioPar.addMissingXsecs:
        commands.getoutput("./tools/xsecComputer.py %s -S -e %i -s %d -N -f %s " %(slhafile,ioPar.nevts, ioPar.sqrts, slhafile))
    elif slhaStatus.xsec:
        log.warning("Input file does not contain cross sections, set computeXsecs = True")
        warnings = warnings + "#Cross sections computed by SModelS\n"
        commands.getoutput("./tools/xsecComputer.py %s -S -e %i -s %d -N -f %s " %(slhafile,ioPar.nevts, ioPar.sqrts, slhafile))
    else: computeXsecs = None

if ioPar.addnlo:
    print "Now computing NLO"
    xsecs_nlo = xsecComputer.computeXSec(sqrts, 1, ioPar.nevts, slhafile, loFromSlha=True)
    xsecComputer.addXSecToFile(xsecs_nlo, slhafile) #FIXME should there be a comment << YES >>
if ioPar.addnll:
    xsecs_nll = xsecComputer.computeXSec(sqrts, 2, ioPar.nevts, slhafile, loFromSlha=True)
    xsecComputer.addXSecToFile(xsecs_nll, slhafile) #FIXME also: comment? << YES >>

sys.exit(10)'''

#set cross section computation bools according to input parameters, slha
if not parser.getboolean("options","doSLHAdec"):
    computeXsecs = True
    writeXsecs = None
else:
    if parser.getboolean("options","addMissingXsecs"):
        computeXsecs = True
    elif slhaStatus.xsec[0]==-1:
        log.warning("Input file does not contain cross sections, set computeXsecs = True")
        warnings = warnings + "#Cross sections computed by SModelS\n"
        computeXsecs = True
        addnll = True
    else: computeXsecs = None

# sqrts from parameter input file
sqrts = addunit(parser.getfloat("options","sqrts"),"TeV")


if computeXsecs:
    #first compute at LO
    xsecs = xsecComputer.computeXSec(sqrts, 0, parser.getint("options","nevts"), slhafile)
    comment = "Nevts: " + parser.get("options","nevts")
    xsecComputer.addXSecToFile(xsecs, slhafile, comment)

if parser.getboolean("options","addnlo"):
    print "Now computing NLO"
    xsecs_nlo = xsecComputer.computeXSec(sqrts, 1, parser.getint("options","nevts"), slhafile, loFromSlha=True)
    xsecComputer.addXSecToFile(xsecs_nlo, slhafile) #FIXME should there be a comment << YES >>

if parser.getboolean("options","addnll"):
    xsecs_nll = xsecComputer.computeXSec(sqrts, 2, parser.getint("options","nevts"), slhafile, loFromSlha=True)
    xsecComputer.addXSecToFile(xsecs_nll, slhafile) #FIXME also: comment? << YES >>

#decomposition
#sigmacut = minimum value of cross-section for an element to be considered eligible for decomposition. Too small sigmacut leads to too large deocmposition time. 
sigmacut = addunit(parser.getfloat("parameters","sigmacut"),"fb")

try:
    # Decompose input SLHA file, store the output elements in smstoplist
    smstoplist = slhaDecomposer.decompose(slhafile, sigmacut, doCompress=parser.getboolean("options","doCompress"), doInvisible=parser.getboolean("options","doInvisible"), minmassgap=addunit(parser.getfloat("parameters","minmassgap"),"GeV"))
except:
    status = -1 #<< is this a decomposition status? >>
    outputStatus = ioObjects.OutputStatus(status, slhastat, warnings)
    outputStatus.printout("file", args.outputfile)
    slhaStatus.printout("file", args.outputfile)
    sys.exit()

# Safety switch in case decomposition fails
if not smstoplist:
    status = -3 #<< is this a decomposition status? >>
    outputStatus = ioObjects.OutputStatus(status, slhastat, warnings)
    outputStatus.printout("file", args.outputfile)
    slhaStatus.printout("file", args.outputfile)
    sys.exit()

if parser.getboolean("stdout","printGtop"):
    smstoplist.printout()

# This is my porposed format for element tabel
if parser.getboolean("stdout","printThEl"):
    for (i,topo) in enumerate(smstoplist):
        print '\n'
        print "A new global topoloy starts here" 
        print "====================================================================="
        for j, el in enumerate(topo.elementList):
            print "\t ........................................................................."
            el.printout()
    print "====================================================================="
    print "====================================================================="
    print "The list ends here" 
    print "====================================================================="
    print "====================================================================="

# Set database address
smsHelpers.base = parser.get("path","database")

analyses = parser.get("data","analyses")
#in case that a list of analyses is given, retrieve list
if "," in analyses: analyses = analyses.split(",")

topologies = parser.get("data","topologies")
if "," in topologies: topologies = topologies.split(",")

# Load analyses
listofanalyses = smsAnalysisFactory.load(analyses, topologies)

#<< What is this?>>
results = ioObjects.ResultList(bestresultonly = not parser.getboolean("file","expandedSummary"), describeTopo = parser.getboolean("file","describeTopo"))

constrainedElements = []

if parser.getboolean("stdout","printAnaEl"):
    for analysis in listofanalyses:
        elements = _getElementsFrom(smstoplist, analysis)
        if len(elements) == 0: continue
        # This is my porposed format for analyses elements table
        print "========================================================"
        print "Analysis Name:", analysis.label.split(":")[0]
        print "Analysis Topology:", analysis.label.split(":")[1]
        print "Analysis Sqrts:", analysis.sqrts
        print "========================================================"
        ref_el = None
        for el in elements:
            el.printout()
            print "........................................................"

#Get theory prediction for each analysis and print basic output
for analysis in listofanalyses:
    theorypredictions = theoryPredictionFor(analysis, smstoplist)
    if not theorypredictions:
        continue
    if parser.getboolean("stdout","printResults"):
        print "================================================================================"
        theorypredictions.printout() # again, check print function
    print "................................................................................"

    # Create a list of results, to determine the best result
    for theoryprediction in theorypredictions:
        res = ioObjects.ExptResults(theoryprediction.analysis.label.split(":")[0], theoryprediction.analysis.label.split(":")[1], rmvunit(theoryprediction.analysis.sqrts,"TeV"), theoryprediction.getmaxCondition(), rmvunit(theoryprediction.value[0].value,"fb"), rmvunit(theoryprediction.analysis.getUpperLimitFor(theoryprediction.mass),"fb"), smsResults.getConstraints(theoryprediction.analysis.label.split(":")[0], theoryprediction.analysis.label.split(":")[1]), rmvunit(theoryprediction.mass[0][0],"GeV"), rmvunit(theoryprediction.mass[0][-1],"GeV"))
        results.addResult(res)

# Find the best result, best result is defined as maximum 
results.findBest()

# If there is no best result, this means that there are no matching experimental results for the point.
# Decomposition status has following flags:
#-1: "#could not run the decomposition",
#-3: "#no cross sections above sigmacut found",
#-2: "#bad input slha, did not run decomposition",
#0: "#no matching experimental results",
#1: "#decomposition was successful". 

if not results.bestresult:
    decompstatus = 0  
    outputStatus = ioObjects.OutputStatus(decompstatus, slhastat, warnings)
    outputStatus.printout("file", args.outputfile)
    slhaStatus.printout("file", args.outputfile)
else:
    decompstatus = 1 
    outputStatus = ioObjects.OutputStatus(decompstatus, slhastat, warnings)
    outputStatus.printout("file", args.outputfile)
    slhaStatus.printout("file", args.outputfile)
    results.printout("file", args.outputfile)

missingtopos = ioObjects.MissingTopoList(sqrts)

missingtopos.findMissingTopos(smstoplist, listofanalyses, sigmacut, addunit(parser.getfloat("parameters","minmassgap"),"GeV"))

missingtopos.printout("file", args.outputfile)
