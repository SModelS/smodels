#!/usr/bin/env python

import sys, os, argparse, logging
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
argparser.add_argument('-f', '--filename', help = 'filename of input slha, necessary input', required = True)
argparser.add_argument('-p', '--parameterfile', help = 'filename of parameter file, optional input, default file = ./parameters.in', default = 'parameters.in')
argparser.add_argument('-o', '--outputfile', help = 'filename of output file, optional input, default outputfile = ./summary.txt', default = 'summary.txt')
args = argparser.parse_args() # pylint: disable-msg=C0103

slhafile = args.filename #get input filename
ioPar = ioObjects.InputParameters()

if not ioPar.setFromFile(args.parameterfile): #read parameters from input file
    log.error("Could not read %s" %args.parameterfile)
    sys.exit()

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

slhaStatus = slhaChecks.SlhaStatus(slhafile, sigmacut = ioPar.sigmacut, maxDisplacement = .001, checkXsec = not ioPar.addMissingXsecs, massgap = ioPar.minmassgap, maxcond = ioPar.maxcond)
slhastat, warnings = slhaStatus.status


if slhastat == -1 or slhastat == -3:
    status = -2
    outputStatus = ioObjects.OutputStatus(status, slhastat, warnings)
    outputStatus.printout("file",args.outputfile)
    slhaStatus.printout("file",args.outputfile)
    sys.exit()

# << there should be some automated mechanism to check this thing?>>
writeXsecs = True #set this true by default, switch to false only in case of lhe decomposition

#set cross section computation bools according to input parameters, slha
if not ioPar.doSLHAdec:
    computeXsecs = True
    writeXsecs = None
else:
    if ioPar.addMissingXsecs:
        computeXsecs = True
    elif not slhaStatus.xsec:
        log.warning("Input file does not contain cross sections, set computeXsecs = True")
        warnings = warnings + "#Cross sections computed by SModelS\n"
        computeXsecs = True
    else: computeXsecs = None

# sqrts from parameter input file
sqrts = addunit(ioPar.sqrts,"TeV")

if computeXsecs:
    #first compute at LO
    xsecs = xsecComputer.computeXSec(sqrts, 0, int(ioPar.nevts), slhafile)
    comment = "Nevts: " + str(ioPar.nevts)
    xsecComputer.addXSecToFile(xsecs, slhafile, comment)

if ioPar.addnlo:
    xsecs_nlo = xsecComputer.computeXSec(sqrts, 1, ioPar.nevts, slhafile, loFromSlha=True)
    xsecComputer.addXSecToFile(xsecs_nlo, slhafile) #FIXME should there be a comment << YES >>

if ioPar.addnll:
    xsecs_nll = xsecComputer.computeXSec(sqrts, 2, ioPar.nevts, slhafile, loFromSlha=True)
    xsecComputer.addXSecToFile(xsecs_nll, slhafile) #FIXME also: comment? << YES >>

#decomposition
#sigmacut = minimum value of cross-section for an element to be considered eligible for decomposition. Too small sigmacut leads to too large deocmposition time. 
sigmacut = addunit(ioPar.sigmacut,"fb")

try:
    # Decompose input SLHA file, store the output elements in smstoplist
    smstoplist = slhaDecomposer.decompose(slhafile, sigmacut, doCompress=ioPar.doCompress, doInvisible=ioPar.doInvisible, minmassgap=addunit(ioPar.minmassgap,"GeV"))
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

smstoplist.printout()
# This is my porposed format for element tabel

for (i,topo) in enumerate(smstoplist):
    print '\n'
    print "====================================================================="
    print "====================================================================="
    print "A new global topoloy starts here" 
    print "====================================================================="
    print "====================================================================="
    for j, el in enumerate(topo.elementList):
        print "........................................................................."
        print "........................................................................."
        #el.printout()
        print "Particles in topology:", el.getParticles()
        print 'The element masses are'
        for item in range(len(el.getMasses())):
            print "Masses branch %i:" %item, el.getMasses()[item]
        print "The element weights are:"
        for k in el.weight.getDictionary():
            print "Sqrts:", k, "\t Weights:", el.weight.getDictionary()[k]
    print "====================================================================="
    print "====================================================================="
    print "The list ends here" 
    print "====================================================================="
    print "====================================================================="

# Set database address
smsHelpers.base = ioPar.database

# Load analyses
listofanalyses = smsAnalysisFactory.load()

#<< What is this?>>
results = ioObjects.ResultList(bestresultonly = not ioPar.expandedSummary, describeTopo = ioPar.describeTopo)

constrainedElements = []

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
        if el != ref_el:
            print "\t Element is:", el
        ref_el = el
        for item in range(len(el.getMasses())):
            print "\t Element masses in branch %i are:" %item, el.getMasses()[item]
        print "\t The element weights are:"
        for k in el.weight.getDictionary():
            print "\t Sqrts:", k, "\t Weights:", el.weight.getDictionary()[k]
        print "\t ........................................................"
    sys.exit(10)

sys.exit(10)

#Get theory prediction for each analysis and print basic output
for analysis in listofanalyses:
    theorypredictions = theoryPredictionFor(analysis, smstoplist)
    if not theorypredictions:
        continue
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

missingtopos.findMissingTopos(smstoplist, listofanalyses, sigmacut, addunit(ioPar.minmassgap,"GeV"))

missingtopos.printout("file", args.outputfile)
