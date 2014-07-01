#!/usr/bin/env python

import sys, os, argparse, logging
import setPath
from tools.physicsUnits import rmvunit, addunit
from tools import slhaChecks, ioObjects, xsecComputer
from experiment import smsHelpers, smsAnalysisFactory, smsResults
from theory import slhaDecomposer
from smodels.theory.theoryPrediction import theoryPredictionFor

log = logging.getLogger(__name__)

#get name of input slha file (and parameter file)
argparser = argparse.ArgumentParser()
argparser.add_argument('-f', '--filename', help = 'filename of input slha', required=True)
argparser.add_argument('-p', '--parameterfile', help = 'filename of parameter file')
args = argparser.parse_args() # pylint: disable-msg=C0103

slhafile = args.filename #get input filename
#ioPar = ioObjects.InputParameters()
#if not ioPar.setFromFile(args.parameterfile): #read parameters from input file
#    log.error("Could not read %s" %args.parameterfile)
#    sys.exit()

try:
    if not args.parameterfile: import ioPar
    else: exec("import %s as ioPar" %args.parameterfile.split(".")[0]) #use parameter file given as argument
except: 
    log.error("Could not read parameter file")
    sys.exit()

# if all anas, topos are selected, set parameter to None
if not hasattr(ioPar,"analyses") or ioPar.analyses == all: ioPar.analyses = None
if not hasattr(ioPar,"topologies") or ioPar.topologies == all: ioPar.topologies = None

if os.path.exists("summary.txt"): #remove old output file
    log.warning("Remove old output file in summary.txt")
    os.remove("summary.txt")

#lists to store results
bestresult = []
outputarray = []

#get status, warnings from slhaChecks
slhaStatus = slhaChecks.SlhaStatus(slhafile, sigmacut=ioPar.sigmacut, maxDisplacement=.001, checkXsec=not ioPar.addMissingXsecs)
slhastat, warnings = slhaStatus.status

if slhastat == -1 or slhastat == -3:
    status = -2
    outputStatus = ioObjects.OutputStatus(status, slhastat, warnings)
    outputStatus.printout("file","summary.txt") #FIXME wie in ein bestimmtes file anhaengen?
    sys.exit()

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


#for now define some input parameters form xsection calcultaion, maybe later put them in parameter.in also?
sqrts = addunit(8,"TeV")
maxOrder = 2 # comput at NLL FIXME is not used at the moment! how to handle all different input options?

if computeXsecs:
    #first compute at LO
    xsecs = xsecComputer.computeXSec(sqrts, 0, int(ioPar.nevts), slhafile)
    comment = "Nevts: " + str(ioPar.nevts)
    xsecComputer.addXSecToFile(xsecs, slhafile, comment)

if ioPar.addnlo:
    xsecs_nlo = xsecComputer.computeXSec(sqrts, 1, ioPar.nevts, slhafile, loFromSlha=True)
    xsecComputer.addXSecToFile(xsecs_nlo, slhafile) #FIXME should there be a comment

if ioPar.addnll:
    xsecs_nll = xsecComputer.computeXSec(sqrts, 2, ioPar.nevts, slhafile, loFromSlha=True)
    xsecComputer.addXSecToFile(xsecs_nll, slhafile) #FIXME also: comment?

#set database adress
smsHelpers.base = "/home/laa/smodels-database" #FIXME check if this is right syntax in develop

#decomposition
sigmacut = addunit(ioPar.sigmacut,"fb")
try:
    smstoplist = slhaDecomposer.decompose(slhafile, sigmacut, doCompress=ioPar.doCompress,doInvisible=ioPar.doInvisible, minmassgap=addunit(ioPar.minmassgap,"GeV"))
except:
    status = -1
    outputStatus = ioObjects.OutputStatus(status, slhastat, warnings)
    outputStatus.printout("file","summary.txt")
    sys.exit()

if not smstoplist:
    status = -3
    outputStatus = ioObjects.OutputStatus(status, slhastat, warnings)
    outputStatus.printout("file","summary.txt")
    sys.exit()

# Print decomposition summary
if ioPar.printThEl: smstoplist.printout() #not my print function, check what this is doing, FIXME add function for printing GTops?

# Load analyses
listofanalyses = smsAnalysisFactory.load()

results = ioObjects.ResultList(bestresultonly = not ioPar.expandedSummary, describeTopo = ioPar.describeTopo)

#Get theory prediction for each analysis and print basic output
for analysis in listofanalyses:
    theorypredictions = theoryPredictionFor(analysis, smstoplist)
    if not theorypredictions:
        continue
    theorypredictions.printout() # again, check print function

    for theoryprediction in theorypredictions:
        res = ioObjects.ExptResults(theoryprediction.analysis.label.split(":")[0], theoryprediction.analysis.label.split(":")[1], rmvunit(theoryprediction.analysis.sqrts,"TeV"), theoryprediction.getmaxCondition(), rmvunit(theoryprediction.value[0].value,"fb"), rmvunit(theoryprediction.analysis.getUpperLimitFor(theoryprediction.mass),"fb"),smsResults.getConstraints(theoryprediction.analysis.label.split(":")[0], theoryprediction.analysis.label.split(":")[1]),rmvunit(theoryprediction.mass[0][0],"GeV"),rmvunit(theoryprediction.mass[0][-1],"GeV"))
        results.addResult(res)

results.findBest()

if not results.bestresult:
    status = 0
    outputStatus = ioObjects.OutputStatus(status, slhastat, warnings)
    outputStatus.printout("file","summary.txt")
else:
    status = 1
    outputStatus = ioObjects.OutputStatus(status, slhastat, warnings)
    outputStatus.printout("file","summary.txt")
    results.printout("file","summary.txt")

mt = ioObjects.MissingTopoList()

mt.findMissingTopos(smstoplist, listofanalyses, sigmacut, addunit(ioPar.minmassgap,"GeV"))

mt.printout("file","summary.txt")
