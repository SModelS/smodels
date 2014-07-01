#!/usr/bin/env python

import sys, os, argparse
import setPath
from tools.physicsUnits import rmvunit, addunit
from tools import slhaChecks, ioObjects, xsecComputer
from experiment import smsHelpers, smsAnalysisFactory
from theory import slhaDecomposer
from smodels.theory.theoryPrediction import theoryPredictionFor

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

if not args.parameterfile: import ioPar
else: exec("import %s as ioPar" %args.parameterfile.split(".")[0]) #use parameter file given as argument

# if all anas, topos are selected, set parameter to None
if ioPar.analyses == all: ioPar.analyses = None
if ioPar.topologies == all: ioPar.topologies = None


if os.path.exists("summary.txt"): #remove old output file
    log.warning("Remove old output file in summary.txt")
    os.remove("summary.txt")

#lists to store results
bestresult = []
outputarray = []

#get status, warnings from slhaChecks
slhaStatus = slhaChecks.SlhaStatus(slhafile, sigmacut=ioPar.sigmacut, maxDisplacement=.001, checkXsec=not ioPar.addMissingXsecs)
slhastat, warnings = slhaStatus.status

if slhastat == -1:
    status = -2
    outputStatus = ioObjects.OutputStatus(status, slhastat, warnings)
    outputStatus.printout() #FIXME wie in ein bestimmtes file anhaengen?
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

    if ioPar.addnll:
        xsecs_nll = xsecComputer.computeXSec(sqrts, 2, ioPar.nevts, slhafile, loFromSlha=True)
        xsecComputer.addXSecToFile(xsecs_nll, slhafile, comment)

#set database adress
smsHelpers.base = "/home/laa/smodels-database" #FIXME check if this is right syntax in develop

#decomposition
smstoplist = slhaDecomposer.decompose(slhafile, ioPar.sigmacut, doCompress=ioPar.doCompress,doInvisible=ioPar.doInvisible, minmassgap=addunit(ioPar.minmassgap,"GeV"))

# Print decomposition summary
if ioPar.printThEl: smstoplist.printout() #not my print function, check what this is doing, FIXME add function for printing GTops?

# Load analyses
listofanalyses = smsAnalysisFactory.load()


#Get theory prediction for each analysis and print basic output
for analysis in listofanalyses:
    theorypredictions = theoryPredictionFor(analysis, smstoplist)
    if not theorypredictions:
        continue
    theorypredictions.printout() # again, check print function



