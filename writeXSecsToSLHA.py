#!/usr/bin/python

from tools.PhysicsUnits import addunit
from tools import xsecComputer


#Compute cross-sections using pythia and NLLfast (if needed) and add cross-section blocks to SLHA file:
slhafile = 'slha/andrePT4bb.slha'
sqrts = addunit(8.,'TeV')
order = 2  #maximum cross-section order available (2 = up to NLL+NLO, 1 = up to NLO, 0 = LO) 
nevts = 10000  #Number of events generated with pythia
xsecComputer.addXSecToFile(sqrts,order,nevts,slhafile)
sqrts = addunit(7.,'TeV')
xsecComputer.addXSecToFile(sqrts,order,nevts,slhafile)


