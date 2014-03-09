#!/usr/bin/python

import sys, glob, os
from theory import SLHADecomposer, CrossSection, SLHATools
from tools.PhysicsUnits import addunit
from experiment import SMSAnalysisFactory, LimitGetter
from numpy import arange

#ListOfAnalyses = SMSAnalysisFactory.load(sys.argv[1],sys.argv[2])
ListOfAnalyses = SMSAnalysisFactory.load('SUS13006','TChiChipmSlepStau')
mass = [[addunit(7.50E+02,'GeV'), addunit(7.36E+02,'GeV'), addunit(3.75E+02,'GeV')], [addunit(7.50E+02,'GeV'), addunit(7.36E+02,'GeV'), addunit(3.75E+02,'GeV')]]

for x in arange(0.05,0.95,0.01):
  mI = x*(mass[0][0]-mass[0][-1]) + mass[0][-1]
  mass = [[mass[0][0],mI,mass[0][2]],[mass[0][0],mI,mass[0][2]]]
  print x,LimitGetter.GetPlotLimit(mass,ListOfAnalyses[0],complain=False).asNumber()


sys.exit()
